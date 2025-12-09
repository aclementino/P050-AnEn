#ifndef PROCESS_NETCDF
#define PROCESS_NETCDF

#include "structs.h"
#include "kdtree.h"

// =============================================================================
// MACROS PARA PROCESSAMENTO DE DADOS
// =============================================================================

#define MONACHE(TYPE)                                                    \
    case TYPE:                                                           \
    {                                                                    \
        for (int x = 0; x < ds->win_size; x++)                           \
        {                                                                \
            double diff =                                                \
                ((TYPE_VAR_##TYPE *)var->data)[(forecast - ds->k) + x] - \
                ((TYPE_VAR_##TYPE *)var->data)[(analog - ds->k) + x];    \
            sum += diff * diff;                                          \
        }                                                                \
    };                                                                   \
    break;

#define MONACHESW(TYPE)                                                            \
    case TYPE:                                                                     \
    {                                                                              \
        for (int x = 0; x < ds->win_size; x++)                                     \
        {                                                                          \
            double diff =                                                          \
                ((TYPE_VAR_##TYPE *)file[f].var[i].data)[(forecast - ds->k) + x] - \
                ((TYPE_VAR_##TYPE *)file[f].var[i].data)[(analog - ds->k) + x];    \
            sum += diff * diff;                                                    \
        }                                                                          \
    };                                                                             \
    break;

#define ALLOCATE_MEMORY_REC_DATA(TYPE, LENGTH)                                          \
    case TYPE:                                                                          \
        predicted_file->var[n].created_data = malloc(LENGTH * sizeof(TYPE_VAR_##TYPE)); \
        break;

#define REC_DATA(TYPE)                                                                                       \
    case TYPE:                                                                                               \
    {                                                                                                        \
        for (int i = 0; i < ds->num_Na; i++)                                                                 \
        {                                                                                                    \
            if (closest[i].window_index < 0 || closest[i].window_index >= file->dim->len)                    \
            {                                                                                                \
                continue;                                                                                    \
            }                                                                                                \
            if (isnan((double)((TYPE_VAR_##TYPE *)file->var[n].data)[closest[i].window_index]))              \
                continue;                                                                                    \
            else                                                                                             \
            {                                                                                                \
                if (!isnan((double)((TYPE_VAR_##TYPE *)file->var[n].data)[closest[i].window_index]))         \
                {                                                                                            \
                    sum_values += (double)(((TYPE_VAR_##TYPE *)file->var[n].data)[closest[i].window_index]); \
                    count++;                                                                                 \
                }                                                                                            \
            }                                                                                                \
        }                                                                                                    \
        ((TYPE_VAR_##TYPE *)file->var[n].created_data)[f_count] = (TYPE_VAR_##TYPE)(sum_values / count);     \
    }                                                                                                        \
    break;

#define CALC_RMSE(TYPE)                                                                        \
    case TYPE:                                                                                 \
        for (int i = ds->start_prediction; i <= ds->end_prediction; i++)                       \
        {                                                                                      \
            if (isnan((double)((TYPE_VAR_##TYPE *)file->var[n].data)[i]))                      \
            {                                                                                  \
                ii++;                                                                          \
                continue;                                                                      \
            }                                                                                  \
            else                                                                               \
            {                                                                                  \
                if (isnan((double)((TYPE_VAR_##TYPE *)file->var[n].created_data)[ii]))         \
                {                                                                              \
                    ii++;                                                                      \
                    continue;                                                                  \
                }                                                                              \
                else                                                                           \
                {                                                                              \
                    double error = (double)((TYPE_VAR_##TYPE *)file->var[n].data)[i] -         \
                                   (double)((TYPE_VAR_##TYPE *)file->var[n].created_data)[ii]; \
                    sum_error += error * error;                                                \
                    ii++;                                                                      \
                    count++;                                                                   \
                }                                                                              \
            }                                                                                  \
        }                                                                                      \
        break;

// =============================================================================
// ESTRUTURAS PARA ALGORITMOS OTIMIZADOS
// =============================================================================

/**
 * @brief Estrutura para dados pré-filtrados (ANEN)
 *
 * Armazena arrays de índices válidos para forecasts e analogs,
 * eliminando a necessidade de validação durante o processamento.
 */
typedef struct
{
    int *valid_forecasts;    // Array de índices de forecasts válidos
    int *valid_analogs;      // Array de índices de analogs válidos
    int num_valid_forecasts; // Quantidade de forecasts válidos
    int num_valid_analogs;   // Quantidade de analogs válidos
} PreFilteredData;

/**
 * @brief Dados compartilhados entre threads para algoritmo ANEN
 *
 * Estrutura thread-safe que permite múltiplas threads processarem
 * diferentes forecasts usando os mesmos dados pré-filtrados.
 */
typedef struct
{
    NetCDF *predicted_file;         // Arquivo de dados preditos
    NetCDF *predictor_file;         // Arquivo de dados preditores
    DataSegment *ds;                // Configurações do algoritmo
    int n;                          // Índice da variável sendo processada
    PreFilteredData *filtered_data; // Dados pré-filtrados (read-only)
} ANENSharedData;

/**
 * @brief Dados específicos de cada worker thread para ANEN
 *
 * Cada thread recebe uma faixa de forecasts para processar
 * e mantém suas próprias estatísticas de performance.
 */
typedef struct
{
    ANENSharedData *shared; // Dados compartilhados
    int thread_id;          // ID da thread (0 a num_threads-1)
    int start_forecast_idx; // Índice inicial no array valid_forecasts
    int end_forecast_idx;   // Índice final no array valid_forecasts
    int processed_count;    // Contador local de forecasts processados
    double processing_time; // Tempo de processamento desta thread
} ANENWorkerData;

// =============================================================================
// ESTRUTURAS PARA ALGORITMO KD-ANEN (KD-TREE + ANALOG ENSEMBLE)
// =============================================================================

/**
 * @brief Dados compartilhados entre threads para algoritmo KD-ANEN
 *
 * Similar ao ANEN mas usa KD-Tree pré-construída para busca eficiente.
 * A árvore é construída uma vez e usada por todas as threads.
 */
typedef struct
{
    NetCDF *predicted_file;  // Arquivo de dados preditos (escrita thread-safe)
    NetCDF *predictor_file;  // Arquivo de dados preditores (read-only)
    DataSegment *ds;         // Configurações do algoritmo (read-only)
    int n;                   // Índice da variável sendo processada
    KDTree *root;            // Raiz da KD-Tree (read-only, thread-safe)
    int *valid_forecasts;    // Array de forecasts válidos (read-only)
    int num_valid_forecasts; // Quantidade de forecasts válidos
} KDANENSharedData;

/**
 * @brief Dados específicos de cada worker thread para KD-ANEN
 *
 * Cada thread processa uma faixa de forecasts usando a KD-Tree
 * compartilhada para busca logarítmica de vizinhos.
 */
typedef struct
{
    KDANENSharedData *shared; // Dados compartilhados
    int thread_id;            // ID da thread (0 a num_threads-1)
    int start_forecast_idx;   // Índice inicial no array valid_forecasts
    int end_forecast_idx;     // Índice final no array valid_forecasts
    int processed_count;      // Contador local de forecasts processados
    double processing_time;   // Tempo de processamento desta thread
} KDANENWorkerData;

// =============================================================================
// FUNÇÕES PRINCIPAIS DOS ALGORITMOS
// =============================================================================

/**
 * @brief Algoritmo ANEN Paralelo - Analog Ensemble com pré-filtro
 *
 * Implementa busca exaustiva otimizada com pré-filtro de dados válidos.
 * Elimina validações durante o processamento principal através de
 * pré-computação sequencial dos índices válidos.
 *
 * Características:
 * - Complexidade: O(n×m/t) onde n=forecasts, m=analogs, t=threads
 * - Speedup típico: 3-4x com 4-8 threads
 * - Overhead: ~5% para pré-filtro
 * - Ideal para: datasets médios/grandes com boa distribuição de dados válidos
 *
 * @param file Array de arquivos NetCDF [predicted, predictor]
 * @param ds Configurações do algoritmo (períodos, janelas, threads)
 */
void anen_dependent_parallel(NetCDF *file, DataSegment *ds);

/**
 * @brief Algoritmo KD-ANEN Paralelo - KD-Tree + Analog Ensemble
 *
 * Implementa busca inteligente usando KD-Tree para reduzir complexidade
 * de O(n×m) para O(n×log m). Ideal para datasets grandes onde o custo
 * de construção da árvore é compensado pela eficiência da busca.
 *
 * Características:
 * - Complexidade: O(n×log m/t) + O(m×log m) para construção
 * - Speedup típico: 8-25x com 4-8 threads
 * - Overhead: ~10-20% para construção da árvore
 * - Ideal para: datasets grandes com muitos analogs (>500k)
 *
 * @param file Array de arquivos NetCDF [predicted, predictor]
 * @param ds Configurações do algoritmo (períodos, janelas, threads)
 */
void kdanen_independent_parallel(NetCDF *file, DataSegment *ds);

// =============================================================================
// FUNÇÕES AUXILIARES PARA ANEN
// =============================================================================

/**
 * @brief Inicializa estrutura de dados pré-filtrados
 *
 * Executa validação sequencial de todos os forecasts e analogs,
 * armazenando apenas os índices válidos em arrays compactos.
 *
 * @param filtered Estrutura a ser inicializada
 * @param var Variável a ser validada
 * @param ds Configurações do algoritmo
 * @param data_length Tamanho total dos dados
 */
void init_prefiltered_data(PreFilteredData *filtered, Variable *var,
                           DataSegment *ds, int data_length);

/**
 * @brief Libera memória dos dados pré-filtrados
 *
 * @param filtered Estrutura a ser liberada
 */
void free_prefiltered_data(PreFilteredData *filtered);

/**
 * @brief Worker thread para processamento ANEN
 *
 * Função executada por cada thread no algoritmo ANEN paralelo.
 * Processa uma faixa de forecasts usando dados pré-filtrados.
 *
 * @param arg Ponteiro para ANENWorkerData
 * @return NULL
 */
void *anen_parallel_worker(void *arg);

// =============================================================================
// FUNÇÕES AUXILIARES GERAIS
// =============================================================================

/**
 * @brief Validação simples de janela de dados
 *
 * Verifica se todos os valores em uma janela centrada na posição
 * especificada são válidos (não NaN).
 *
 * @param var Variável a ser validada
 * @param position Posição central da janela
 * @param k Raio da janela (tamanho = 2k+1)
 * @param win_size Tamanho total da janela
 * @param data_length Tamanho total dos dados
 * @return true se janela é válida, false caso contrário
 */
bool validate_window_simple(Variable *var, int position, int k,
                            int win_size, int data_length);

/**
 * @brief Função de processamento genérica
 *
 * Interface comum para todos os algoritmos de processamento.
 *
 * @param file Array de arquivos NetCDF
 * @param ds Configurações do algoritmo
 * @param func Ponteiro para função específica do algoritmo
 */
void processing_data(NetCDF *file, DataSegment *ds, process_func func);

/**
 * @brief Cálculo da métrica de distância Monache
 *
 * Calcula distância euclidiana entre duas janelas de dados.
 *
 * @param var Variável sendo processada
 * @param ds Configurações do algoritmo
 * @param forecast Posição do forecast
 * @param analog Posição do analog
 * @param i Índice da variável
 * @return Distância calculada
 */
double monache_metric(Variable *var, DataSegment *ds, int forecast, int analog, int i);

/**
 * @brief Cálculo da métrica Monache para múltiplas séries
 *
 * Versão estendida que considera múltiplas séries temporais.
 *
 * @param file Array de arquivos NetCDF
 * @param ds Configurações do algoritmo
 * @param forecast Posição do forecast
 * @param analog Posição do analog
 * @param i Índice da variável
 * @return Distância calculada
 */
double monache_metric_super_window(NetCDF *file, DataSegment *ds,
                                   int forecast, int analog, int i);

/**
 * @brief Comparador para ordenação de pontos próximos
 *
 * Função de comparação para qsort, ordena por distância (maior primeiro).
 *
 * @param x Primeiro elemento
 * @param y Segundo elemento
 * @return Resultado da comparação
 */
int compare_closest_point_ord_const(const void *x, const void *y);

/**
 * @brief Reconstrói dados usando vizinhos mais próximos
 *
 * Versão corrigida que mapeia corretamente os índices entre
 * dados originais e dados reconstruídos.
 *
 * @param file Arquivo NetCDF
 * @param ds Configurações do algoritmo
 * @param closest Array de pontos mais próximos
 * @param forecast_position Posição no array de dados reconstruídos
 * @param n Índice da variável
 * @param num_found Número de vizinhos encontrados
 */
void recreate_data_fixed(NetCDF *file, DataSegment *ds, ClosestPoint *closest,
                         int forecast_position, int n, int num_found);

/**
 * @brief Calcula RMSE com mapeamento correto de índices
 *
 * Versão corrigida que compara corretamente dados originais
 * com dados reconstruídos.
 *
 * @param file Arquivo NetCDF
 * @param ds Configurações do algoritmo
 * @param n Índice da variável
 */
void calculate_rmse_fixed(NetCDF *file, DataSegment *ds, int n);

/**
 * @brief Aloca array de pontos próximos com inicialização segura
 *
 * Previne uso de dados não inicializados que podem causar
 * comportamento indefinido.
 *
 * @param num_points Número de pontos a alocar
 * @return Ponteiro para array alocado ou NULL se erro
 */
ClosestPoint *allocate_closest_points_safe(int num_points);

/**
 * @brief Valida processo completo de reconstrução
 *
 * Verifica se a reconstrução foi bem-sucedida e se há
 * dados válidos para cálculo de métricas.
 *
 * @param file Arquivo NetCDF
 * @param ds Configurações do algoritmo
 * @param n Índice da variável
 * @return true se processo foi bem-sucedido
 */
bool validate_reconstruction_process(NetCDF *file, DataSegment *ds, int n);

// =============================================================================
// FUNÇÕES LEGACY (COMPATIBILIDADE)
// =============================================================================

/**
 * @brief Algoritmo exaustivo independente (legacy)
 *
 * Implementação original para compatibilidade.
 * Recomenda-se usar anen_dependent_parallel para novos projetos.
 */
void exhaustive_processing_independent(NetCDF *file, DataSegment *ds);

/**
 * @brief Algoritmo exaustivo dependente (legacy)
 *
 * Implementação original para compatibilidade.
 * Recomenda-se usar anen_dependent_parallel para novos projetos.
 */
void exhaustive_processing_dependent(NetCDF *file, DataSegment *ds);

/**
 * @brief Versão paralela corrigida do algoritmo exaustivo (legacy)
 *
 * Implementação com correções de bugs mas sem otimizações avançadas.
 * Recomenda-se usar anen_dependent_parallel para melhor performance.
 */
void exhaustive_processing_dependent_fixed_parallel(NetCDF *file, DataSegment *ds);

/**
 * @brief Versão paralela corrigida do algoritmo parcial (legacy)
 *
 * Implementação com KD-Tree mas interface legacy.
 * Recomenda-se usar kdanen_independent_parallel para melhor interface.
 */
void partial_processing_independent_fixed_parallel(NetCDF *file, DataSegment *ds);

/* Init - Multi-Series KD-Tree structures */

/**
 * @brief Estrutura estendida para KD-Tree com múltiplas séries
 */
typedef struct KDTreeMultiSeries
{
    unsigned int window_id;
    int total_dimensions;
    struct KDTreeMultiSeries *left;
    struct KDTreeMultiSeries *right;
} KDTreeMultiSeries;

/**
 * @brief Pool de nós para KD-Tree de múltiplas séries
 */
typedef struct
{
    KDTreeMultiSeries nodes[NODE_POOL_SIZE];
    int next_available;
} MultiSeriesNodePool;

/**
 * @brief Dados compartilhados para threads no algoritmo KD-ANEN dependent
 */
typedef struct
{
    NetCDF *predicted_file;
    NetCDF *predictor_file;
    DataSegment *ds;
    int n;
    KDTreeMultiSeries *root;
    int *valid_forecasts;
    int num_valid_forecasts;
    int total_dimensions;
} KDANENDependentSharedData;

/**
 * @brief Dados específicos de cada worker thread para KD-ANEN dependent
 */
typedef struct
{
    KDANENDependentSharedData *shared;
    int thread_id;
    int start_forecast_idx;
    int end_forecast_idx;
    int processed_count;
    double reconstruct_time;
    double processing_time;
} KDANENDependentWorkerData;

/* End - Multi-Series KD-Tree structures */

// =============================================================================
// FUNÇÕES PARA KD-ANEN DEPENDENT (MÚLTIPLAS SÉRIES)
// =============================================================================

/**
 * @brief Algoritmo KD-ANEN Dependent Paralelo - KD-Tree para múltiplas séries
 */
void kdanen_dependent_parallel(NetCDF *file, DataSegment *ds);

/**
 * @brief Cria pool de nós para KD-Tree de múltiplas séries
 */
MultiSeriesNodePool *create_multiseries_node_pool(void);

/**
 * @brief Aloca nó da pool para múltiplas séries
 */
KDTreeMultiSeries *allocate_multiseries_node_from_pool(MultiSeriesNodePool *pool, int window_id, int total_dims);

/**
 * @brief Obtém valor de uma dimensão específica da super janela
 */
double get_multiseries_value(NetCDF *file, DataSegment *ds, int window_id, int dimension, int var_idx);

/**
 * @brief Constrói KD-Tree balanceada para múltiplas séries
 */
KDTreeMultiSeries *build_multiseries_balanced_kdtree(int *window_ids, int n, NetCDF *file, DataSegment *ds,
                                                     int depth, MultiSeriesNodePool *pool, int var_idx);

/**
 * @brief Calcula distância quadrática entre super janelas (múltiplas séries)
 */
double squared_distance_multiseries(KDTreeMultiSeries *root, NetCDF *file, DataSegment *ds,
                                    int target_id, int node_id, int var_idx);

/**
 * @brief Busca pontos mais próximos na KD-Tree de múltiplas séries
 */
void search_multiseries_closest_points(KDTreeMultiSeries *root, NetCDF *file, DataSegment *ds,
                                       ClosestPoint *closest, int target_id, int depth,
                                       int var_idx, int *found);

/**
 * @brief Worker thread para processamento KD-ANEN dependent
 */
void *kdanen_dependent_parallel_worker(void *arg);

/**
 * @brief Estrutura para passar contexto via variável global thread-safe
 */
typedef struct
{
    NetCDF *file;
    DataSegment *ds;
    int axis;
    int var_idx;
} SortContextMulti;

// Variável global thread-local para contexto
static __thread SortContextMulti *current_sort_context = NULL;

/**
 * @brief Função de comparação para qsort padrão
 */
int compare_multiseries_points(const void *a, const void *b);

/**
 * @brief Ordena pontos por eixo para múltiplas séries (versão compatível)
 */
void sort_multiseries_points_by_axis(int *points, int n, NetCDF *file, DataSegment *ds, int axis, int var_idx);

// =============================================================================
// KD-ANEN DEPENDENT PARALLEL - VERSÃO ENTRELAÇADA (INTERLEAVED)
// =============================================================================

/**
 * @brief Obtém valor de uma dimensão específica - LAYOUT ENTRELAÇADO
 *
 * Layout entrelaçado organiza dimensões como:
 * Dim 0→S1P0, Dim 1→S2P0, Dim 2→S3P0, Dim 3→S1P1, Dim 4→S2P1, Dim 5→S3P1, ...
 */
double get_multiseries_value_interleaved(NetCDF *file, DataSegment *ds, int window_id, int dimension, int var_idx);

/**
 * @brief Função de comparação para qsort padrão - LAYOUT ENTRELAÇADO
 */
int compare_multiseries_standard_interleaved(const void *a, const void *b);

/**
 * @brief Ordena pontos por eixo para múltiplas séries - LAYOUT ENTRELAÇADO
 */
void sort_multiseries_points_by_axis_interleaved(int *points, int n, NetCDF *file, DataSegment *ds, int axis, int var_idx);

/**
 * @brief Constrói KD-Tree balanceada para múltiplas séries - LAYOUT ENTRELAÇADO
 */
KDTreeMultiSeries *build_multiseries_balanced_kdtree_interleaved(int *window_ids, int n, NetCDF *file, DataSegment *ds,
                                                                 int depth, MultiSeriesNodePool *pool, int var_idx);

/**
 * @brief Worker thread para processamento KD-ANEN dependent - LAYOUT ENTRELAÇADO
 */
void *kdanen_dependent_parallel_worker_interleaved(void *arg);

/**
 * @brief Algoritmo KD-ANEN Dependent Paralelo - VERSÃO ENTRELAÇADA
 *
 * Implementa busca eficiente em espaço multidimensional com layout entrelaçado
 * para melhor cache locality e poda mais eficiente.
 */
void kdanen_dependent_parallel_interleaved(NetCDF *file, DataSegment *ds);

/**
 * @brief Função de benchmark para comparar layouts
 *
 * Executa ambas as versões e compara performance mantendo correção.
 */
void benchmark_kdanen_layouts(NetCDF *file, DataSegment *ds);

#endif
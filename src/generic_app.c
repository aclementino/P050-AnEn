// #include "structs.h"
#include "kdtree.h"
#include "randw.h"
#include "preprocess.h"
#include "process.h"

// =============================================================================
// CONFIGURACOES DE PERIODO
// =============================================================================

// DATE_ISO - format "YYYY-MM-DDTHH:MM:00"
#define PREDICTION_INIT "2019-01-01T00:00:00"
#define PREDICTION_END "2019-12-31T23:54:00"
#define TRAINING_INIT(x) x
#define TRAINING_END(x) x

// =============================================================================
// VARIAVEIS GLOBAIS
// =============================================================================

// Start measuring time
struct timeval begin, end;
long seconds = 0;
long microseconds = 0;
double elapsed = 0;

KDTree *kdtree;

/**
 * @brief Função principal do programa
 *
 * Processa dados NetCDF usando algoritmos de Analog Ensemble.
 * Suporta diferentes períodos de treino e algoritmos otimizados.
 *
 * Argumentos:
 * argv[1] - Número de threads (1, 2, 4, 8, etc.)
 * argv[2] - Anos de treino (1, 2, 4, 8)
 * argv[3...] - Arquivos NetCDF (primeiro = predito, demais = preditores)
 *
 * Exemplo de uso:
 * ./programa 4 8 dados_preditos.nc dados_preditores.nc
 */
int main(int argc, char *argv[])
{
    char *T_INIT = NULL;
    char *T_END = NULL;

    // =============================================================================
    // CONFIGURAÇÃO DE PERÍODOS DE TREINO
    // =============================================================================

    switch (strtol(argv[2], NULL, 10))
    {
    case 1:
        T_INIT = "2018-01-01T00:00:00"; // 1 ano de treino
        T_END = "2018-12-31T00:00:00";
        break;
    case 2:
        T_INIT = "2017-01-01T00:00:00"; // 2 anos de treino
        T_END = "2018-12-31T00:00:00";
        break;
    case 4:
        T_INIT = "2015-01-01T00:00:00"; // 4 anos de treino
        T_END = "2018-12-31T00:00:00";
        break;
    case 8:
        T_INIT = "2011-01-01T00:00:00"; // 8 anos de treino  4003020
        T_END = "2018-12-31T00:00:00";
        break;
    default:
        fprintf(stderr, "Erro: Período de treino inválido. Use 1, 2, 4 ou 8 anos.\n");
        fprintf(stderr, "Uso: %s <threads> <anos_treino> <arquivo_predito> <arquivo_preditor>\n", argv[0]);
        return EXIT_FAILURE;
        break;
    }

    // =============================================================================
    // CABEÇALHO DE SAÍDA
    // =============================================================================

    // printf("=== ANALOG ENSEMBLE - VERSÃO OTIMIZADA ===\n");
    // printf("Período de predição: %s a %s\n", PREDICTION_INIT, PREDICTION_END);
    // printf("Período de treino: %s a %s (%s ano(os))\n", T_INIT, T_END, argv[2]);
    // printf("Threads: %s\n", argv[1]);
    // printf("Arquivos: %d\n", argc - 3);
    // printf("\n");

    // Cabeçalho CSV para resultados
    printf("n_files,n_threads,t_rdfiles,s_training,e_training,s_prediction,e_prediction,algorithm,t_process,rmse,t_total\n");

    // =============================================================================
    // CONFIGURAÇÃO DO ALGORITMO
    // =============================================================================

    DataSegment ds;

    // Parâmetros do algoritmo
    ds.k = 5;                                   // Raio da janela (janela = 2k+1 = 11 pontos)
    ds.win_size = (ds.k * 2) + 1;               // Tamanho da janela
    ds.win_size_interpolation = (ds.k * 2) - 1; // Tamanho máximo para interpolação
    ds.num_Na = 25;                             // Número de vizinhos mais próximos
    ds.num_thread = strtol(argv[1], NULL, 10);  // Número de threads
    ds.argc = argc - 3;                         // Número de arquivos
    ds.indice_generic = 0;                      // Índice genérico para processamento

    printf("%i,%i,", ds.argc, ds.num_thread);

    // =============================================================================
    // CARREGAMENTO DE DADOS
    // =============================================================================

    GET_START;
    NetCDF *file = create_struct(&ds, &argv[3]);
    GET_END;

    if (!file)
    {
        fprintf(stderr, "Erro: Falha ao carregar estruturas NetCDF.\n");
        return EXIT_FAILURE;
    }

    // =============================================================================
    // CONFIGURAÇÃO DE PERÍODOS
    // =============================================================================

    ds.start_training = binary_search(file, convert_time(TRAINING_INIT(T_INIT))) + ds.k;
    ds.end_training = binary_search(file, convert_time(TRAINING_END(T_END)));
    ds.start_prediction = binary_search(file, convert_time(PREDICTION_INIT));
    ds.end_prediction = binary_search(file, convert_time(PREDICTION_END)) - ds.k;

    printf("%i,%i,%i,%i,",
           ds.start_training,
           ds.end_training,
           ds.start_prediction,
           ds.end_prediction);

    // Validação dos períodos
    if (ds.start_training < 0 || ds.end_training < 0 ||
        ds.start_prediction < 0 || ds.end_prediction < 0)
    {
        fprintf(stderr, "Erro: Períodos inválidos encontrados.\n");
        deallocate_memory(file, ds.argc);
        return EXIT_FAILURE;
    }

    // printf("Período de treino: %d a %d (%d pontos)\n",
    //        ds.start_training, ds.end_training, ds.end_training - ds.start_training + 1);
    // printf("Período de predição: %d a %d (%d pontos)\n",
    //        ds.start_prediction, ds.end_prediction, ds.end_prediction - ds.start_prediction + 1);

    // =============================================================================
    // PRÉ-PROCESSAMENTO DOS DADOS
    // =============================================================================

    // printf("\n=== PRÉ-PROCESSAMENTO ===\n");

    // Análise de dados inválidos
    analyze_data(file, &ds, count_invalid_values);

    GET_START;
    // Interpolação de valores faltantes
    analyze_data(file, &ds, interpolation_values);
    GET_END;

    // Recontagem após interpolação
    analyze_data(file, &ds, count_invalid_values);

    // Contagem de janelas válidas
    analyze_data(file, &ds, count_valid_window);

    // Estatísticas de qualidade dos dados
        // for (int i = 0; i < ds.argc; i++)
        // {
        //     printf("Arquivo %d: %d janelas válidas (treino), %d janelas válidas (predição)\n",
        //            i, file[i].var[1].num_win_valid_training, file[i].var[1].num_win_valid_prediction);
        // }

    // =============================================================================
    // SELEÇÃO E EXECUÇÃO DO ALGORITMO
    // =============================================================================

    // printf("\n=== PROCESSAMENTO ===\n");

    GET_START;

    // processing_data(file, &ds, kdanen_independent_parallel);
    // processing_data(file, &ds, anen_dependent_parallel);
    processing_data(file, &ds, kdanen_dependent_parallel);
    // processing_data(file, &ds, kdanen_dependent_parallel_interleaved);

    GET_END;

    // =============================================================================
    // FINALIZAÇÃO
    // =============================================================================

    // printf("\n=== LIMPEZA E FINALIZAÇÃO ===\n");

    // Liberar memória
    deallocate_memory(file, ds.argc);
    file = NULL;

    // Finalizar NetCDF
    nc_finalize();

    // printf("Processamento concluído com sucesso.\n");

    return EXIT_SUCCESS;
}
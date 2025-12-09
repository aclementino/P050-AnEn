#ifndef PROCESS_KDTREE
#define PROCESS_KDTREE

#include "structs.h"

#define IFNAN_KDTREE(TYPE)                                                           \
    case TYPE:                                                                       \
    {                                                                                \
        if (!isnan((float)((TYPE_VAR_##TYPE *)predictor_file->var[n].data)[analog])) \
            root = insert_kdt_node(root, &predictor_file->var[n], ds, analog, 0);    \
    };                                                                               \
    break;

#define IF_AXIS(TYPE)                                                         \
    case TYPE:                                                                \
    {                                                                         \
        if (((TYPE_VAR_##TYPE *)var->data)[(window_id - ds->k) + axis] <      \
            ((TYPE_VAR_##TYPE *)var->data)[(root->window_id - ds->k) + axis]) \
            root->left =                                                      \
                insert_kdt_node(root->left, var, ds, window_id, depth + 1);   \
        else                                                                  \
            root->right =                                                     \
                insert_kdt_node(root->right, var, ds, window_id, depth + 1);  \
        return root;                                                          \
    };                                                                        \
    break;

#define MONACHE_KDTREE(TYPE)                                                             \
    case TYPE:                                                                           \
    {                                                                                    \
        for (int x = 0; x < ds->win_size; x++)                                           \
        {                                                                                \
            double diff =                                                                \
                (double)(((TYPE_VAR_##TYPE *)var->data)[(target_id - ds->k) + x] -       \
                         ((TYPE_VAR_##TYPE *)var->data)[(root->window_id - ds->k) + x]); \
            sum += diff * diff;                                                          \
        }                                                                                \
    };                                                                                   \
    break;

#define MONACHESW_KDTREE(TYPE)                                                           \
    case TYPE:                                                                           \
    {                                                                                    \
        for (int x = 0; x < ds->win_size; x++)                                           \
        {                                                                                \
            double diff =                                                                \
                ((TYPE_VAR_##TYPE *)file[f].var[i].data)[(target_id - ds->k) + x] -      \
                ((TYPE_VAR_##TYPE *)file[f].var[i].data)[(root->window_id - ds->k) + x]; \
            sum += diff * diff;                                                          \
        }                                                                                \
    };                                                                                   \
    break;

#define IF_STORED(TYPE)                                                              \
    case TYPE:                                                                       \
    {                                                                                \
        if (((TYPE_VAR_##TYPE *)var->data)[(target_id - ds->k) + axis] <             \
            ((TYPE_VAR_##TYPE *)var->data)[(root->window_id - ds->k) + axis])        \
        {                                                                            \
            next = root->left;                                                       \
            other = root->right;                                                     \
        }                                                                            \
        else                                                                         \
        {                                                                            \
            next = root->right;                                                      \
            other = root->left;                                                      \
        }                                                                            \
        search_closest_points(                                                       \
            next, var, ds, closest, target_id, (depth + 1), found);                  \
        if (*found < ds->num_Na ||                                                   \
            fabs(((TYPE_VAR_##TYPE *)var->data)[(target_id - ds->k) + axis] -        \
                 ((TYPE_VAR_##TYPE *)var->data)[(root->window_id - ds->k) + axis]) < \
                closest[0].distance)                                                 \
            search_closest_points(                                                   \
                other, var, ds, closest, target_id, (depth + 1), found);             \
    };                                                                               \
    break;

#define IF_STORED_SW(TYPE)                                                                   \
    case TYPE:                                                                               \
    {                                                                                        \
        if (((TYPE_VAR_##TYPE *)file->var[i].data)[(target_id - ds->k) + axis] <             \
            ((TYPE_VAR_##TYPE *)file->var[i].data)[(root->window_id - ds->k) + axis])        \
        {                                                                                    \
            next = root->left;                                                               \
            other = root->right;                                                             \
        }                                                                                    \
        else                                                                                 \
        {                                                                                    \
            next = root->right;                                                              \
            other = root->left;                                                              \
        }                                                                                    \
        search_closest_points_super_window(                                                  \
            next, file, ds, closest, target_id, depth + 1, i, found);                        \
        if (*found < ds->num_Na ||                                                           \
            fabs(((TYPE_VAR_##TYPE *)file->var[i].data)[(target_id - ds->k) + axis] -        \
                 ((TYPE_VAR_##TYPE *)file->var[i].data)[(root->window_id - ds->k) + axis]) < \
                closest[0].distance)                                                         \
            search_closest_points_super_window(                                              \
                other, file, ds, closest, target_id, depth + 1, i, found);                   \
    };                                                                                       \
    break;

// K-D Tree node structure
typedef struct KDTree
{
    unsigned int window_id;
    struct KDTree *left;
    struct KDTree *right;
} KDTree;

typedef struct
{
    NetCDF *predicted_file;
    NetCDF *predictor_file;
    DataSegment *ds;
    int n;                   // Índice da variável
    KDTree *root;            // Árvore KD compartilhada (somente leitura)
    int *valid_forecasts;    // Array de forecasts válidos
    int num_valid_forecasts; // Total de forecasts válidos
} ThreadData;

typedef struct
{
    ThreadData *shared;
    int thread_id;
    int start_idx;
    int end_idx;
} WorkerData;

// /* --- kdtree 3 ---
// Node pool for efficient memory management
#define NODE_POOL_SIZE 1000
typedef struct
{
    KDTree nodes[NODE_POOL_SIZE];
    int next_available;
} NodePool;
NodePool *create_node_pool();
KDTree *allocate_node_from_pool(NodePool *pool, int window_id);
void reset_node_pool(NodePool *pool);
void free_node_pool(NodePool *pool);
// Standard node creation function (fallback)
KDTree *create_kdt_node(int window_id);
// Improved comparison context
typedef struct
{
    Variable *var;
    int axis;
    int k; // Offset for window indexing
} SortContext;
// Optimized comparison function
int compare_wrapper(const void *a, const void *b, void *context);
// Optimized sorting function using qsort_r where available
void sort_points_by_axis(int *points, int n, Variable *var, int axis, int k);
// Optimized median finding using median-of-medians algorithm for better pivot selection
int select_median(int *arr, int n, Variable *var, int axis, int k);
// Improved balanced KD-tree builder
KDTree *build_balanced_kdtree(int *window_ids, int n, Variable *var, DataSegment *ds, int depth, NodePool *pool);
// Optimized KD-tree creation
KDTree *create_balanced_kdtree(Variable *var, DataSegment *ds);
// Optimized insert function that maintains better balance
KDTree *insert_kdt_node(KDTree *root, Variable *var, DataSegment *ds, int window_id, int depth);
// Optimized distance calculation (squared distance to avoid sqrt)
double squared_distance_kdtree(KDTree *root, Variable *var, DataSegment *ds, int window_id, int target_id);
// Calculate distance and apply sqrt only when needed
double monache_metric_kdtree(KDTree *root, Variable *var, DataSegment *ds, int window_id, int target_id);
double monache_metric_super_window_kdtree(KDTree *root, NetCDF *file, DataSegment *ds, int window_id, int target_id, int i);
// Optimized minimum distance between point and hyperrectangle (for pruning)
double min_distance_to_hyperrect(Variable *var, DataSegment *ds, int target_id, KDTree *node, int depth);
// Optimized closest points search with a priority queue
void search_closest_points(KDTree *root, Variable *var, DataSegment *ds, ClosestPoint *closest, int target_id, int depth, int *found);
// Fixed and optimized version of search_closest_points_super_window
void search_closest_points_super_window(KDTree *root, NetCDF *file, DataSegment *ds, ClosestPoint *closest, int target_id, int depth, int i, int *found);
// Helper function to collect window IDs from the tree
void collect_window_ids(KDTree *root, int *window_ids, int *index);
// Efficient tree deallocation
void deallocate_kdtree(KDTree *root);
// Improved rebalance function with node pool
KDTree *rebalance_kdtree(KDTree *root, Variable *var, DataSegment *ds);
// Function to check if the k-d tree is balanced
bool is_kdtree_balanced(KDTree *root);
// Function to get the height of a k-d tree
int get_tree_height(KDTree *node);
// Function to count the number of nodes in the tree
int count_nodes(KDTree *node);
// Optimized comparison function for qsort
int compare_near_point(const void *x, const void *y);
// Enhanced diagnostic function
void diagnose_tree_balance(KDTree *root);
// Function to visualize the tree structure
void visualize_kdtree(KDTree *root, int depth);
// Batch processing function for multiple queries
void batch_nearest_neighbors(KDTree *root, Variable *var, DataSegment *ds,
                             int *query_points, int num_queries,
                             ClosestPoint **results);
// A hybrid approach combining KD-tree and local brute force search
void hybrid_nearest_neighbor_search(KDTree *root, Variable *var, DataSegment *ds,
                                    int target_id, ClosestPoint *closest,
                                    int *candidate_ids, int *num_candidates);

int partition_around_value(int *arr, int n, int pivot_value, Variable *var, int axis, int k);
int compare_values(int window_a, int window_b, Variable *var, int axis, int k);
KDTree *build_optimized_balanced_kdtree(int *window_ids, int n, Variable *var,
                                        DataSegment *ds, int depth, NodePool *pool);
void search_closest_points_optimized(KDTree *root, Variable *var, DataSegment *ds,
                                     ClosestPoint *closest, int target_id,
                                     int depth, int *found);
// */
#endif
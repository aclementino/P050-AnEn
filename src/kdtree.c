// /* --- kdtree 3 ---
#define _GNU_SOURCE
#include "kdtree.h"

NodePool *create_node_pool()
{
    NodePool *pool = (NodePool *)malloc(sizeof(NodePool));
    pool->next_available = 0;
    return pool;
}

KDTree *allocate_node_from_pool(NodePool *pool, int window_id)
{
    if (pool->next_available >= NODE_POOL_SIZE)
    {
        // Pool is full, fall back to malloc
        return create_kdt_node(window_id);
    }

    KDTree *node = &(pool->nodes[pool->next_available++]);
    node->window_id = window_id;
    node->left = NULL;
    node->right = NULL;
    return node;
}

void reset_node_pool(NodePool *pool)
{
    pool->next_available = 0;
}

void free_node_pool(NodePool *pool)
{
    free(pool);
}

// Standard node creation function (fallback)
KDTree *create_kdt_node(int window_id)
{
    KDTree *node = (KDTree *)malloc(sizeof(KDTree));
    node->window_id = window_id;
    node->left = NULL;
    node->right = NULL;
    return node;
}

// // Improved comparison context
// typedef struct
// {
//     Variable *var;
//     int axis;
//     int k; // Offset for window indexing
// } SortContext;

// Optimized comparison function
int compare_wrapper(const void *a, const void *b, void *context)
{
    SortContext *ctx = (SortContext *)context;
    const int *window_a = (const int *)a;
    const int *window_b = (const int *)b;
    int axis = ctx->axis;
    Variable *var = ctx->var;
    int k = ctx->k;

    // Optimize by directly comparing values based on type
    switch (var->type)
    {
    case NC_FLOAT:
    {
        float val_a = ((float *)var->data)[(*window_a - k) + axis];
        float val_b = ((float *)var->data)[(*window_b - k) + axis];
        if (val_a < val_b)
            return -1;
        if (val_a > val_b)
            return 1;
        return 0;
    }
    case NC_DOUBLE:
    {
        double val_a = ((double *)var->data)[(*window_a - k) + axis];
        double val_b = ((double *)var->data)[(*window_b - k) + axis];
        if (val_a < val_b)
            return -1;
        if (val_a > val_b)
            return 1;
        return 0;
    }
    // Add cases for other types
    default:
        return 0;
    }
}

// Optimized sorting function using qsort_r where available
void sort_points_by_axis(int *points, int n, Variable *var, int axis, int k)
{
    SortContext ctx = {var, axis, k};

#if defined(__GLIBC__) || defined(__APPLE__)
#ifdef __GLIBC__
    // GNU libc version
    qsort_r(points, n, sizeof(int), compare_wrapper, &ctx);
#else
    // BSD/macOS version has different argument order
    qsort_r(points, n, sizeof(int), &ctx, compare_wrapper);
#endif
#else
    // Fallback to an optimized quicksort implementation
    // This is a simplified implementation - a full quicksort would be better
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < n - i - 1; j++)
        {
            if (compare_wrapper(&points[j], &points[j + 1], &ctx) > 0)
            {
                // Swap
                int temp = points[j];
                points[j] = points[j + 1];
                points[j + 1] = temp;
            }
        }
    }
#endif
}

// Optimized median finding using median-of-medians algorithm for better pivot selection
int select_median(int *arr, int n, Variable *var, int axis, int k)
{
    if (n <= 5)
    {
        // For small arrays, just sort and return the median
        sort_points_by_axis(arr, n, var, axis, k);
        return arr[n / 2];
    }

    // For larger arrays, use the median-of-medians approach
    int num_medians = (n + 4) / 5; // Ceiling of n/5
    int *medians = (int *)malloc(num_medians * sizeof(int));

    // Find medians of groups of 5
    for (int i = 0; i < num_medians; i++)
    {
        int start = i * 5;
        int end = (i + 1) * 5 > n ? n : (i + 1) * 5;
        int *subarray = &arr[start];
        int subarray_size = end - start;

        sort_points_by_axis(subarray, subarray_size, var, axis, k);
        medians[i] = subarray[subarray_size / 2];
    }

    // Find the median of medians
    int median_of_medians = select_median(medians, num_medians, var, axis, k);
    free(medians);

    // Partition the array around the median of medians
    int i = 0;
    for (int j = 0; j < n; j++)
    {
        if (arr[j] == median_of_medians)
        {
            // Swap arr[i] and arr[j]
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
            i++;
            break;
        }
    }

    return median_of_medians;
}

// Improved balanced KD-tree builder
KDTree *build_balanced_kdtree(int *window_ids, int n, Variable *var, DataSegment *ds, int depth, NodePool *pool)
{
    if (n <= 0)
        return NULL;

    int axis = depth % ds->win_size; // Select axis based on depth

    // Select median using an optimized approach
    int median_idx = n / 2;
    sort_points_by_axis(window_ids, n, var, axis, ds->k);

    // Create node with median point
    KDTree *node;
    if (pool)
    {
        node = allocate_node_from_pool(pool, window_ids[median_idx]);
    }
    else
    {
        node = create_kdt_node(window_ids[median_idx]);
    }

    // Recursively build subtrees
    node->left = build_balanced_kdtree(window_ids, median_idx, var, ds, depth + 1, pool);
    node->right = build_balanced_kdtree(&window_ids[median_idx + 1], n - median_idx - 1, var, ds, depth + 1, pool);

    return node;
}

// Otimização 1: Usar select_median ao invés de sort completo
KDTree *build_optimized_balanced_kdtree(int *window_ids, int n, Variable *var,
                                        DataSegment *ds, int depth, NodePool *pool)
{
    if (n <= 0)
        return NULL;

    int axis = depth % ds->win_size;

    // OTIMIZAÇÃO: Usar select_median O(n) ao invés de sort O(n log n)
    int median_value = select_median(window_ids, n, var, axis, ds->k);

    // Particionar em torno da mediana
    int median_idx = partition_around_value(window_ids, n, median_value, var, axis, ds->k);

    KDTree *node = allocate_node_from_pool(pool, window_ids[median_idx]);

    // Recursão
    node->left = build_optimized_balanced_kdtree(window_ids, median_idx,
                                                 var, ds, depth + 1, pool);
    node->right = build_optimized_balanced_kdtree(&window_ids[median_idx + 1],
                                                  n - median_idx - 1,
                                                  var, ds, depth + 1, pool);
    return node;
}

// Optimized KD-tree creation
KDTree *create_balanced_kdtree(Variable *var, DataSegment *ds)
{
    int total_windows = ds->win_count;
    int *window_ids = (int *)malloc(total_windows * sizeof(int));

    // Initialize window IDs
    for (int i = 0; i < total_windows; i++)
    {
        window_ids[i] = i + ds->k; // Assuming window IDs start from ds->k
    }

    // Create a node pool for efficient memory allocation
    NodePool *pool = create_node_pool();

    // Build the balanced tree
    KDTree *root = build_balanced_kdtree(window_ids, total_windows, var, ds, 0, pool);

    free(window_ids);
    // Note: We don't free the pool here as it's needed for the nodes
    return root;
}

// Optimized insert function that maintains better balance
KDTree *insert_kdt_node(KDTree *root, Variable *var, DataSegment *ds, int window_id, int depth)
{
    if (root == NULL)
        return create_kdt_node(window_id);

    int axis = depth % ds->win_size;

    // Optimized value comparison
    int comparison = 0;
    switch (var->type)
    {
    case NC_FLOAT:
    {
        float val_a = ((float *)var->data)[(window_id - ds->k) + axis];
        float val_b = ((float *)var->data)[(root->window_id - ds->k) + axis];
        comparison = (val_a < val_b) ? -1 : 1;
        break;
    }
    case NC_DOUBLE:
    {
        double val_a = ((double *)var->data)[(window_id - ds->k) + axis];
        double val_b = ((double *)var->data)[(root->window_id - ds->k) + axis];
        comparison = (val_a < val_b) ? -1 : 1;
        break;
    }
        // Add implementations for other types
    }

    if (comparison < 0)
        root->left = insert_kdt_node(root->left, var, ds, window_id, depth + 1);
    else
        root->right = insert_kdt_node(root->right, var, ds, window_id, depth + 1);

    // Check if tree needs rebalancing after insertion
    if (!is_kdtree_balanced(root))
        // If tree is severely unbalanced, trigger rebalancing
        if (get_tree_height(root) > log2(count_nodes(root) + 1) * 1.5)
            return rebalance_kdtree(root, var, ds);

    return root;
}

// Optimized distance calculation (squared distance to avoid sqrt)
double squared_distance_kdtree(KDTree *root, Variable *var, DataSegment *ds, int window_id, int target_id)
{
    double sum = 0.0;

    // Optimize for the most common type (NC_FLOAT)
    if (var->type == NC_FLOAT)
    {
        float *data = (float *)var->data;
        int window_offset = window_id - ds->k;
        int target_offset = target_id - ds->k;

        for (int i = 0; i < ds->win_size; i++)
        {
            float diff = data[window_offset + i] - data[target_offset + i];
            sum += diff * diff;

            // Early termination if sum exceeds the current best distance
            if (sum > ds->current_best_distance && ds->current_best_distance > 0)
                return INFINITY;
        }
    }
    else
    {
        // Handle other types
        switch (var->type)
        {
        case NC_DOUBLE:
        {
            double *data = (double *)var->data;
            int window_offset = window_id - ds->k;
            int target_offset = target_id - ds->k;

            for (int i = 0; i < ds->win_size; i++)
            {
                double diff = data[window_offset + i] - data[target_offset + i];
                sum += diff * diff;
                if (sum > ds->current_best_distance && ds->current_best_distance > 0)
                {
                    return INFINITY;
                }
            }
            break;
        }
            // Add cases for other types
        }
    }

    return sum;
}

// Calculate distance and apply sqrt only when needed
double monache_metric_kdtree(KDTree *root, Variable *var, DataSegment *ds, int window_id, int target_id)
{
    double squared_dist = squared_distance_kdtree(root, var, ds, window_id, target_id);

    if (isnan(squared_dist) || isinf(squared_dist))
    {
        return squared_dist;
    }

    return sqrt(squared_dist);
}

double monache_metric_super_window_kdtree(KDTree *root, NetCDF *file, DataSegment *ds, int window_id, int target_id, int i)
{
    double sum = 0.0;

    for (int f = 0; f < (ds->argc - 1); f++) // Soma sobre séries temporais
    {
        switch (file[f].var[i].type)
        {
            MONACHESW_KDTREE(NC_INT);
            MONACHESW_KDTREE(NC_FLOAT);
            MONACHESW_KDTREE(NC_DOUBLE);
            MONACHESW_KDTREE(NC_UINT);
            MONACHESW_KDTREE(NC_INT64);
            MONACHESW_KDTREE(NC_UINT64);
        }
    }

    // return isnan(sum) ? NAN : sum;
    return isnan(sum) ? NAN : sqrt(sum);
}

// Optimized minimum distance between point and hyperrectangle (for pruning)
double min_distance_to_hyperrect(Variable *var, DataSegment *ds, int target_id, KDTree *node, int depth)
{
    if (node == NULL)
        return INFINITY;

    int axis = depth % ds->win_size;
    double sum = 0.0;

    // Calculate minimum distance to the hyperrectangle of this subtree
    if (var->type == NC_FLOAT)
    {
        float *data = (float *)var->data;
        float node_val = data[(node->window_id - ds->k) + axis];
        float target_val = data[(target_id - ds->k) + axis];

        if (target_val < node_val)
        {
            // Target is on the left side
            double diff = target_val - node_val;
            sum += diff * diff;
        }
    }
    else if (var->type == NC_DOUBLE)
    {
        double *data = (double *)var->data;
        double node_val = data[(node->window_id - ds->k) + axis];
        double target_val = data[(target_id - ds->k) + axis];

        if (target_val < node_val)
        {
            // Target is on the left side
            double diff = target_val - node_val;
            sum += diff * diff;
        }
    }

    return sum;
}

// Optimized closest points search with a priority queue
void search_closest_points(KDTree *root, Variable *var, DataSegment *ds, ClosestPoint *closest, int target_id, int depth, int *found)
{
    if (root == NULL)
        return;

    int axis = depth % ds->win_size;

    // Calculate actual distance
    double squared_dist = squared_distance_kdtree(root, var, ds, root->window_id, target_id);

    // If distance is finite, consider this point
    if (!isnan(squared_dist) && !isinf(squared_dist))
    {
        double distance = sqrt(squared_dist);

        if (*found < ds->num_Na)
        {
            closest[*found].window_index = root->window_id;
            closest[*found].distance = distance;
            (*found)++;

            if (*found == ds->num_Na)
            {
                // When we've found enough points, sort them
                qsort(closest, ds->num_Na, sizeof(ClosestPoint), compare_near_point);
                // Save the current best distance for early termination
                ds->current_best_distance = closest[0].distance * closest[0].distance;
            }
        }
        else if (distance < closest[0].distance)
        {
            // Found a better point
            closest[0].window_index = root->window_id;
            closest[0].distance = distance;
            qsort(closest, ds->num_Na, sizeof(ClosestPoint), compare_near_point);
            // Update the current best distance
            ds->current_best_distance = closest[0].distance * closest[0].distance;
        }
    }

    // Determine which child to visit first
    KDTree *first_child = NULL;
    KDTree *second_child = NULL;

    // Optimized axis-value comparison
    if (var->type == NC_FLOAT)
    {
        float target_val = ((float *)var->data)[(target_id - ds->k) + axis];
        float node_val = ((float *)var->data)[(root->window_id - ds->k) + axis];

        if (target_val < node_val)
        {
            first_child = root->left;
            second_child = root->right;
        }
        else
        {
            first_child = root->right;
            second_child = root->left;
        }
    }
    else if (var->type == NC_DOUBLE)
    {
        // Similar implementation for double
        double target_val = ((double *)var->data)[(target_id - ds->k) + axis];
        double node_val = ((double *)var->data)[(root->window_id - ds->k) + axis];

        if (target_val < node_val)
        {
            first_child = root->left;
            second_child = root->right;
        }
        else
        {
            first_child = root->right;
            second_child = root->left;
        }
    }

    // Visit the first child
    search_closest_points(first_child, var, ds, closest, target_id, depth + 1, found);

    // Visit the second child only if necessary
    if (*found < ds->num_Na)
    {
        // We need more points, so visit the second child
        search_closest_points(second_child, var, ds, closest, target_id, depth + 1, found);
    }
    else
    {
        // Calculate the minimum possible distance to the hyperrectangle
        double min_dist = 0.0;

        if (var->type == NC_FLOAT)
        {
            float target_val = ((float *)var->data)[(target_id - ds->k) + axis];
            float node_val = ((float *)var->data)[(root->window_id - ds->k) + axis];
            float diff = target_val - node_val;
            min_dist = diff * diff;
        }
        else if (var->type == NC_DOUBLE)
        {
            // Similar implementation for double
            double target_val = ((double *)var->data)[(target_id - ds->k) + axis];
            double node_val = ((double *)var->data)[(root->window_id - ds->k) + axis];
            double diff = target_val - node_val;
            min_dist = diff * diff;
        }

        // If the minimum distance is less than our current best distance, we need to search there too
        if (min_dist < ds->current_best_distance)
        {
            search_closest_points(second_child, var, ds, closest, target_id, depth + 1, found);
        }
    }
}

// Fixed and optimized version of search_closest_points_super_window
void search_closest_points_super_window(KDTree *root, NetCDF *file, DataSegment *ds, ClosestPoint *closest, int target_id, int depth, int i, int *found)
{
    if (root == NULL)
        return;

    int axis = depth % ds->win_size;

    double distance = monache_metric_super_window_kdtree(root, file, ds, root->window_id, target_id, i);

    if (!isnan(distance))
    {
        if (*found < ds->num_Na)
        {
            closest[*found].window_index = root->window_id;
            closest[*found].distance = distance;
            (*found)++;
            if (*found == ds->num_Na)
            {
                qsort(closest, ds->num_Na, sizeof(ClosestPoint), compare_near_point);
                ds->current_best_distance = closest[0].distance * closest[0].distance;
            }
        }
        else if (distance < closest[0].distance)
        {
            closest[0].window_index = root->window_id;
            closest[0].distance = distance;
            qsort(closest, ds->num_Na, sizeof(ClosestPoint), compare_near_point);
            ds->current_best_distance = closest[0].distance * closest[0].distance;
        }
    }

    KDTree *first_child = NULL;
    KDTree *second_child = NULL;

    // Optimized comparison based on type
    switch (file->var[i].type)
    {
    case NC_FLOAT:
    {
        float target_val = ((float *)file->var[i].data)[(target_id - ds->k) + axis];
        float node_val = ((float *)file->var[i].data)[(root->window_id - ds->k) + axis];

        if (target_val < node_val)
        {
            first_child = root->left;
            second_child = root->right;
        }
        else
        {
            first_child = root->right;
            second_child = root->left;
        }
        break;
    }
    case NC_DOUBLE:
    {
        double target_val = ((double *)file->var[i].data)[(target_id - ds->k) + axis];
        double node_val = ((double *)file->var[i].data)[(root->window_id - ds->k) + axis];

        if (target_val < node_val)
        {
            first_child = root->left;
            second_child = root->right;
        }
        else
        {
            first_child = root->right;
            second_child = root->left;
        }
        break;
    }
        // Add cases for other types
    }

    // Visit first child
    search_closest_points_super_window(first_child, file, ds, closest, target_id, depth + 1, i, found);

    // Visit second child only if necessary
    if (*found < ds->num_Na || min_distance_to_hyperrect(&file->var[i], ds, target_id, root, depth) < ds->current_best_distance)
    {
        search_closest_points_super_window(second_child, file, ds, closest, target_id, depth + 1, i, found);
    }
}

// Otimização 2: Busca com poda avançada
void search_closest_points_optimized(KDTree *root, Variable *var, DataSegment *ds,
                                     ClosestPoint *closest, int target_id,
                                     int depth, int *found)
{
    if (root == NULL)
        return;

    // OTIMIZAÇÃO: Poda precoce usando distância ao hiper-retângulo
    if (*found >= ds->num_Na)
    {
        double min_hyperrect_dist = min_distance_to_hyperrect(var, ds, target_id, root, depth);
        if (min_hyperrect_dist >= ds->current_best_distance)
        {
            return; // Pode podar toda esta subárvore
        }
    }

    // Continuar com a busca normal...
    double squared_dist = squared_distance_kdtree(root, var, ds, root->window_id, target_id);

    if (!isnan(squared_dist) && !isinf(squared_dist))
    {
        double distance = sqrt(squared_dist);

        if (*found < ds->num_Na)
        {
            closest[*found].window_index = root->window_id;
            closest[*found].distance = distance;
            (*found)++;

            if (*found == ds->num_Na)
            {
                qsort(closest, ds->num_Na, sizeof(ClosestPoint), compare_near_point);
                ds->current_best_distance = closest[0].distance * closest[0].distance;
            }
        }
        else if (distance < closest[0].distance)
        {
            closest[0].window_index = root->window_id;
            closest[0].distance = distance;
            qsort(closest, ds->num_Na, sizeof(ClosestPoint), compare_near_point);
            ds->current_best_distance = closest[0].distance * closest[0].distance;
        }
    }

    // Determinar ordem de visita dos filhos
    int axis = depth % ds->win_size;
    KDTree *first_child, *second_child;

    if (var->type == NC_FLOAT)
    {
        float target_val = ((float *)var->data)[(target_id - ds->k) + axis];
        float node_val = ((float *)var->data)[(root->window_id - ds->k) + axis];

        if (target_val < node_val)
        {
            first_child = root->left;
            second_child = root->right;
        }
        else
        {
            first_child = root->right;
            second_child = root->left;
        }
    }

    // Visitar primeiro filho
    search_closest_points_optimized(first_child, var, ds, closest, target_id, depth + 1, found);

    // Visitar segundo filho apenas se necessário (com poda otimizada)
    if (*found < ds->num_Na)
    {
        search_closest_points_optimized(second_child, var, ds, closest, target_id, depth + 1, found);
    }
    else
    {
        // OTIMIZAÇÃO: Poda mais precisa usando hiper-retângulo
        double hyperrect_dist = min_distance_to_hyperrect(var, ds, target_id, second_child, depth + 1);
        if (hyperrect_dist < ds->current_best_distance)
        {
            search_closest_points_optimized(second_child, var, ds, closest, target_id, depth + 1, found);
        }
    }
}

// Helper function to collect window IDs from the tree
void collect_window_ids(KDTree *root, int *window_ids, int *index)
{
    if (root == NULL)
        return;

    collect_window_ids(root->left, window_ids, index);
    window_ids[(*index)++] = root->window_id;
    collect_window_ids(root->right, window_ids, index);
}

// Efficient tree deallocation
void deallocate_kdtree(KDTree *root)
{
    if (root == NULL)
        return;

    deallocate_kdtree(root->left);
    deallocate_kdtree(root->right);
    free(root);
}

// Improved rebalance function with node pool
KDTree *rebalance_kdtree(KDTree *root, Variable *var, DataSegment *ds)
{
    if (root == NULL)
        return NULL;

    // Count nodes in the tree
    int node_count = count_nodes(root);
    int *window_ids = (int *)malloc(node_count * sizeof(int));

    // Collect all window IDs through inorder traversal
    int index = 0;
    collect_window_ids(root, window_ids, &index);

    // Create a node pool for efficient memory allocation
    NodePool *pool = create_node_pool();

    // Deallocate old tree
    deallocate_kdtree(root);

    // Build a new balanced tree
    KDTree *new_root = build_balanced_kdtree(window_ids, node_count, var, ds, 0, pool);

    free(window_ids);
    return new_root;
}

// Function to check if the k-d tree is balanced
bool is_kdtree_balanced(KDTree *root)
{
    if (root == NULL)
        return true;

    int height = get_tree_height(root);
    int total_nodes = count_nodes(root);
    double theoretical_min_height = log2(total_nodes + 1);

    if (height <= (theoretical_min_height * 1.5))
        return true;

    return false;
}

// Function to get the height of a k-d tree
int get_tree_height(KDTree *node)
{
    if (node == NULL)
        return 0;

    int left_height = get_tree_height(node->left);
    int right_height = get_tree_height(node->right);

    return 1 + (left_height > right_height ? left_height : right_height);
}

// Function to count the number of nodes in the tree
int count_nodes(KDTree *node)
{
    if (node == NULL)
        return 0;

    return 1 + count_nodes(node->left) + count_nodes(node->right);
}

// Optimized comparison function for qsort
int compare_near_point(const void *x, const void *y)
{
    ClosestPoint *point1 = (ClosestPoint *)x;
    ClosestPoint *point2 = (ClosestPoint *)y;

    if (point1->distance < point2->distance)
        return 1;
    if (point1->distance > point2->distance)
        return -1;
    return 0;
}

// Enhanced diagnostic function
void diagnose_tree_balance(KDTree *root)
{
    if (root == NULL)
    {
        printf("Tree is empty\n");
        return;
    }

    int height = get_tree_height(root);
    int total_nodes = count_nodes(root);
    double theoretical_min_height = log2(total_nodes + 1);

    printf("Tree Statistics:\n");
    printf("  Total nodes: %d\n", total_nodes);
    printf("  Height: %d\n", height);
    printf("  Theoretical minimum height: %.2f\n", theoretical_min_height);
    printf("  Height ratio: %.2f\n", (double)height / theoretical_min_height);

    if (is_kdtree_balanced(root))
    {
        printf("  Status: BALANCED\n");
    }
    else
    {
        printf("  Status: UNBALANCED\n");
        printf("  Recommendation: Call rebalance_kdtree() function\n");
    }

    int left_height = get_tree_height(root->left);
    int right_height = get_tree_height(root->right);
    int left_nodes = count_nodes(root->left);
    int right_nodes = count_nodes(root->right);

    printf("  Left subtree: %d nodes, height %d\n", left_nodes, left_height);
    printf("  Right subtree: %d nodes, height %d\n", right_nodes, right_height);
    printf("  Height difference: %d\n", abs(left_height - right_height));
}

// Function to visualize the tree structure
void visualize_kdtree(KDTree *root, int depth)
{
    if (root == NULL)
        return;

    for (int i = 0; i < depth; i++)
        printf("  ");
    printf("|- Node: %d\n", root->window_id);

    visualize_kdtree(root->left, depth + 1);
    visualize_kdtree(root->right, depth + 1);
}

// Função auxiliar para particionamento
int partition_around_value(int *arr, int n, int pivot_value, Variable *var, int axis, int k)
{
    int pivot_idx = -1;

    // Encontrar o índice do pivot
    for (int i = 0; i < n; i++)
    {
        if (arr[i] == pivot_value)
        {
            pivot_idx = i;
            break;
        }
    }

    if (pivot_idx == -1)
        return n / 2; // Fallback

    // Mover pivot para o final
    int temp = arr[pivot_idx];
    arr[pivot_idx] = arr[n - 1];
    arr[n - 1] = temp;

    // Particionar
    int store_idx = 0;
    for (int i = 0; i < n - 1; i++)
    {
        if (compare_values(arr[i], pivot_value, var, axis, k) <= 0)
        {
            temp = arr[i];
            arr[i] = arr[store_idx];
            arr[store_idx] = temp;
            store_idx++;
        }
    }

    // Mover pivot para posição final
    temp = arr[store_idx];
    arr[store_idx] = arr[n - 1];
    arr[n - 1] = temp;

    return store_idx;
}

int compare_values(int window_a, int window_b, Variable *var, int axis, int k)
{
    if (var->type == NC_FLOAT)
    {
        float val_a = ((float *)var->data)[(window_a - k) + axis];
        float val_b = ((float *)var->data)[(window_b - k) + axis];
        if (val_a < val_b)
            return -1;
        if (val_a > val_b)
            return 1;
        return 0;
    }
    // Adicionar outros tipos conforme necessário
    return 0;
}

#include <preprocess.h>

time_t convert_time(char *rawtime)
{
    struct tm time_info;
    memset(&time_info, 0, sizeof(struct tm));

    if (sscanf(rawtime, "%4d-%2d-%2dT%2d:%2d:%2d",
               &time_info.tm_year, &time_info.tm_mon, &time_info.tm_mday,
               &time_info.tm_hour, &time_info.tm_min, &time_info.tm_sec) != 6)
    {
        printf("Application: Error converting ISO 8601 string.\n");
        return -1;
    }

    time_info.tm_year -= 1900;
    time_info.tm_mon -= 1;

    // mktime => local time; timegm => UTC time
    return timegm(&time_info) / 60;
}

int binary_search(NetCDF *file, int target)
{
    int *values = file[0].var[0].data;

    int init = 0, end = file[0].dim->len - 1;

    while (init <= end)
    {
        int middle = init + (end - init) / 2;

        if (values[middle] == target)
            return middle;

        if (values[middle] < target)
            init = middle + 1;
        else
            end = middle - 1;
    }

    return -1;
}

void analyze_data(NetCDF *file, DataSegment *ds, process_func func)
{
    for (int i = 0; i < (ds->argc); i++)
    {
        if (file == NULL || file[i].var == NULL)
        {
            printf("Application(analyze_data): No data in file %i.\n", (i + 1));
            return;
        }

        ds->indice_generic = i;

        func(&file[i], ds); // Chama a função de processamento para cada item
    }
}

void print_data_values(NetCDF *file, DataSegment *ds)
{
    Variable *var = NULL;

    // for (int i = 0; i < file->dim->len; i++)
    for (int i = 703440; i < 704109; i++)
    {
        if (i == 0)
        {
            printf("file data %i:\n", i + 1);

            for (int t = 0; t < file->nvars; t++)
            {
                t == 0 ? printf("%s\t\t", file->var[t].name)
                       : printf("%s\t", file->var[t].name);
            }
            printf("\n");
        }

        for (int j = 0; j < file->nvars; j++)
        {
            var = &file->var[j];
            switch (var->type)
            {
                PRINT_LN_VALUES(NC_BYTE, "%hhd");
                PRINT_LN_VALUES(NC_CHAR, "%c");
                PRINT_LN_VALUES(NC_SHORT, "%hd");
                PRINT_LN_VALUES(NC_INT, "%i");
                PRINT_LN_VALUES(NC_FLOAT, "%.1f");
                PRINT_LN_VALUES(NC_DOUBLE, "%lf");
                PRINT_LN_VALUES(NC_UBYTE, "%hhu");
                PRINT_LN_VALUES(NC_USHORT, "%hu");
                PRINT_LN_VALUES(NC_UINT, "%u");
                PRINT_LN_VALUES(NC_INT64, "%lld");
                PRINT_LN_VALUES(NC_UINT64, "%llu");
                PRINT_LN_VALUES(NC_STRING, "%s");
            default:
            {
                fprintf(stderr, "Tipo não suportado: %i", var->type);
                exit(EXIT_FAILURE);
            }
            break;
            }
        }
        printf("\n");
    }
}

void count_invalid_values(NetCDF *file, DataSegment *ds)
{
    Variable *var = NULL;

    for (int i = 1; i < file->nvars; i++)
    {
        int invalid_count = 0;

        var = &file->var[i];
        switch (var->type)
        {
            VALIDATE_DATA(NC_BYTE);
            VALIDATE_DATA(NC_CHAR);
            VALIDATE_DATA(NC_SHORT);
            VALIDATE_DATA(NC_INT);
            VALIDATE_DATA(NC_FLOAT);
            VALIDATE_DATA(NC_DOUBLE);
            VALIDATE_DATA(NC_UBYTE);
            VALIDATE_DATA(NC_USHORT);
            VALIDATE_DATA(NC_UINT);
            VALIDATE_DATA(NC_INT64);
            VALIDATE_DATA(NC_UINT64);
        default:
        {
            fprintf(stderr, "Tipo não suportado: %i", var->type);
            exit(EXIT_FAILURE);
        }
        break;
        }

        var->invalid_count = invalid_count;

        var->invalid_percentage = 100 - (((double)(file->dim->len - invalid_count) / file->dim->len) * 100);
    }
}

void print_info_percentage(NetCDF *file, DataSegment *ds)
{
    for (int i = 0; i < file->nvars; i++)
        printf("%i. invalid_count: %i  \tinvalid_percentage: %.2lf\n",
               ds->indice_generic, file->var[i].invalid_count, file->var[i].invalid_percentage);
}

void interpolation_values(NetCDF *file, DataSegment *ds)
{
    Variable *var_zero = NULL;
    Variable *var = NULL;
    int d_changes = 0;
    double d_time[2] = {0, 0};
    double d_interp[2] = {0, 0};
    int d_index[2] = {0, 0};

    // Initialize the interpolation objects
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_interp *interp = gsl_interp_alloc(gsl_interp_linear, 2);

    var_zero = &file->var[0];

    // printf("file: %i\n", ds->indice_generic);
    for (int i = 1; i < file->nvars; i++)
    {
        // printf("var: %i | ", i);

        var = &file->var[i];
        int count = 0;
        if ((int)var->invalid_percentage != 0 && (int)var->invalid_percentage != 100)
        {
            if (ds->indice_generic == 0)
                for (int j = 0; j < ds->start_prediction; j++)
                {
                    switch (var->type)
                    {
                        IF_ISNAN(NC_BYTE);
                        IF_ISNAN(NC_CHAR);
                        IF_ISNAN(NC_SHORT);
                        IF_ISNAN(NC_INT);
                        IF_ISNAN(NC_FLOAT);
                        IF_ISNAN(NC_DOUBLE);
                        IF_ISNAN(NC_UBYTE);
                        IF_ISNAN(NC_USHORT);
                        IF_ISNAN(NC_UINT);
                        IF_ISNAN(NC_INT64);
                        IF_ISNAN(NC_UINT64);
                    }
                    if (d_time[1] && (d_time[0] < d_time[1]))
                        switch (var->type)
                        {
                            INTERP_DATA(NC_BYTE);
                            INTERP_DATA(NC_CHAR);
                            INTERP_DATA(NC_SHORT);
                            INTERP_DATA(NC_INT);
                            INTERP_DATA(NC_FLOAT);
                            INTERP_DATA(NC_DOUBLE);
                            INTERP_DATA(NC_UBYTE);
                            INTERP_DATA(NC_USHORT);
                            INTERP_DATA(NC_UINT);
                            INTERP_DATA(NC_INT64);
                            INTERP_DATA(NC_UINT64);
                        }
                }
            else
            {
                for (int j = 0; j < ds->start_prediction; j++)
                {
                    switch (var->type)
                    {
                        IF_ISNAN(NC_BYTE);
                        IF_ISNAN(NC_CHAR);
                        IF_ISNAN(NC_SHORT);
                        IF_ISNAN(NC_INT);
                        IF_ISNAN(NC_FLOAT);
                        IF_ISNAN(NC_DOUBLE);
                        IF_ISNAN(NC_UBYTE);
                        IF_ISNAN(NC_USHORT);
                        IF_ISNAN(NC_UINT);
                        IF_ISNAN(NC_INT64);
                        IF_ISNAN(NC_UINT64);
                    }
                    if (d_time[1] && (d_time[0] < d_time[1]))
                        switch (var->type)
                        {
                            INTERP_DATA(NC_BYTE);
                            INTERP_DATA(NC_CHAR);
                            INTERP_DATA(NC_SHORT);
                            INTERP_DATA(NC_INT);
                            INTERP_DATA(NC_FLOAT);
                            INTERP_DATA(NC_DOUBLE);
                            INTERP_DATA(NC_UBYTE);
                            INTERP_DATA(NC_USHORT);
                            INTERP_DATA(NC_UINT);
                            INTERP_DATA(NC_INT64);
                            INTERP_DATA(NC_UINT64);
                        }
                }
                for (int j = ds->start_prediction; j < ds->end_prediction + ds->k; j++)
                {
                    switch (var->type)
                    {
                        IF_ISNAN(NC_BYTE);
                        IF_ISNAN(NC_CHAR);
                        IF_ISNAN(NC_SHORT);
                        IF_ISNAN(NC_INT);
                        IF_ISNAN(NC_FLOAT);
                        IF_ISNAN(NC_DOUBLE);
                        IF_ISNAN(NC_UBYTE);
                        IF_ISNAN(NC_USHORT);
                        IF_ISNAN(NC_UINT);
                        IF_ISNAN(NC_INT64);
                        IF_ISNAN(NC_UINT64);
                    }
                    if (d_time[1] && (d_time[0] < d_time[1]))
                        switch (var->type)
                        {
                            INTERP_DATA(NC_BYTE);
                            INTERP_DATA(NC_CHAR);
                            INTERP_DATA(NC_SHORT);
                            INTERP_DATA(NC_INT);
                            INTERP_DATA(NC_FLOAT);
                            INTERP_DATA(NC_DOUBLE);
                            INTERP_DATA(NC_UBYTE);
                            INTERP_DATA(NC_USHORT);
                            INTERP_DATA(NC_UINT);
                            INTERP_DATA(NC_INT64);
                            INTERP_DATA(NC_UINT64);
                        }
                }
            }
            d_changes = 0;
            continue;
        }
    }

    // Free resource
    gsl_interp_free(interp);
    gsl_interp_accel_free(acc);
}

void count_valid_window(NetCDF *file, DataSegment *ds)
{
    Variable *var = NULL;

    int invalid_count = 0;
    int num_isvalid = 0;
    int start = 0;
    int end = 0;

    for (int i = 1; i < file->nvars; i++)
    {
        var = &file->var[i];

        /* ------------------ Training ------------------ */
        start = ds->start_training;
        end = ds->end_training;

        switch (var->type)
        {
            COUNT_VALID_WIN(NC_BYTE);
            COUNT_VALID_WIN(NC_CHAR);
            COUNT_VALID_WIN(NC_SHORT);
            COUNT_VALID_WIN(NC_INT);
            COUNT_VALID_WIN(NC_FLOAT);
            COUNT_VALID_WIN(NC_DOUBLE);
            COUNT_VALID_WIN(NC_UBYTE);
            COUNT_VALID_WIN(NC_USHORT);
            COUNT_VALID_WIN(NC_UINT);
            COUNT_VALID_WIN(NC_INT64);
            COUNT_VALID_WIN(NC_UINT64);
        default:
        {
            fprintf(stderr, "Tipo não suportado: %i", var->type);
            exit(EXIT_FAILURE);
        }
        break;
        }

        var->num_win_valid_training = num_isvalid;

        /* ----------------- Prediction ----------------- */
        start = ds->start_prediction;
        end = ds->end_prediction;
        num_isvalid = 0;

        switch (var->type)
        {
            COUNT_VALID_WIN(NC_BYTE);
            COUNT_VALID_WIN(NC_CHAR);
            COUNT_VALID_WIN(NC_SHORT);
            COUNT_VALID_WIN(NC_INT);
            COUNT_VALID_WIN(NC_FLOAT);
            COUNT_VALID_WIN(NC_DOUBLE);
            COUNT_VALID_WIN(NC_UBYTE);
            COUNT_VALID_WIN(NC_USHORT);
            COUNT_VALID_WIN(NC_UINT);
            COUNT_VALID_WIN(NC_INT64);
            COUNT_VALID_WIN(NC_UINT64);
        default:
        {
            fprintf(stderr, "Tipo não suportado: %i", var->type);
            exit(EXIT_FAILURE);
        }
        break;
        }

        var->num_win_valid_prediction = num_isvalid;
    }
}
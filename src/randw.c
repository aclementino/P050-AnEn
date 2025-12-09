#include "randw.h"

/*
 * Function to handle NetCDF errors, by printing an error message and exiting
 * with a non-zero status.
 */
void handle_error(int status)
{
    if (status != NC_NOERR)
    {
        fprintf(stderr, "%s\n", nc_strerror(status));
        exit(1);
    }
}

NetCDF *create_struct(DataSegment *ds, char *argv[])
{
    NetCDF *file = (NetCDF *)malloc(ds->argc * sizeof(NetCDF));

    if (file == NULL)
    {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        return NULL;
    }

    for (int i = 0; i < (ds->argc); i++)
    {
        /*
         * Open the input NetCDF file
         * NC_NOWRITE tells NetCDF we want read-only access to the file
         */
        handle_error(nc_open(argv[i], NC_NOWRITE, &file[i].ncid_in));

        read_header_file(&file[i], ds);
        read_data_file(&file[i], ds);
        // printf("file: %s.\n", argv[i]);
    }

    return file;
}
void *allocate_memory(char type, size_t len)
{
    if (type < NC_BYTE || type > NC_STRING)
    {
        fprintf(stderr, "Error: Unknown data type.\n");
        exit(EXIT_FAILURE);
    }

    void *data = NULL;

    switch (type)
    {
        ALLOCATE_MEMORY(NC_BYTE, len);
        ALLOCATE_MEMORY(NC_CHAR, len);
        ALLOCATE_MEMORY(NC_SHORT, len);
        ALLOCATE_MEMORY(NC_INT, len);
        ALLOCATE_MEMORY(NC_FLOAT, len);
        ALLOCATE_MEMORY(NC_DOUBLE, len);
        ALLOCATE_MEMORY(NC_UBYTE, len);
        ALLOCATE_MEMORY(NC_USHORT, len);
        ALLOCATE_MEMORY(NC_UINT, len);
        ALLOCATE_MEMORY(NC_INT64, len);
        ALLOCATE_MEMORY(NC_UINT64, len);
        ALLOCATE_MEMORY(NC_STRING, len);
    default:
        break;
    }

    if (data == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory.\n");
        exit(EXIT_FAILURE);
    }
    return data;
}
void deallocate_memory(NetCDF *file, int argc)
{
    Variable *var = NULL;
    Dimension *dim = NULL;

    for (int i = 0; i < argc; i++)
    {
        dim = file[i].dim;
        var = file[i].var;
        
        free(dim);
        dim = NULL;
        file[i].dim = NULL;

        for (int j = 0; j < file[i].nvars; j++)
        {
            free(var[j].data);
            var[j].data = NULL;
        }
        free(var);
        var = NULL;
        file[i].var = NULL;
        
        /* Close the file, freeing all resources. */
        handle_error(nc_close(file[i].ncid_in));
    }
    free(file);
}
void read_dimensions(NetCDF *file)
{
    file->dim = (Dimension *)malloc(file->ndims * sizeof(Dimension));
    Dimension *dim = file->dim;
    /* Get dimensions */
    for (int i = 0; i < file->ndims; i++)
    {
        handle_error(nc_inq_dim(file->ncid_in, i, dim[i].name, &dim[i].len));
    }
}
void read_variables(NetCDF *file)
{
    file->var = (Variable *)malloc(file->nvars * sizeof(Variable));
    Variable *var = file->var;
    /* Get variables and yours attributes */
    for (int i = 0; i < file->nvars; i++)
    {
        /* Get var name */
        handle_error(nc_inq_varname(file->ncid_in, i, var[i].name));
        /* Get id type var */
        handle_error(nc_inq_varid(file->ncid_in, var[i].name, &var[i].id));
        /* Get type var */
        handle_error(nc_inq_vartype(file->ncid_in, i, &var[i].type));
    }
}
void read_header_file(NetCDF *file, DataSegment *ds)
{
    /* Get information from the input file */
    handle_error(nc_inq(file->ncid_in,
                        &file->ndims,
                        &file->nvars,
                        NULL,
                        NULL));

    /* Get base dimension for variables */
    read_dimensions(file);

    /* Get the number of variables in the file */
    read_variables(file);
}
void read_data_file(NetCDF *file, DataSegment *ds)
{
    Variable *var = file->var;

    // Copiando os dados das variáveis
    for (int i = 0; i < file->nvars; i++)
    {
        var[i].data = allocate_memory(var[i].type, file->dim->len);

        if (var[i].data == NULL)
        {
            fprintf(stderr, "Error: Variable memory allocation failed.\n");
            exit(1);
        }

        // Alocando memória para os dados e copiando
        switch (var[i].type)
        {
        case NC_BYTE:
        case NC_CHAR:
            nc_get_var_text(file->ncid_in, var[i].id, var[i].data);
            break;
        case NC_SHORT:
            nc_get_var_short(file->ncid_in, var[i].id, var[i].data);
            break;
        case NC_INT:
            nc_get_var_int(file->ncid_in, var[i].id, var[i].data);
            break;
        case NC_FLOAT:
            nc_get_var_float(file->ncid_in, var[i].id, var[i].data);
            break;
        case NC_DOUBLE:
            nc_get_var_double(file->ncid_in, var[i].id, var[i].data);
            break;
        case NC_UBYTE:
            nc_get_var_uchar(file->ncid_in, var[i].id, var[i].data);
            break;
        case NC_USHORT:
            nc_get_var_ushort(file->ncid_in, var[i].id, var[i].data);
            break;
        case NC_UINT:
            nc_get_var_uint(file->ncid_in, var[i].id, var[i].data);
            break;
        case NC_INT64:
            nc_get_var_longlong(file->ncid_in, var[i].id, var[i].data);
            break;
        case NC_UINT64:
            nc_get_var_ulonglong(file->ncid_in, var[i].id, var[i].data);
            break;
        case NC_STRING:
            nc_get_var_string(file->ncid_in, var[i].id, var[i].data);
            break;
        default:
            fprintf(stderr, "Unknown variable type.\n");
            exit(1);
        }
    }
}
void write_file(NetCDF *file, char *argv)
{
    time_t t = time(NULL); // Get the current time
    int ncid_out;
    int ndims, natts, unlimdimid;
    int dimids[NC_MAX_DIMS];
    char buffer[20];
    Variable *var = file->var;

    /* Convert to local date/time format */
    struct tm *tm_info = localtime(&t);

    /* Format the date and time */
    strftime(buffer, sizeof(buffer), "%Y%m%d%H%M", tm_info);

    /* Get information from the input file */
    handle_error(nc_inq(file->ncid_in, &ndims, NULL, &natts, &unlimdimid));

    /* Create the output NetCDF file */
    handle_error(nc_create(buffer, NC_CLOBBER, &ncid_out));

    /* Copy and Write the dimensions */
    for (int i = 0; i < ndims; i++)
    {
        char dim_name[NC_MAX_NAME + 1];
        size_t dimlen;
        handle_error(nc_inq_dim(file->ncid_in, i, dim_name, &dimlen));
        handle_error(nc_def_dim(ncid_out, dim_name, dimlen, &dimids[i]));
    }

    /* Copy and Write the variables and your attributes*/
    for (int i = 0; i < file->nvars; i++)
    {
        int var_ndims, var_natts, var_dimids[NC_MAX_VAR_DIMS];

        handle_error(nc_inq_var(file->ncid_in, i, NULL, NULL, &var_ndims, var_dimids, &var_natts));
        handle_error(nc_def_var(ncid_out, var[i].name, var[i].type, var_ndims, var_dimids, NULL));

        for (int j = 0; j < var_natts; j++)
        {
            char att_name[NC_MAX_NAME + 1];
            handle_error(nc_inq_attname(file->ncid_in, i, j, att_name));
            handle_error(nc_copy_att(file->ncid_in, i, att_name, ncid_out, i));
        }
    }

    /* Write the global attributes */
    for (int i = 0; i < natts; i++)
    {
        char att_name[NC_MAX_NAME + 1];
        handle_error(nc_inq_attname(file->ncid_in, NC_GLOBAL, i, att_name));
        handle_error(nc_copy_att(file->ncid_in, NC_GLOBAL, att_name, ncid_out, NC_GLOBAL));
    }

    /* Finalizing the output file definition */
    handle_error(nc_enddef(ncid_out));

    /* Copy the variables data */
    for (int i = 0; i < file->nvars; i++)
        handle_error(nc_put_var(ncid_out, i, var[i].data));

    /* Closing NetCDF file*/
    handle_error(nc_close(ncid_out));
    printf("Arquivo NetCDF copiado com sucesso para %s!\n", buffer);
}
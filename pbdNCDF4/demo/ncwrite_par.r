### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

library(pbdMPI, quiet = TRUE)
library(pbdNCDF4, quiet = TRUE)

### default of pbdMPI
rank <- comm.rank()
size <- comm.size()

### define data dimension.
n <- size * 4
p <- 5

### set up slab for this process
st <- c(rank * n/ size + 1, 1)
co <- c(n / size, p)

### generate local matrix portion
x <- matrix(rank + 1, nrow = n / size, ncol = p) +
     rep(1:p, each = n / size)

### set up dimensions for full matrix (not just local slab)
rdim <- ncdim_def("rows", "number", vals = 1:n)
cdim <- ncdim_def("columns", "number", vals = 1:p)

### define matrix variable in file (must not create full storage!!)
x.nc_var <- ncvar_def("testMatrix", "count",
                      list(rows = rdim, columns = cdim),
                      missval = -1, prec = "integer")

### create (collectively) in parallel a file with given dimensions
nc <- nc_create_par("test_par.nc", x.nc_var)
nc_var_par_access(nc, x.nc_var)

### write variable values to file
ncvar_put(nc, "testMatrix", as.vector(x), start = st, count = co)
nc_sync(nc)

### close file
nc_close(nc)

finalize()

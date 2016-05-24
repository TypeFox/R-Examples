### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

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

### create (collectively) a file with given dimensions
if(comm.rank() == 0){
  nc <- nc_create("test_ser.nc", x.nc_var)
  nc_close(nc)
}
barrier()

### write variable values to file
for(i in 0:(comm.size() - 1)){
  if(comm.rank() == i){
    nc <- nc_open("test_ser.nc", write = TRUE)
    ncvar_put(nc, "testMatrix", as.vector(x), start = st, count = co)
    nc_close(nc)
  }
  barrier()
}

finalize()

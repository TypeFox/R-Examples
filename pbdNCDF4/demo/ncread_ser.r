### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

library(pbdNCDF4, quiet = TRUE)

### default of pbdMPI
rank <- comm.rank()
size <- comm.size()

### serial open
if(! file.exists("test_ser.nc")){
  stop("test_ser.nc does not exist. Run ncwrite_ser to create.")
}
nc <- nc_open("test_ser.nc")

### get data dimension.
n <- length(nc$dim$rows$vals)
p <- length(nc$dim$columns$vals)

### set up slab for this process
st <- c(rank*n/size + 1, 1)
co <- c(n/size         , p)

### read variable values to file
x <- ncvar_get(nc, "testMatrix", start = st, count = co)

### Print results.
comm.cat("n = ", n, " p = ", p, "\n", sep = "")
comm.print(x, all.rank = TRUE)

### close file
nc_close(nc)
ncdump("test_ser.nc")

finalize()

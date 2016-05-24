library(pbdDEMO, quietly = TRUE)

# -------------------------------------
# Read example file
# -------------------------------------

### default of pbdMPI
rank <- comm.rank()
size <- comm.size()

file.path <- system.file("extra/data/x.nc", package = "pbdDEMO")
nc <- nc_open_par(file.path)
nc_var_par_access(nc, "testMatrix")

### get data dimension.
n <- length(nc$dim$rows$vals)
p <- length(nc$dim$columns$vals)

### set up slab for this process
st <- c(rank*n/size + 1, 1)
co <- c(n/size         , p)

### read variable values to file
x <- ncvar_get(nc, "testMatrix", start = st, count = co)

### close file
nc_close(nc)


# -------------------------------------
# Distributed matrix computations
# -------------------------------------

init.grid()

x <- matrix(as.double(x), dim(x)[1L], dim(x)[2L])

dx <- new("ddmatrix", Data=x, dim=c(n, p), ldim=dim(x), bldim=dim(x), ICTXT=2)

dx <- redistribute(dx, bldim=c(2,2), ICTXT=0)


### some matrix algebra
dx <- t(dx-2) %*% dx

dx <- redistribute(dx, ICTXT=2)


# -------------------------------------
# Write output to file
# -------------------------------------

### define data dimension.
n <- dim(dx)[1L]
p <- dim(dx)[2L]

### set up slab for this process
st <- c(rank * floor(n/ size) + 1, 1)
co <- dx@ldim

### set up dimensions for full matrix (not just local slab)
rdim <- ncdim_def("rows", "number", vals = 1:n)
cdim <- ncdim_def("columns", "number", vals = 1:p)

### define matrix variable in file (must not create full storage!!)
x.nc_var <- ncvar_def("testMatrix", "count",
                      list(rows = rdim, columns = cdim),
                      missval = -1, prec = "integer")

### create (collectively) in parallel a file with given dimensions
nc <- nc_create_par("./test_par.nc", x.nc_var)
nc_var_par_access(nc, x.nc_var)

### write variable values to file
ncvar_put(nc, "testMatrix", as.vector(submatrix(dx)), start = st, count = co)
nc_sync(nc)

### close file
nc_close(nc)

ncdump("./test_par.nc")

finalize()


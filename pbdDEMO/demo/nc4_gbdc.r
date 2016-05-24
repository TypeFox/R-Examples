suppressMessages(library(pbdDEMO, quietly = TRUE))
suppressMessages(library(pbdNCDF4, quietly = TRUE))

# -------------------------------------
# Write and read NetCDF4 file in gbd matrix
# -------------------------------------

### default of pbdMPI
rank <- comm.rank()
size <- comm.size()

### divide data into pieces by rank
X <- TREFHT$data
ncol <- ncol(X)
ncol.per.rank <- ceiling(ncol / size)
st <- 1 + ncol.per.rank * rank
en <- min(c(ncol, st + ncol.per.rank - 1))
X.gbdc <- X[, st:en]

# define dimension and variable
lon <- ncdim_def("lon", "degree_east", vals = TREFHT$def$dim[[1]]$vals)
lat <- ncdim_def("lat", "degree_north", vals = TREFHT$def$dim[[2]]$vals)
var.def <- ncvar_def("TREFHT", "K", list(lon = lon, lat = lat), NULL)

### parallel write
file.name <- "nc4_gbdc.nc"
nc <- nc_create_par(file.name, var.def)
demo.ncvar_put_gbd(nc, "TREFHT", X.gbdc, gbd.major = 2)
nc_close(nc)
if(comm.rank() == 0){
  ncdump(file.name)
}

### parallel read (everyone owns a portion)
nc <- nc_open_par(file.name)
if(comm.rank() == 0){
  print(nc)
}
new.X.gbdc <- demo.ncvar_get_gbd(nc, "TREFHT", gbd.major = 2)
nc_close(nc)

finalize()

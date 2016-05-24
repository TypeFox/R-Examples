library(pbdDEMO, quietly = TRUE)
library(pbdNCDF4, quietly = TRUE)

# -------------------------------------
# Serial write and read NCDF4 file
# -------------------------------------

### prepare data
X <- TREFHT$data

### define dimension and variable
lon <- ncdim_def("lon", "degree_east", vals = TREFHT$def$dim[[1]]$vals)
lat <- ncdim_def("lat", "degree_north", vals = TREFHT$def$dim[[2]]$vals)
var.def <- ncvar_def("TREFHT", "K", list(lon = lon, lat = lat), NULL)

### serial write
file.name <- "nc4_serial.nc"
if(comm.rank() == 0){
  nc <- nc_create(file.name, var.def)
  ncvar_put(nc, "TREFHT", X)
  nc_close(nc)
  ncdump(file.name)
}
barrier()

### serial read (everyone owns the same copy)
nc <- nc_open(file.name)
if(comm.rank() == 0){
  print(nc)
}
new.X <- ncvar_get(nc, "TREFHT")
nc_close(nc)

finalize()


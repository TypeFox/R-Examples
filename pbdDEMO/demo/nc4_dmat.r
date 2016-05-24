suppressMessages(library(pbdDEMO, quietly = TRUE))
suppressMessages(library(pbdNCDF4, quietly = TRUE))

# -------------------------------------
# Write and read NetCDF4 file in ddmatrix
# -------------------------------------

### initial grid
init.grid()
if(comm.size() != 4){
  comm.stop("This example requries 4 processors.")
}

### divide data into ddmatrix
x <- TREFHT$data
dx <- as.ddmatrix(x)

# define dimension and variable
lon <- ncdim_def("lon", "degree_east", vals = TREFHT$def$dim[[1]]$vals)
lat <- ncdim_def("lat", "degree_north", vals = TREFHT$def$dim[[2]]$vals)
var.def <- ncvar_def("TREFHT", "K", list(lon = lon, lat = lat), NULL)

### parallel write
file.name <- "nc4_dmat.nc"
nc <- nc_create_par(file.name, var.def)
demo.ncvar_put_dmat(nc, "TREFHT", dx)
nc_close(nc)
if(comm.rank() == 0){
  ncdump(file.name)
}

### parallel read (everyone owns a portion)
nc <- nc_open_par(file.name)
if(comm.rank() == 0){
  print(nc)
}
new.dx <- demo.ncvar_get_dmat(nc, "TREFHT", bldim = bldim(dx),
                              ICTXT = dmat.ictxt(dx))
nc_close(nc)

finalize()

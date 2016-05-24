suppressMessages(library(pbdDEMO, quietly = TRUE))
suppressMessages(library(pbdNCDF4, quietly = TRUE))

# -------------------------------------
# Read 2D NetCDF4 parallel file in ddmatrix
# -------------------------------------

### initial grid
init.grid()

### check directory
flag.file <- TRUE
if(comm.rank() == 0){
  if(!file.exists("nc4_data")){
    flag.file <- FALSE
  }
}
if(!comm.all(flag.file)){
  comm.stop("nc4_data does not exist.")
}

### divide data into ddmatrix
for(i.row in 1:9){
  for(i.col in 1:9){
    nrow <- i.row
    ncol <- i.col
    file.name <- paste("./nc4_data/2d_", nrow, "_", ncol, ".nc", sep = "")

    ### parallel read
    nc <- nc_open_par(file.name)
    new.dx <- demo.ncvar_get_dmat(nc, "data")
    nc_close(nc)
    # print(new.dx)
    x <- as.matrix(new.dx)
    comm.print(x)
  }
}
finalize()

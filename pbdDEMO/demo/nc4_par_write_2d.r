suppressMessages(library(pbdDEMO, quietly = TRUE))
suppressMessages(library(pbdNCDF4, quietly = TRUE))

# -------------------------------------
# Write 2D NetCDF4 parallel file in ddmatrix
# -------------------------------------

### initial grid
init.grid()

### check directory
if(comm.rank() == 0){
  if(!file.exists("nc4_data")){
    dir.create("nc4_data")
  }
}
barrier()

### divide data into ddmatrix
for(i.row in 1:9){
  for(i.col in 1:9){
    nrow <- i.row
    ncol <- i.col
    x <- matrix(1:(nrow * ncol), nrow = nrow, ncol = ncol)
    dx <- as.ddmatrix(x)

    ### define variables
    rowname <- ncdim_def("rowname", "1", vals = seq(1, nrow))
    colname <- ncdim_def("colname", "1", vals = seq(1, ncol))
    var.def <- ncvar_def("data", "1",
                         list(rowname = rowname, colname = colname),
                         NULL)

    ### parallel write
    file.name <- paste("./nc4_data/2d_", nrow, "_", ncol, ".nc", sep = "")
    nc <- nc_create_par(file.name, var.def)
    demo.ncvar_put_dmat(nc, "data", dx,verbose = FALSE)
    nc_close(nc)
  }
}

finalize()

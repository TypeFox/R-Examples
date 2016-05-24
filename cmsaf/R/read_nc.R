read_nc <- function(var,infile){

  id <- nc_open(infile)
    data <- ncvar_get(id,var)
  nc_close(id)

  return(data)
}

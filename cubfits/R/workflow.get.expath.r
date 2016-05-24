### Get locations of example files.

get.expath <- function(file.name, path.root = "./ex_data/", pkg = "cubfits"){
  file.name <- paste(path.root, file.name, sep = "")
  file.path <- tools::file_path_as_absolute(
                 system.file(file.name, package = "cubfits"))
  file.path
} # get.expath()


upload_dico <- function(file){
  if(!grepl(".csv$", file) && !grepl(".txt$", file)){
    stop("Uploaded file must be a .csv or .txt file!")
  }
  file.copy(file,paste(.libPaths()[1],"x.ent/dico/",sep="/"))
}
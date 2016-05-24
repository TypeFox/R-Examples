dir.exists <- function(path){
  system(paste("test -d", path), intern=FALSE) == 0 ; ## 0 for exists!
}

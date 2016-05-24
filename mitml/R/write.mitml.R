write.mitml <- function(x, filename, drop=FALSE){
# write mitml class object to file

  if(drop){
    x <- x[!names(x)%in%c("par.burnin","par.imputation")]
  }

  save(x,file=filename)
  invisible()
}


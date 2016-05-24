toNA<-function(x,values=NA) {
 if(missing(x)) 
  stop("Usage: toNA(x,values=NA)\n\twhere value is one or more values to be set to NA")
 x[x %in% values] <- NA
 return(x)
}
doubles<-function(x) {
  # Remove duplicates lines of a matrix
  # 
  # Recursive function for the removal of duplicate lines of a table
  # x is a boolean table, n rows par k columns
  # xr is the reduced table
  #
  zi1<-array(x[1,],c(1,dim(x)[-1])) ## Extrait la 1re ligne du tableau
  zi2<-dotprod(zi1,t(x),g="&",f="==")  ## determining all positions of the tested line in x
  xr<-x[(c(1-zi2))*c(1:dim(x)[1]),,drop=FALSE] ## remove duplicates
  res<- if (dim(xr)[1] <1) x else rbind(x[1,1:dim(x)[2]],doubles(xr))## on cycle until no duplicates in xr
}
tabresulSingl<-function(x,remove = FALSE) {
  # Prepare a table of results limited to the singletons
  # 
  # Input: x is the list of two elements resulting of nzdsr 
  # $DempsterRule: 1 col of masses plus boolean matrix;
  # $con: measure of conflicy between beliefs
  #
  # x : result from function nzdsr
r<-tabresul(x, remove = remove)
r<-r$mbp
z2<-r[,c(1:(ncol(r)-4))]
  if (!is.null(dim(z2))) {
    r1<-r[apply(z2,1,sum)==1,]
  } else {
    z2<-matrix(r[,c(1:(ncol(x$DempsterRule)-1))],ncol=length(z2))
    r1<-r[apply(z2,1,sum)==1,]
    if (is.null(dim(r1))) {
      r1<-t(as.matrix(r1))
      }
  }
r1<-list(mbp=r1, Conflict=x$con)
return(r1)
}

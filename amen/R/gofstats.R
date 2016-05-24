#' Goodness of fit statistics
#' 
#' Goodness of fit statistics evaluating second and third-order dependence
#' patterns
#' 
#' 
#' @usage gofstats(Y)
#' @param Y a relational data matrix
#' @return a vector of gof statistics
#' @author Peter Hoff
#' @examples
#' 
#' data(YX_nrm) 
#' 
#' gofstats(YX_nrm$Y) 
#' 
#' 
#' @export gofstats
gofstats<-function(Y)
{
  sd.rowmean<-sd(rowMeans(Y,na.rm=TRUE) ,na.rm=TRUE) 
  sd.colmean<-sd(colMeans(Y,na.rm=TRUE) ,na.rm=TRUE)
  
  dyad.dep<- suppressWarnings( cor( c(Y),c(t(Y)) , use="complete.obs") ) 
 
  E<-Y-mean(Y,na.rm=TRUE) ;  D<-1*(!is.na(E)) ; E[is.na(E)]<-0
  triad.dep<- sum(diag(E%*%E%*%E))/( sum(diag(D%*%D%*%D)) * sd(c(Y),na.rm=TRUE)^3)

  gof<-c(sd.rowmean,sd.colmean, dyad.dep , triad.dep ) 

  gof[is.na(gof)]<-0 

  names(gof)<-c("sd.rowmean","sd.colmean","dyad.dep","triad.dep")
  gof
}





#' Simulate an relational matrix based on a relative rank nomination scheme
#' 
#' Simulate an relational matrix based on a relative rank nomination scheme
#' 
#' 
#' @usage simY_rrl(EZ, rho, odobs, YO)
#' @param EZ a square matrix giving the expected value of the latent Z matrix
#' @param rho dyadic correlation
#' @param odobs a scalar or vector giving the observed number of nominations
#' for each node
#' @param YO a square matrix identifying where missing values should be
#' maintained
#' @return a square matrix, where higher values represent stronger
#' relationships
#' @author Peter Hoff
#' @export simY_rrl
simY_rrl <-
function(EZ,rho,odobs,YO=NULL)
{ 
  ZS<-simZ(EZ,rho)  
  diag(ZS)<- -Inf
  if(!is.null(YO)) { ZS[is.na(YO)]<- -Inf } 

  YS<-ZS*0  
  for(i in 1:nrow(EZ))
  {
    ri<-order( -ZS[i,] )[seq(1,odobs[i],length=odobs[i]) ]
    YS[i,ri]<- seq(odobs[i],1,length=odobs[i])
  }
  diag(YS)<-NA 
  if(!is.null(YO)) { YS[is.na(YO)]<- NA }
YS
}




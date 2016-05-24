#' Computes the design socioarray of covariate values
#' 
#' Computes the design socioarray of covariate values for an AME fit
#' 
#' 
#' @usage design_array(Xrow=NULL,Xcol=NULL,Xdyad=NULL,intercept=TRUE,n)
#' @param Xrow an n x pr matrix of row covariates
#' @param Xcol an n x pc matrix of column covariates
#' @param Xdyad an n x n x pd array of dyadic covariates
#' @param intercept logical
#' @param n number of rows/columns
#' @return an n x n x (pr+pc+pd+intercept) 3-way array
#' @author Peter Hoff
#' @export design_array
design_array<-function(Xrow=NULL,Xcol=NULL,Xdyad=NULL,intercept=TRUE,n)
{ 


### covariate array
pr<-length(Xrow)/n
pc<-length(Xcol)/n
pd<-length(Xdyad)/n^2
X<-array(dim=c(n,n,pr+pc+pd))
dnX<-NULL
###


### row covariates
if(pr>0)
{
  Xrow<-as.matrix(Xrow)
  Xrowa<-array(dim=c(n,n,pr))
  for( j in 1:pr ){ Xrowa[,,j]<-matrix( Xrow[,j], n,n) }
  X[,,1:pr]<- Xrowa
  dnX<-c(dnX,paste0(colnames(Xrow),rep(".row" ,pr)))
}
###


### column covariates
if(pc>0)
{
  Xcol<-as.matrix(Xcol)
  Xcola<-array(dim=c(n,n,pc))
  for( j in 1:pc ){ Xcola[,,j]<-t(matrix( Xcol[,j], n,n)) }
  X[,,pr+1:pc]<- Xcola
  dnX<-c(dnX,paste0(colnames(Xcol),rep(".col" ,pc)))
}
###


### dyadic covariates
if(pd>0)                                               
{  
  if(pd==1){ Xdyad<-array(Xdyad,dim=c(n,n,pd)) }                      
  X[,,pr+pc+1:pd]<-Xdyad                            
  dnX<-c(dnX,paste0(dimnames(Xdyad)[[3]],rep(".dyad",pd)))
}      
###


### add intercept 
if(!any(apply(X,3,function(x){var(c(x),na.rm=TRUE)})==0) & intercept )
{
  X1<-array(dim=c(0,0,1)+dim(X))
  X1[,,1]<-1 ; X1[,,-1]<-X
  X<-X1
  dnX<-c("intercept",dnX)
} 
###


### set variable names
if(dim(X)[[3]]>1) { dimnames(X)[[3]]<- dnX }
if(dim(X)[[3]]==1){ dimnames(X)[[3]]<- list(dnX) }
###


### missing values
if( sum(is.na(X)) > sum( is.na(apply(X,3,diag)) ) )
{
  cat("WARNING: replacing NAs in design matrix with zeros","\n")
} 
X[is.na(X)]<-0
###


X
}

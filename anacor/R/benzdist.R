benzdist <- function(scaling, x, y, z, tab, n, m, r, c, row.covariates, col.covariates)
{
# computes Benzecri distances for rows/columns and RMSE
  
 bdobs.col <- bdobs.row <- bdfit.col <- bdfit.row <- NULL   #initialize bd matrices  
 rmse.row <- rmse.col <- NULL
 if (any(scaling == "Benzecri"))
 {    

  #z <- as.matrix(tab/sqrt(outer(r,c)))
  if (is.matrix(row.covariates)) {                    #computation for rowcov
    xr<-cbind(rep(1,n),row.covariates)
    cx<-crossprod(xr,r*xr)
    nr<-dim(xr)[2]
    or<-symSqrt(cx,inv=TRUE)
    ROW<-TRUE
  } else ROW <- FALSE

  if (is.matrix(col.covariates)) {                    #computation for colcov
    yr<-cbind(rep(1,m),col.covariates)
    cy<-crossprod(yr,c*yr)
    mr<-dim(yr)[2]
    oc<-symSqrt(cy,inv=TRUE)         
    COL<-TRUE
  } else COL <- FALSE

  cc<-z%*%t(z)
  tau<-sum(tab)
 
  if (ROW) v<-xr%*%or%*%cc%*%or%*%t(xr) else v<-cc/sqrt(outer(r,r))   #n x n matrix of D^(-1/2)ZZ'D^(-1/2) (observed)
  
  #--------------------- row distances -------------------------
  if (scaling[1] == "Benzecri")
  {  
    w1<-diag(v)
    nn <- tau
    bdobs.row <- (outer(w1,w1,"+")-2*v)

    xrow <- x
    n<-dim(xrow)[1]
    c1 <- (xrow/nn)%*%(t(xrow)/nn)
    w2<-diag(c1)
    bdfit.row <-nn*(outer(w2,w2,"+")-2*c1)

    rmse.row <- sqrt((1/(n*(n-1)))*sum((bdobs.row - bdfit.row)^2))
  }

  if (scaling[2] == "Benzecri")
  {
    cc<-crossprod(z)
    if (COL) v<-yr%*%oc%*%cc%*%oc%*%t(yr) else v<-cc/sqrt(outer(c,c))  #m x m matrix
    nn <- tau
    w1<-diag(v)
    bdobs.col <- (outer(w1,w1,"+")-2*v)
  
    ycol <- y
    m<-dim(ycol)[1]
    c1<-(ycol/nn)%*%(t(ycol)/nn)
    w2<-diag(c1)
    bdfit.col <-nn*(outer(w2,w2,"+")-2*c1)
    
    rmse.col <- sqrt((1/(m*(m-1)))*sum((bdobs.col - bdfit.col)^2))
  }

  colnames(bdfit.col) <- colnames(bdobs.col)
  colnames(bdfit.row) <- colnames(bdobs.row)
  rownames(bdfit.col) <- rownames(bdobs.col)
  rownames(bdfit.row) <- rownames(bdobs.row)
  
 } 
  bdres <- list(bdobs.row = bdobs.row, bdfit.row = bdfit.row, bdobs.col = bdobs.col, bdfit.col = bdfit.col,
                rmse.row = rmse.row, rmse.col = rmse.col)    
}

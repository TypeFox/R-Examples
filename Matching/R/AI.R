#This function returns the homogeneous AI SE.
#This standalone function only works for ATT, and no weights
#The Match() function includes the version which works for everything
#  and which calculates heterogeneous variance estimates

AIse <- function(Y, Tr, index.treated, index.control, weights, est=NULL)
  {
    N <- length(Y)

    if(is.null(est))
      {
        est <- sum(Y[index.treated]*weights)/sum(weights)-
          sum(Y[index.control]*weights)/sum(weights)
      } else{
        est <- as.double(est)
      }

#    ret <- est.func(N=N, All=0, Tr=Tr,
#                    indx=cbind(index.treated,index.control,weights),
#                    weight=rep(1,N),
#                    BiasAdj=FALSE, Kz=NULL)

    ret <- .Call("EstFuncC", as.integer(N), as.integer(0), as.integer(length(index.treated)),
                 as.double(Y), as.double(Tr),
                 as.double(rep(1,N)), as.double(cbind(index.treated,index.control,weights)),
                 PACKAGE="Matching")

#    YCAUS <- ret$YCAUS
#    Kcount <- ret$Kcount
#    KKcount <- ret$KKcount

    YCAUS <- ret[,1]
    Kcount <- ret[,2]
    KKcount <- ret[,3]    

    Yt <- Y[index.treated]
    Yc <- Y[index.control]
    Tau <- Yt - Yc

    eps <- Tau - est
    eps.sq <- eps*eps
    Sigs <- 0.5 * matrix(1, N, 1) %*% (t(eps.sq) %*% weights)/sum(weights)

    SN <- sum(Tr)
    var.pop=sum((Sigs*((1-Tr)*Kcount*Kcount-(1-Tr)*KKcount))/(SN*SN))
    
    dvar.pop <- sum(Tr*(YCAUS-est)*(YCAUS-est))/(SN*SN)
    
    var <- var.pop + dvar.pop
    
    se <- sqrt(var)
    return(se)
  }

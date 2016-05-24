vcov.lgarch <-
function(object, arma=FALSE, ...)
{
  #vcovarma:
  if(is.null(object$vcov.arma)){
    cat("vcov cannot be extracted (vcov=FALSE during estimation)\n")
    m.dim <- length(object$par.arma)
    vcovarma <- matrix(NA,m.dim,m.dim)
  }else{
    vcovarma <- object$vcov.arma
    if(object$aux$mean.correction){
      vcovarma <- cbind(rep(NA,NCOL(vcovarma)), vcovarma)
      vcovarma <- rbind(rep(NA,NCOL(vcovarma)), vcovarma)
    }
  }
  colnames(vcovarma) <- names(object$par.arma)
  rownames(vcovarma) <- names(object$par.arma)

  #vcovlgarch:
  if(arma){
    result <- vcovarma
  }else{

    if(!is.null(object$vcov.lgarch)){
      result <- object$vcov.lgarch
    }else{

      #make vcovlgarch matrix:
      vcovlgarch <- matrix(NA,length(object$par),
        length(object$par))
      colnames(vcovlgarch) <- names(object$par)
      rownames(vcovlgarch) <- names(object$par)

      #if garch:
      if(object$aux$ma > 0){
        vcovlgarch[object$aux$ma.indx, object$aux$ma.indx] <- vcovarma[object$aux$ma.indx, object$aux$ma.indx]
      } #end if(..ma > 0)

      #if arch:
      if(object$aux$ar > 0){
        ar.ma.diff <- object$aux$ar - object$aux$ma
        if(ar.ma.diff==0){
          varalpha <- vcovarma[2,2] + vcovarma[3,3] + 2*vcovarma[2,3]
          vcovlgarch[2,2] <- varalpha
          covalphabeta <- -vcovarma[2,3] - vcovarma[3,3]
          vcovlgarch[2,3] <- covalphabeta
          vcovlgarch[3,2] <- covalphabeta
        }else{
          vcovlgarch[2,2] <- vcovarma[2,2]
        } #end if(..diff==0)
      } #end if(..ar > 0)

      #if xreg:
      if(object$aux$xreg.k > 0){
        vcovlgarch[object$aux$xreg.indx, object$aux$xreg.indx] <- vcovarma[object$aux$xreg.indx, object$aux$xreg.indx]
        if(object$aux$ma > 0){
          vcovlgarch[object$aux$ma.indx, object$aux$xreg.indx] <- -vcovarma[object$aux$ma.indx,object$aux$xreg.indx]
          vcovlgarch[object$aux$xreg.indx, object$aux$ma.indx] <- -vcovarma[object$aux$xreg.indx, object$aux$ma.indx]
        }
        if(object$aux$ar > 0){
          vcovlgarch[object$aux$ar.indx, object$aux$xreg.indx] <- vcovarma[object$aux$ar.indx,object$aux$xreg.indx]
          if(object$aux$ma > 0){
            vcovlgarch[object$aux$ar.indx, object$aux$xreg.indx] <- vcovlgarch[object$aux$ar.indx, object$aux$xreg.indx] + vcovlgarch[object$aux$ma.indx, object$aux$xreg.indx]
          }
          vcovlgarch[object$aux$xreg.indx, object$aux$ar.indx] <- vcovlgarch[object$aux$ar.indx, object$aux$xreg.indx]
        } #end if(..ar > 0)
      } #end if(..k > 0)

      #Var(Elnz2^hat):
      if(object$aux$method=="cex2"){
        vcovlgarch[NROW(vcovlgarch), NCOL(vcovlgarch)] <- vcovarma[NROW(vcovarma), NCOL(vcovarma)]
      }else{

        zhat <- coredata(residuals.lgarch(object, verbose=TRUE))
        if(object$aux$yzeron > 0){
          zhat <- zhat[-object$aux$yzerowhere]
        }  #end if..
        zhat2 <- zhat^2
        avar <- var(zhat2 - log(zhat2))
        vcovlgarch[NROW(vcovlgarch), NCOL(vcovlgarch)] <- avar/length(zhat)

#        old:
#        uadj <- lgarchRecursion1(as.numeric(object$par.arma), object$aux)
#        if(object$aux$yzeron > 0){
#          uadj <- uadj[-object$aux$yzerowhere]
#        }  #end if..
#        expuadj <- exp(uadj)
#        uexpuadj <- uadj*exp(uadj)
#        avaruadj <- var(uadj) + var(expuadj)/mean(expuadj)^2 - 2*mean(uexpuadj)/mean(expuadj)
#        vcovlgarch[NROW(vcovlgarch), NCOL(vcovlgarch)] <- avaruadj/length(uadj)

      } #end if("cex2")else(..)

      #TO DO:
      #intercept:
      #if(object$aux$method=="cex2"){
      #} #end if("cex2")else(..)

      result <- vcovlgarch

    } #end if(!is.null(vcovlgarch))else(..)
  } #end if(arma)else(..)

  #out:
  return(result)
}

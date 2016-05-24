vcov.mlgarch <-
function(object, varma=FALSE, ...)
{
  #if varma:
  if(varma){
    if(is.null(object$vcov.varma)){
      cat("vcov cannot be extracted (vcov=FALSE during estimation)\n")
      m.dim <- length(object$par.varma)
      result <- matrix(NA,m.dim,m.dim)
      colnames(result) <- names(object$par.varma)
      rownames(result) <- names(object$par.varma)
    }else{ result <- object$vcov.varma }
  }

  #if mlgarch:
  if(varma==FALSE){
    if(!is.null(object$vcov.mlgarch)){
      result <- object$vcov.mlgarch
    }else{

      #extract vcov.varma:
      if(is.null(object$vcov.varma)){
        cat("vcov cannot be extracted (vcov=FALSE during estimation)\n")
        m.dim <- length(object$par)
        result <- matrix(NA,m.dim,m.dim)
        colnames(result) <- names(object$par)
        rownames(result) <- names(object$par)
      }else{
        vcov.varma <- object$vcov.varma

        #make vcov.mlgarch:
        vcov.mlgarch <- matrix(NA,length(object$par),length(object$par))
        colnames(vcov.mlgarch) <- names(object$par)
        rownames(vcov.mlgarch) <- names(object$par)

        #make garch part:
        if(object$aux$ma > 0){
        vcov.mlgarch[object$aux$ma.indx, object$aux$ma.indx] <- vcov.varma[object$aux$ma.indx, object$aux$ma.indx]
        } #end if(..ma > 0)

        #make arch part:
        if(object$aux$ar > 0){
          ar.ma.diff <- object$aux$ar - object$aux$ma
          if(ar.ma.diff==0){
            vcov.ar <- diag(vcov.varma[object$aux$ar.indx,object$aux$ar.indx])
            vcov.ma <- diag(vcov.varma[object$aux$ma.indx,object$aux$ma.indx])
            vcov.ar.ma <- diag(vcov.varma[object$aux$ar.indx,object$aux$ma.indx])
            vcov.tmp <- matrix(NA,length(vcov.ar),length(vcov.ar))
            diag(vcov.tmp) <- vcov.ar + vcov.ma + 2*vcov.ar.ma
            vcov.mlgarch[object$aux$ar.indx,object$aux$ar.indx] <- vcov.tmp
          }else{
            vcov.mlgarch[object$aux$ar.indx,object$aux$ar.indx] <- vcov.varma[object$aux$ar.indx,object$aux$ar.indx]
          } #end if(..diff==0)
        } #end if(..ar > 0)

        #make xreg part:
        if(object$aux$xreg.k > 0){
          vcov.mlgarch[object$aux$xreg.indx, object$aux$xreg.indx] <- vcov.varma[object$aux$xreg.indx, object$aux$xreg.indx]
#           #univariate code:
#         if(object$aux$ma > 0){
#           vcov.lgarch[object$aux$ma.indx, object$aux$xreg.indx] <- -vcov.arma[object$aux$ma.indx,object$aux$xreg.indx]
#           vcov.lgarch[object$aux$xreg.indx, object$aux$ma.indx] <- -vcov.arma[object$aux$xreg.indx, object$aux$ma.indx]
#         }
#         if(object$aux$ar > 0){
#           vcov.lgarch[object$aux$ar.indx, object$aux$xreg.indx] <- vcov.arma[object$aux$ar.indx,object$aux$xreg.indx]
#           if(object$aux$ma > 0){
#             vcov.lgarch[object$aux$ar.indx, object$aux$xreg.indx] <- vcov.lgarch[object$aux$ar.indx, object$aux$xreg.indx] + vcov.lgarch[object$aux$ma.indx, object$aux$xreg.indx]
#           }
#           vcov.lgarch[object$aux$xreg.indx, object$aux$ar.indx] <- vcov.lgarch[object$aux$ar.indx, object$aux$xreg.indx]
#         } #end if(..ar > 0)
        } #end if(..k > 0)

        #Var(Elnz2^hat):
        uadj <- mlgarchRecursion1(as.numeric(object$par.varma), object$aux)
        uadj.m <- NCOL(uadj)
        vcov.m <- NROW(vcov.mlgarch)
        for(i in 1:uadj.m){
          uadjtmp <- uadj[,i]
          if(object$aux$yzeron[i] > 0){
            uadjtmp <- uadjtmp[ -object$aux$yzerowhere[[i]] ]
          }
          expuadj <- exp(uadjtmp)
          uexpuadj <- uadjtmp*exp(uadjtmp)
          avaruadj <- var(expuadj)/mean(expuadj)^2 + var(uadjtmp) - 2*mean(uexpuadj)/mean(expuadj)
          dim1 <- vcov.m - uadj.m + i
          dim2 <- dim1
          vcov.mlgarch[dim1, dim2] <- avaruadj/length(uadjtmp)
        } #end for(i in..)

      result <- vcov.mlgarch
      } #end if(is.null(vcov.varma))else(..)
    } #end if(!is.null(vcov.mlgarch))else(..)
  } #end if(varma==FALSE)else(..)

  #out:
  return(result)
}

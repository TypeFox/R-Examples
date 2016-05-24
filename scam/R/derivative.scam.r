#################################################################
## function to get derivatives of the smooth model terms ....  ##
## (currently only univariate smooths)                         ##
## analytic derivatives for P-, SCOP-splines                   ##
## finite differencing using predict.gam() for others          ##                          ##
################################################################# 


derivative.scam <- function(object, smooth.number=1,deriv=1){
 ## object - fitted scam object
 ## smooth.number - number of the smooth model term which derivative is needed to be calculated
 ## deriv - either 1 if the 1st derivative is required, or 2 if the 2nd

 sn <- smooth.number
 if (length(object$smooth[[sn]]$term)!=1) stop("derivative.smooth() currently handles only 1D smooths")
 if (deriv!=1 & deriv!=2) 
       stop(paste("deriv can be either 1 or 2"))
 n <- length(object$y)
 q <- object$smooth[[sn]]$bs.dim ## smooth basis dimension
 first <- object$smooth[[sn]]$first.para
 last <- object$smooth[[sn]]$last.para
 beta.t <- object$coefficients.t[first:last]
 Vp <- object$Vp.t[first:last,first:last] 

 if (inherits(object$smooth[[sn]], c("mpi.smooth","mpd.smooth", "cv.smooth", "cx.smooth",   
      "mdcv.smooth","mdcx.smooth","micv.smooth","micx.smooth"))){
      Sig <- object$smooth[[sn]]$Sigma 
      if (deriv==1){ ## 1st order derivatives
             P <- diff(diag(q-1),difference=1)
             Xd <- object$smooth[[sn]]$Xdf1%*%P%*%Sig   
          } 
          else { ## 2nd order derivative
             P <- diff(diag(q-1),difference=2)
             Xd <- object$smooth[[sn]]$Xdf2%*%P%*%Sig   
          }
 } # else if (inherits(object$smooth[[sn]], "pspline.smooth")){ ## unconstrained P-splines
   #         if (deriv==1)  ## 1st order derivatives
   #                Xd <- object$smooth[[sn]]$Xdf1%*%diff(diag(q-1),difference=1)
   #         else ## 2nd order derivative
   #               Xd <- object$smooth[[sn]]$Xdf2%*%diff(diag(q-1),difference=2) 
   # }
   else{   ## for other unconstrained smooths via finite differencing using predict.gam() 
          xx <- object$model[object$smooth[[sn]]$term]   ## find the data
          if (object$smooth[[sn]]$by != "NA") {
                  by <- rep(1, n)
                  newd <- data.frame(x = xx, by = by)
                  names(newd) <- c(object$smooth[[sn]]$term, object$smooth[[sn]]$by)
                }
                else {
                  newd <- data.frame(x = xx)
                  names(newd) <- object$smooth[[sn]]$term
                }
          X0 <- PredictMat(object$smooth[[sn]], newd)
          eps <- 1e-7 ## finite difference interval
          xx <- xx + eps 
          if (object$smooth[[sn]]$by != "NA") {
                  newd <- data.frame(x = xx, by = by)
                  names(newd) <- c(object$smooth[[sn]]$term, object$smooth[[sn]]$by)
                }
                else {
                  newd <- data.frame(x = xx)
                  names(newd) <- object$smooth[[sn]]$term
                }
          X1 <- PredictMat(object$smooth[[sn]], newd)
          Xd <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives
          if (deriv==2){ ## for 2nd oder derivative
                  xx <- xx + eps 
                  if (object$smooth[[sn]]$by != "NA") {
                           newd <- data.frame(x = xx, by = by)
                           names(newd) <- c(object$smooth[[sn]]$term, object$smooth[[sn]]$by)
                        }
                        else {
                           newd <- data.frame(x = xx)
                           names(newd) <- object$smooth[[sn]]$term
                        }
                   X2 <- PredictMat(object$smooth[[sn]], newd)
                   Xd <- (X2-2*X1+X0)/eps^2
              }
      }
 df <- Xd%*%beta.t ## values of the derivatives 
 df.sd <- rowSums(Xd%*%Vp*Xd)^.5  ## standard errors of the derivative of smooth 
 list(d=df,se.d=df.sd)
}  ## end of derivative.scam





#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# iqLearn_optTx1 : Estimate the recommended optimal first-stage treatment      #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# mainObj : object of class iqLearnFS_ME created by a call to iqLearnFSM       #
#                                                                              #
# cmObj   : object of class iqLearnFS_C created by a call to iqLearnFSC        #
#                                                                              #
# sigObj  : object of class iqLearnFS_VHom or iqLearnFS_VHet created by a call #
#           to function iqLearnVar                                             #
#                                                                              #
# dens    : either "norm" or "nonpar;" density estimator for the conditional   #
#           density of the contrast function                                   #
#                                                                              #
# newdata : data.frame of covariate information for new patients               #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns a list :                                                           =#
#=  $values : numeric, value functions at each treatment option               =#
#=  $optTx  : optimal tx (+1 if pos >= neg; -1 if pos < neg)                  =#
#=                                                                            =#
#==============================================================================#
iqLearn_optTx1 <- function(mainObj, 
                           cmObj, 
                           sigObj, 
                           dens, 
                           newdata){

  #--------------------------------------------------------------------------#
  # Recast dens input to lower case.                                         #
  #--------------------------------------------------------------------------#
  dens <- tolower(dens)

  #--------------------------------------------------------------------------#
  # Retrieve treatment information from contrast object                      #
  #--------------------------------------------------------------------------#
  tx <- TxVec(cmObj)

  if( !missing(newdata) ) {
    #----------------------------------------------------------------------#
    # If new data is provided obtain predicted values for the main effect  #
    # and contrast regressions at each treatment option.                   #
    #----------------------------------------------------------------------#
    lhat <- iqLearn_pm(object = mainObj, newdata = newdata)
    mu <- iqLearn_pm(object = cmObj, newdata = newdata)

    #----------------------------------------------------------------------#
    # Retrieve the variance information.                                   #
    #----------------------------------------------------------------------#
    if( is(sigObj,'IQLearnFS_VHom') ){
      sig <- StdDev(sigObj)
    } else if( is(sigObj,'IQLearnFS_VHet') ){
      sig <- iqLearn_pm(object = sigObj, newdata = newdata)
      sig <- exp(sig)*exp(Scale(sigObj)/2.0)
    }

  } else {
    #----------------------------------------------------------------------#
    # If no new data is provided retrieve fitted values.                   #
    #----------------------------------------------------------------------#
    lhat <- mainObj@qFunctions
    mu   <- cmObj@qFunctions
    sig  <- sigObj@qFunctions
    if( is(sigObj,'IQLearnFS_VHet') ){
      sig <- exp(sig)*exp(Scale(sigObj)/2.0)
    }

  }

  #--------------------------------------------------------------------------#
  # Estimate Q1 using either normal density or empirical estimate            #
  #--------------------------------------------------------------------------#
  q1Hat <- matrix(data = 0.0,  
                  nrow = nrow(lhat),  
                  ncol = 2L,  
                  dimnames=list(NULL, c("-1","1")))

  if( dens=="norm" ) {

    q1Hat[,2L] <- mu[,2L]*(1.0 - 2.0*pnorm(-mu[,2L]/sig[,2L])) +
                  sqrt(2.0/pi)*sig[,2L]*exp(-mu[,2L]^2/(2.0*sig[,2L]^2)) 
    q1Hat[,1L] <- mu[,1L]*(1.0 - 2.0*pnorm(-mu[,1L]/sig[,1L])) +
                  sqrt(2.0/pi)*sig[,1L]*exp(-mu[,1L]^2/(2.0*sig[,1L]^2)) 

    q1Hat <- lhat + q1Hat

  } else {
    func <- function(x,y,r,t,a){
              tma <- as.numeric(t==a)
              sum(abs(x + y*r)*tma)/sum(tma)
            }

    q1Hat[,2L] <- mapply(func, mu[,2L], sig[,2L], 
                  MoreArgs=list(r=Residuals(sigObj),t=tx,a=1L))
    q1Hat[,1L] <- mapply(func, mu[,1L], sig[,1L], 
                  MoreArgs=list(r=Residuals(sigObj),t=tx,a=-1L))

    q1Hat <- lhat + q1Hat

  }

  colnames(q1Hat) <- c("-1","1")

  #--------------------------------------------------------------------------#
  # Optimal tx is that with the largest value function                       #
  #--------------------------------------------------------------------------#
  optTx <- max.col(q1Hat)
  optTx <- as.integer(colnames(q1Hat)[optTx])
  
  #--------------------------------------------------------------------------#
  # Return q-functions and optimal treatment                                 #
  #--------------------------------------------------------------------------#
  return( list("qFunctions" = q1Hat, 
                "optimalTx" = optTx) )
}

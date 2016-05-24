"print.glmmNPML" <-
function(x,digits=max(3,getOption("digits")-3), ...){
  np <- length(x$coefficients)
  if (is.na(x$coef[length(x$coef)])){np<-np-1} ## Sep 2014
  # print(names(x))
  if (np > 0){   
    m <- seq(1,np)[substr(attr(x$coefficients,"names"),1,4)=="MASS"]
    mass.points <- x$coefficients[m]
    cat("\nCall: ",deparse(x$call),"\n\n")
    cat("Coefficients")
    cat(":\n")
    print.default(format(x$coefficients[1:np],digits = digits), print.gap = 2,quote = FALSE);cat("\n")
  } else {
    cat("\nCall: ",deparse(x$call),"\n\n")
    cat("No coefficients. \n")
  }
  
  if (x$family$family=='gaussian' && x$Misc$lambda<=1/length(x$masses)){ # print sigma only if it is constant over components
        cat("Component distribution - MLE of sigma:\t  ",                # 03/09/07
        format(signif(x$sdev$sdev,digits)),'\n')
    }
  if (x$family$family=='Gamma'&& x$Misc$lambda<=1/length(x$masses)){ # print shape only if it is constant over components
       cat("Component distribution - MLE of shape parameter:\t  ",format(signif(x$shape$shape,digits)),'\n')# 03/09/07
    }

  cat("Random effect distribution - standard deviation:\t  ", format(x$rsdev),"\n\n") # 03/09/07
  
  if (!is.null(x$post.prob)){
      p <- x$masses
      # names(p) <- paste('MASS',seq(1,ncol(x$post.prob)),sep='') # omitted from 0.42 
      cat('Mixture proportions')
      cat(":\n")
      print.default(format(p,digits),print.gap=2,quote=FALSE)
   }
  cat("-2 log L:\t   ", format(round(x$disparity,digits=1)),"\n")
  invisible(x)
}


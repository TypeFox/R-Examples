
#'@export
#'@rdname PPall
#'@param x an object of class \code{gpcm4pl} which is the result of using the \code{PPall()} function
#'@param ... just some points.
#'@method print ppeo
print.ppeo <-
  function(x, ...)
  {
  print(apply(x$resPP$resPP,2,function(qq) round(qq,4)))
  }

    

#'@export
#'@rdname PPall
#'@param object An object of class \code{gpcm4pl} which is the result of using the \code{PPall()} function
#'@param nrowmax When printing the matrix of estimates - how many rows should be shown? Default = 15.
#'@method summary ppeo
summary.ppeo <-
  function(object, nrowmax=15, ...)
  {
  cat("PP Version: ",as.character(attr(object$call,"version")),"\n")
  cat("\n Call:",deparse(object$call),"\n- job started @",attr(object$call,"date"),"\n\n")  
  
  cat("Estimation type:",object$type,"\n\n")
  
  cat("Number of iterations:",  object$resPP$nsteps,"\n")
  
  cat("-------------------------------------\n")
  
  respps <- apply(object$resPP$resPP,2,function(x) round(x,4))
  if(nrow(respps) <= nrowmax)
    {
      print(respps)
    } else
      {
        print(respps[1:nrowmax,])
        cat("--------> output truncated <--------\n")
      }
  
  }


#'@export
#'@rdname Jkpp
#'@param x an object of class \code{jk} which is the result of using the \code{JKpp()} function
#'@method print jk
print.jk <-
  function(x, ...)
  {
    print(apply(x$resjk,2,function(qq) round(qq,4)))
  }






#'@export
#'@rdname Jkpp
#'@param object An object of class \code{jk} which is the result of using the \code{JKpp()} function
#'@param nrowmax When printing the matrix of estimates - how many rows should be shown? Default = 15.
#'@method summary jk
summary.jk <-
  function(object, nrowmax=15, ...)
  {
    cat("PP Version: ",as.character(attr(object$call,"version")),"\n")
    cat("\n Call:",deparse(object$call),"\n- job started @",attr(object$call,"date"),"\n\n")  
    
    cat("Estimation type: Jackknife --> ",object$type,"\n\n")
    
    
    cat("-------------------------------------\n")
    
    respps <- apply(object$resjk,2,function(x) round(x,4))
    if(nrow(respps) <= nrowmax)
    {
      print(respps)
    } else
    {
      print(respps[1:nrowmax,])
      cat("--------> output truncated <--------\n")
    }
    
  }








#'@export
#'@rdname PV
#'@param x An object of class \code{pv} which is the result of using the \code{PV()} function
#'@method print pv
print.pv <-
  function(x, ...)
  {
    print(apply(x$pvdraws,2,function(qq) round(qq,4)))
  }






#'@export
#'@rdname PV
#'@param object An object of class \code{pv} which is the result of using the \code{PV()} function
#'@param nrowmax When printing the matrix of estimates - how many rows should be shown? Default = 15.
#'@method summary pv
summary.pv <-
  function(object, nrowmax=15, ...)
  {
    cat("PP Version: ",as.character(attr(object$call,"version")),"\n")
    cat("\n Call:",deparse(object$call),"\n- job started @",attr(object$call,"date"),"\n\n")  
    
    cat("Estimation type: Plausible values \n\n")
    
    
    cat("-------------------------------------\n")
    
    respps <- apply(object$pvdraws,2,function(x) round(x,4))
    if(nrow(respps) <= nrowmax)
    {
      print(respps)
    } else
    {
      print(respps[1:nrowmax,])
      cat("--------> output truncated <--------\n")
    }
    
  }























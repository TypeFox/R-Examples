##
##  PURPOSE:   Drawing traceplots of selected parameters from fitted objects
##             * method for objects of class GLMM_MCMC
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   16/01/2010
##
##  FUNCTIONS: tracePlots.GLMM_MCMC (16/01/2010)
##
## ======================================================================

## *************************************************************
## tracePlots.GLMM_MCMC
## *************************************************************
tracePlots.GLMM_MCMC <- function(x, param = c("Deviance", "Cond.Deviance", "alpha", "Eb", "SDb", "Corb", "sigma_eps", "w_b", "mu_b", "sd_b", "gammaInv_b", "gammaInv_eps"),
                                 relabel = FALSE, order,                                 
                                 auto.layout = TRUE, xlab = "Iteration", ylab, col = "slateblue", main = "", ...)
{
  param <- match.arg(param)

  ### Determine component in x where to look for required chains
  if (param == "Deviance") obj <- "Deviance"
  else if (param == "Cond.Deviance") obj <- "Cond.Deviance"
       else if (param %in% c("Eb", "SDb", "Corb")) obj <- "mixture_b"
            else if (param == "sd_b") obj <- "Sigma_b"
                 else                 obj <- param

  if (param %in% c("w_b", "mu_b", "sd_b") & x$prior.b$priorK != "fixed") stop("Not implemented for this value of param.")
  
  ### Number of parameters to plot
  if (param %in% c("Deviance", "Cond.Deviance")) nparam <- 1
  else if (param == "alpha") nparam <- sum(x$p) + sum(x$fixed.intercept)
       else if (param %in% c("Eb", "SDb")) nparam <- x$dimb
            else if (param == "Corb") nparam <- (x$dimb * (x$dimb + 1)) / 2 - x$dimb
                 else if (param == "sigma_eps") nparam <- x$R["Rc"]
                      else if (param == "w_b") nparam <- x$K_b[1]
                           else if (param %in% c("mu_b", "sd_b")) nparam <- x$K_b[1] * x$dimb
                                else if (param == "gammaInv_b") nparam <- x$dimb
                                     else if (param == "gammaInv_eps") nparam <- x$R["Rc"]
  if (!nparam){
    cat("\nNothing to plot.\n")
    return(invisible(x))    
  }  
  
  ### Plotting arguments
  if (length(col) == 1)  col <- rep(col, nparam)
  if (length(xlab) == 1) xlab <- rep(xlab, nparam)
  if (length(main) == 1) main <- rep(main, nparam)  
  if (length(col) != nparam)  stop("col must be of length", nparam)
  if (length(xlab) != nparam) stop("xlab must be of length", nparam)
  if (length(main) != nparam) stop("main must be of length", nparam)    

  if (!missing(ylab)){
    if (length(ylab) == 1)  ylab <- rep(ylab, nparam)
    if (length(ylab) != nparam) stop("ylab must be of length", nparam)    
  }
  
  ### Layout
  if (auto.layout){
    oldPar <- par(bty="n")
    layout(autolayout(nparam))
    on.exit(oldPar)
  }  

  ### Iteration index
  itIndex <- (x$nMCMC["burn"] + 1):(x$nMCMC["burn"] + x$nMCMC["keep"])

  ### Traceplots if related to moments of the mixture
  if (param %in% c("Eb", "SDb", "Corb")){
    if (param == "Eb") COLS <- paste("b.Mean.", 1:nparam, sep="")
    else if (param == "SDb") COLS <- paste("b.SD.", 1:nparam, sep="")
         else if (param == "Corb"){
           Imat <- diag(x$dimb)
           rowsI <- row(Imat)[lower.tri(row(Imat), diag=FALSE)]
           colsI <- col(Imat)[lower.tri(col(Imat), diag=FALSE)] 
           COLS <- paste("b.Corr.", rowsI, ".", colsI, sep="")           
         }
    
    if (missing(ylab)) ylab <- COLS
    
    for (i in 1:nparam){
      plot(itIndex, x[[obj]][, COLS[i]], type="l", xlab=xlab[i], ylab=ylab[i], col=col[i], main=main[i], ...)
    }  
  }

  else{

    ### Traceplots of mixture weights, means or standard deviations (possibly re-labeled)
    if (param %in% c("w_b", "mu_b", "sd_b")){
      if (!relabel) order <- matrix(rep(1:x$K_b[1], nrow(x[[obj]])), ncol=x$K_b[1], byrow=TRUE)
      if (relabel & missing(order)) order <- x$order_b
      if (nrow(order) != nrow(x[[obj]])) stop("order has incompatible number of rows")
      if (ncol(order) != x$K_b[1])       stop("order has incompatible number of columns")
      if (any(apply(order, 1, min) != 1))        stop("order must contain 1 in each row")
      if (any(apply(order, 1, max) != x$K_b[1])) stop("order must contain", x$K_b[1], "in each row")

      if (param == "w_b"){
        if (missing(ylab)) ylab <- colnames(x[[obj]])

        for (k in 1:x$K_b[1]){
          plot(itIndex, x[[obj]][cbind(1:nrow(x[[obj]]), order[,k])],
               type="l", xlab=xlab[k], ylab=ylab[k], col=col[k], main=main[k], ...)
        }
      }
      else{
        if (param == "mu_b"){     ### Draw traceplots of shifted and scaled mixture means
          if (missing(ylab)) ylab <- colnames(x[[obj]])

          for (k in 1:x$K_b[1]){
            for (j in 1:x$dimb){
              i <- (k-1)*x$dimb + j
              plot(itIndex, x$scale.b$shift[j] + x$scale.b$scale[j]*x[[obj]][cbind(1:nrow(x[[obj]]), (order[,k]-1)*x$dimb + j)],
                   type="l", xlab=xlab[i], ylab=ylab[i], col=col[i], main=main[i], ...)
            }
          }
        }
        else{
          if (param == "sd_b"){   ### Draw traceplots of scaled mixture standard deviations
            if (missing(ylab)) ylab <- paste("sd.", rep(1:x$K_b[1], each=x$dimb), ".", rep(1:x$dimb, x$K_b[1]), sep="")
            LTp <- (x$dimb * (x$dimb + 1)) / 2
            Idiag <- matrix(0, nrow=x$dimb, ncol=x$dimb)
            Idiag[lower.tri(Idiag, diag=TRUE)] <- 1:LTp
            jdiag <- diag(Idiag)
        
            for (k in 1:x$K_b[1]){
              for (j in 1:x$dimb){
                i <- (k-1)*x$dimb + j
                plot(itIndex, x$scale.b$scale[j]*sqrt(x[[obj]][cbind(1:nrow(x[[obj]]), (order[,k]-1)*LTp + jdiag[j])]),
                     type="l", xlab=xlab[i], ylab=ylab[i], col=col[i], main=main[i], ...)
              }
            }        
          }  
        }  
      }              
    }

    else{

      ### Traceplots for deviances
      if (param %in% c("Deviance", "Cond.Deviance")){
        
        if (missing(ylab)) ylab <- ifelse(param == "Deviance", "Deviance", "Conditional Deviance")
        plot(itIndex, x[[obj]], type="l", xlab=xlab[1], ylab=ylab, col=col[1], main=main[1], ...)
      }

      ### Traceplots of all other parameters
      else{
      
        if (missing(ylab)) ylab <- colnames(x[[obj]])
    
        for (i in 1:nparam){
          plot(itIndex, x[[obj]][, i], type="l", xlab=xlab[i], ylab=ylab[i], col=col[i], main=main[i], ...)
        }
      }  
    }
  }  
      
  return(invisible(x))
}



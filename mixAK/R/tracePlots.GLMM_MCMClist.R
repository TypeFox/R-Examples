##
##  PURPOSE:   Drawing traceplots of selected parameters from fitted objects
##             * method for objects of class GLMM_MCMClist
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   09/07/2013
##
##  FUNCTIONS: tracePlots.GLMM_MCMClist (09/07/2013)
##
## ======================================================================

## *************************************************************
## tracePlots.GLMM_MCMClist
## *************************************************************
tracePlots.GLMM_MCMClist <- function(x, param=c("Deviance", "Cond.Deviance", "alpha", "Eb", "SDb", "Corb", "sigma_eps", "w_b", "mu_b", "sd_b", "gammaInv_b", "gammaInv_eps"),
                                     relabel=FALSE,                             
                                     auto.layout=TRUE, xlab="Iteration", ylab, col=c("blue3", "red3"), main="", ...)
{
  param <- match.arg(param)

  ### Determine component in x where to look for required chains
  if (param == "Deviance") obj <- "Deviance"
  else if (param == "Cond.Deviance") obj <- "Cond.Deviance"
       else if (param %in% c("Eb", "SDb", "Corb")) obj <- "mixture_b"
            else if (param == "sd_b") obj <- "Sigma_b"
                 else                 obj <- param

  if (param %in% c("w_b", "mu_b", "sd_b") & x[[1]]$prior.b$priorK != "fixed") stop("Not implemented for this value of param.")
  
  ### Number of parameters to plot
  if (param %in% c("Deviance", "Cond.Deviance")) nparam <- 1
  else if (param == "alpha") nparam <- sum(x[[1]]$p) + sum(x[[1]]$fixed.intercept)
       else if (param %in% c("Eb", "SDb")) nparam <- x[[1]]$dimb
            else if (param == "Corb") nparam <- (x[[1]]$dimb * (x[[1]]$dimb + 1)) / 2 - x[[1]]$dimb
                 else if (param == "sigma_eps") nparam <- x[[1]]$R["Rc"]
                      else if (param == "w_b") nparam <- x[[1]]$K_b[1]
                           else if (param %in% c("mu_b", "sd_b")) nparam <- x[[1]]$K_b[1] * x[[1]]$dimb
                                else if (param == "gammaInv_b") nparam <- x[[1]]$dimb
                                     else if (param == "gammaInv_eps") nparam <- x[[1]]$R["Rc"]
  if (!nparam){
    cat("\nNothing to plot.\n")
    return(invisible(x))    
  }  
  
  ### Plotting arguments
  if (length(col) != 2) stop("col must be of length 2")
  if (length(xlab) == 1) xlab <- rep(xlab, nparam)
  if (length(main) == 1) main <- rep(main, nparam)  
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
  itIndex <- (x[[1]]$nMCMC["burn"] + 1):(x[[1]]$nMCMC["burn"] + x[[1]]$nMCMC["keep"])

  ### Traceplots if related to moments of the mixture
  if (param %in% c("Eb", "SDb", "Corb")){
    if (param == "Eb") COLS <- paste("b.Mean.", 1:nparam, sep="")
    else if (param == "SDb") COLS <- paste("b.SD.", 1:nparam, sep="")
         else if (param == "Corb"){
           Imat <- diag(x[[1]]$dimb)
           rowsI <- row(Imat)[lower.tri(row(Imat), diag=FALSE)]
           colsI <- col(Imat)[lower.tri(col(Imat), diag=FALSE)] 
           COLS <- paste("b.Corr.", rowsI, ".", colsI, sep="")           
         }
    
    if (missing(ylab)) ylab <- COLS
    
    for (i in 1:nparam){
      YLIM <- range(c(x[[1]][[obj]][, COLS[i]], x[[2]][[obj]][, COLS[i]]), na.rm = TRUE)
      plot(itIndex, x[[1]][[obj]][, COLS[i]], type="l", xlab=xlab[i], ylab=ylab[i], col=col[1], main=main[i], ylim=YLIM, ...)
      lines(itIndex, x[[2]][[obj]][, COLS[i]], col=col[2])
    }  
  }

  else{

    ### Traceplots of mixture weights, means or standard deviations (possibly re-labeled)
    if (param %in% c("w_b", "mu_b", "sd_b")){
      if (relabel){
        order1 <- x[[1]]$order_b
        order2 <- x[[2]]$order_b        
      }else{
        order1 <- matrix(rep(1:x[[1]]$K_b[1], nrow(x[[1]][[obj]])), ncol=x[[1]]$K_b[1], byrow=TRUE)
        order2 <- matrix(rep(1:x[[2]]$K_b[1], nrow(x[[2]][[obj]])), ncol=x[[2]]$K_b[1], byrow=TRUE)        
      }  

      if (param == "w_b"){
        if (missing(ylab)) ylab <- colnames(x[[1]][[obj]])

        for (k in 1:x[[1]]$K_b[1]){
          YLIM <- range(c(x[[1]][[obj]][cbind(1:nrow(x[[1]][[obj]]), order1[,k])], x[[2]][[obj]][cbind(1:nrow(x[[2]][[obj]]), order2[,k])]), na.rm = TRUE)
          plot(itIndex, x[[1]][[obj]][cbind(1:nrow(x[[1]][[obj]]), order1[,k])],
               type="l", xlab=xlab[k], ylab=ylab[k], col=col[1], main=main[k], ylim=YLIM, ...)
          lines(itIndex, x[[2]][[obj]][cbind(1:nrow(x[[2]][[obj]]), order2[,k])], col=col[2])
        }
      }
      else{
        if (param == "mu_b"){     ### Draw traceplots of shifted and scaled mixture means
          if (missing(ylab)) ylab <- colnames(x[[1]][[obj]])

          for (k in 1:x[[1]]$K_b[1]){
            for (j in 1:x[[1]]$dimb){
              i <- (k-1)*x[[1]]$dimb + j
              YLIM <- range(c(x[[1]]$scale.b$shift[j] + x[[1]]$scale.b$scale[j]*x[[1]][[obj]][cbind(1:nrow(x[[1]][[obj]]), (order1[,k]-1)*x[[1]]$dimb + j)],
                              x[[2]]$scale.b$shift[j] + x[[2]]$scale.b$scale[j]*x[[2]][[obj]][cbind(1:nrow(x[[2]][[obj]]), (order2[,k]-1)*x[[2]]$dimb + j)]), na.rm = TRUE)
              plot(itIndex, x[[1]]$scale.b$shift[j] + x[[1]]$scale.b$scale[j]*x[[1]][[obj]][cbind(1:nrow(x[[1]][[obj]]), (order1[,k]-1)*x[[1]]$dimb + j)],
                   type="l", xlab=xlab[i], ylab=ylab[i], col=col[1], main=main[i], ylim=YLIM, ...)
              lines(itIndex, x[[2]]$scale.b$shift[j] + x[[2]]$scale.b$scale[j]*x[[2]][[obj]][cbind(1:nrow(x[[2]][[obj]]), (order2[,k]-1)*x[[2]]$dimb + j)], col=col[2])
            }
          }
        }
        else{
          if (param == "sd_b"){   ### Draw traceplots of scaled mixture standard deviations
            if (missing(ylab)) ylab <- paste("sd.", rep(1:x[[1]]$K_b[1], each=x[[1]]$dimb), ".", rep(1:x[[1]]$dimb, x[[1]]$K_b[1]), sep="")
            LTp <- (x[[1]]$dimb * (x[[1]]$dimb + 1)) / 2
            Idiag <- matrix(0, nrow=x[[1]]$dimb, ncol=x[[1]]$dimb)
            Idiag[lower.tri(Idiag, diag=TRUE)] <- 1:LTp
            jdiag <- diag(Idiag)
        
            for (k in 1:x[[1]]$K_b[1]){
              for (j in 1:x[[1]]$dimb){
                i <- (k-1)*x[[1]]$dimb + j
                YLIM <- range(c(x[[1]]$scale.b$scale[j]*sqrt(x[[1]][[obj]][cbind(1:nrow(x[[1]][[obj]]), (order1[,k]-1)*LTp + jdiag[j])]),
                                x[[2]]$scale.b$scale[j]*sqrt(x[[2]][[obj]][cbind(1:nrow(x[[2]][[obj]]), (order2[,k]-1)*LTp + jdiag[j])])), na.rm = TRUE)
                plot(itIndex, x[[1]]$scale.b$scale[j]*sqrt(x[[1]][[obj]][cbind(1:nrow(x[[1]][[obj]]), (order1[,k]-1)*LTp + jdiag[j])]),
                     type="l", xlab=xlab[i], ylab=ylab[i], col=col[1], main=main[i], ylim=YLIM, ...)
                lines(itIndex, x[[2]]$scale.b$scale[j]*sqrt(x[[2]][[obj]][cbind(1:nrow(x[[2]][[obj]]), (order2[,k]-1)*LTp + jdiag[j])]), col=col[2])
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
        YLIM <- range(c(x[[1]][[obj]], x[[2]][[obj]]), na.rm = TRUE)
        plot(itIndex, x[[1]][[obj]], type="l", xlab=xlab[1], ylab=ylab, col=col[1], main=main[1], ylim=YLIM, ...)
        lines(itIndex, x[[2]][[obj]], col=col[2])
      }

      ### Traceplots of all other parameters
      else{
      
        if (missing(ylab)) ylab <- colnames(x[[1]][[obj]])
    
        for (i in 1:nparam){
          YLIM <- range(c(x[[1]][[obj]][, i], x[[2]][[obj]][, i]), na.rm = TRUE)
          plot(itIndex, x[[1]][[obj]][, i], type="l", xlab=xlab[i], ylab=ylab[i], col=col[1], main=main[i], ylim=YLIM, ...)
          lines(itIndex, x[[2]][[obj]][, i], col=col[2])
        }
      }  
    }
  }  
      
  return(invisible(x))
}  

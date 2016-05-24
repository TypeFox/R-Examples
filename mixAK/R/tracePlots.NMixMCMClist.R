##
##  PURPOSE:   Drawing traceplots of selected parameters from fitted objects
##             * method for objects of class NMixMCMClist
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   01/04/2015
##
##  FUNCTIONS: tracePlots.NMixMCMClist (01/04/2015)
##
## ======================================================================

## *************************************************************
## tracePlots.NMixMCMClist
## *************************************************************
tracePlots.NMixMCMClist <- function(x, param = c("Emix", "SDmix", "Cormix", "K", "w", "mu", "sd", "gammaInv"),
                                    relabel = FALSE,
                                    auto.layout = TRUE, xlab = "Iteration", ylab, col = c("blue3", "red3"), main = "", ...)
{
  param <- match.arg(param)

  ### Determine component in x where to look for required chains  
  if (param %in% c("Emix", "SDmix", "Cormix")) obj <- "mixture"
  else if (param == "sd") obj <- "Sigma"
       else               obj <- param

  if (param %in% c("w", "mu", "sd") & x[[1]]$prior$priorK != "fixed") stop("Not implemented for this value of param.")
  
  ### Number of parameters to plot
  if (param %in% c("Emix", "SDmix")) nparam <- x[[1]]$dim * x[[1]]$nx_w
  else if (param == "Cormix") nparam <- ((x[[1]]$dim * (x[[1]]$dim + 1)) / 2 - x[[1]]$dim) * x[[1]]$nx_w
       else if (param == "K") nparam <- 1
            else if (param == "w") nparam <- x[[1]]$K[1] * x[[1]]$nx_w
                 else if (param %in% c("mu", "sd")) nparam <- x[[1]]$K[1] * x[[1]]$dim
                      else if (param == "gammaInv") nparam <- x[[1]]$dim
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
  if (param %in% c("Emix", "SDmix", "Cormix")){
    if (param == "Emix"){
      if (x[[1]]$nx_w == 1){                
        COLS <- paste("y.Mean.", 1:nparam, sep="")
      }else{
        nparam0 <- nparam / x[[1]]$nx_w
        COLS <- paste(rep(paste("y.Mean.", 1:nparam0, sep=""), x[[1]]$nx_w), "-", rep(x[[1]]$lx_w, each = nparam0), sep = "")
      }    
    }else{
      if (param == "SDmix"){
        if (x[[1]]$nx_w == 1){          
          COLS <- paste("y.SD.", 1:nparam, sep="")
        }else{
          nparam0 <- nparam / x[[1]]$nx_w
          COLS <- paste(rep(paste("y.SD.", 1:nparam0, sep=""), x[[1]]$nx_w), "-", rep(x[[1]]$lx_w, each = nparam0), sep = "")
        }    
      }else{
        if (param == "Cormix"){
          Imat <- diag(x[[1]]$dim)
          rowsI <- row(Imat)[lower.tri(row(Imat), diag=FALSE)]
          colsI <- col(Imat)[lower.tri(col(Imat), diag=FALSE)]
          if (x[[1]]$nx_w == 1){
            COLS <- paste("y.Corr.", rowsI, ".", colsI, sep="")
          }else{
            nparam0 <- nparam / x[[1]]$nx_w              
            COLS <- paste(rep(paste("y.Corr.", rowsI, ".", colsI, sep=""), x[[1]]$nx_w), "-", rep(x[[1]]$lx_w, each = nparam0), sep = "")
          }    
        }
      }
    }
    if (missing(ylab)) ylab <- COLS

    for (i in 1:nparam){
      YLIM <- range(c(x[[1]][[obj]][, COLS[i]], x[[2]][[obj]][, COLS[i]]), na.rm = TRUE)      
      plot(itIndex, x[[1]][[obj]][, COLS[i]], type = "l", xlab = xlab[i], ylab = ylab[i], col = col[1], main = main[i], ylim=YLIM, ...)
      lines(itIndex, x[[2]][[obj]][, COLS[i]], col=col[2])      
    }  
  }

  else{

    ### Traceplots of mixture weights, means or standard deviations (possibly re-labeled)
    if (param %in% c("w", "mu", "sd")){###%%%
      if (relabel){
        order1 <- x[[1]]$order
        order2 <- x[[2]]$order
      }else{
        order1 <- matrix(rep(1:x[[1]]$K[1], nrow(x[[1]][[obj]])), ncol=x[[1]]$K[1], byrow=TRUE)
        order2 <- matrix(rep(1:x[[2]]$K[1], nrow(x[[2]][[obj]])), ncol=x[[2]]$K[1], byrow=TRUE)        
      }    

      if (param == "w"){
        if (missing(ylab)) ylab <- colnames(x[[1]][[obj]])

        if (x[[1]]$nx_w == 1){
          for (k in 1:x[[1]]$K[1]){
            YLIM <- range(c(x[[1]][[obj]][cbind(1:nrow(x[[1]][[obj]]), order1[,k])], x[[2]][[obj]][cbind(1:nrow(x[[2]][[obj]]), order2[,k])]), na.rm = TRUE)              
            plot(itIndex, x[[1]][[obj]][cbind(1:nrow(x[[obj]]), order1[,k])],
                 type = "l", xlab = xlab[k], ylab = ylab[k], col = col[1], main = main[k], ...)
            lines(itIndex, x[[2]][[obj]][cbind(1:nrow(x[[2]][[obj]]), order2[,k])], col=col[2])            
          }
        }else{
          for (ixw in 1:x[[1]]$nx_w){
            wixw1 <- x[[1]]$w[, (ixw-1)*x[[1]]$K[1] + (1:x[[1]]$K[1])]
            wixw2 <- x[[2]]$w[, (ixw-1)*x[[2]]$K[1] + (1:x[[2]]$K[1])]            
            for (k in 1:x[[1]]$K[1]){
              YLIM <- range(c(wixw1[cbind(1:nrow(wixw1), order1[,k])], wixw2[cbind(1:nrow(wixw2), order2[,k])]), na.rm = TRUE)              
              plot(itIndex, wixw1[cbind(1:nrow(wixw1), order1[,k])],
                   type = "l", xlab = xlab[(ixw - 1) * x[[1]]$K[1] + k], ylab = ylab[(ixw - 1) * x[[1]]$K[1] + k], col = col[1], main = main[(ixw - 1) * x[[1]]$K[1] + k], ylim = YLIM, ...)
              lines(itIndex, wixw2[cbind(1:nrow(wixw2), order2[,k])], col = col[2])
            }            
          }
        }    
      }
      else{
        if (param == "mu"){     ### Draw traceplots of shifted and scaled mixture means
          if (missing(ylab)) ylab <- colnames(x[[1]][[obj]])

          for (k in 1:x[[1]]$K[1]){
            for (j in 1:x[[1]]$dim){
              i <- (k-1)*x[[1]]$dim + j
              YLIM <- range(c(x[[1]]$scale$shift[j] + x[[1]]$scale$scale[j]*x[[1]][[obj]][cbind(1:nrow(x[[1]][[obj]]), (order1[,k]-1)*x[[1]]$dim + j)],
                              x[[2]]$scale$shift[j] + x[[2]]$scale$scale[j]*x[[2]][[obj]][cbind(1:nrow(x[[2]][[obj]]), (order2[,k]-1)*x[[2]]$dim + j)]), na.rm = TRUE)              
              plot(itIndex, x[[1]]$scale$shift[j] + x[[1]]$scale$scale[j]*x[[1]][[obj]][cbind(1:nrow(x[[1]][[obj]]), (order1[,k]-1)*x[[1]]$dim + j)],
                   type = "l", xlab = xlab[i], ylab = ylab[i], col = col[1], main = main[i], ylim = YLIM, ...)
              lines(itIndex, x[[2]]$scale$shift[j] + x[[2]]$scale$scale[j]*x[[2]][[obj]][cbind(1:nrow(x[[2]][[obj]]), (order2[,k]-1)*x[[2]]$dim + j)], col=col[2])              
            }
          }
        }
        else{
          if (param == "sd"){   ### Draw traceplots of scaled mixture standard deviations
            if (missing(ylab)) ylab <- paste("sd.", rep(1:x[[1]]$K[1], each=x[[1]]$dim), ".", rep(1:x[[1]]$dim, x[[1]]$K[1]), sep="")
            LTp <- (x[[1]]$dim * (x[[1]]$dim + 1)) / 2
            Idiag <- matrix(0, nrow=x[[1]]$dim, ncol=x[[1]]$dim)
            Idiag[lower.tri(Idiag, diag=TRUE)] <- 1:LTp
            jdiag <- diag(Idiag)
        
            for (k in 1:x[[1]]$K[1]){
              for (j in 1:x[[1]]$dim){
                i <- (k-1)*x[[1]]$dim + j
                YLIM <- range(c(x[[1]]$scale$scale[j]*sqrt(x[[1]][[obj]][cbind(1:nrow(x[[1]][[obj]]), (order1[,k]-1)*LTp + jdiag[j])]),
                                x[[2]]$scale$scale[j]*sqrt(x[[2]][[obj]][cbind(1:nrow(x[[2]][[obj]]), (order2[,k]-1)*LTp + jdiag[j])])), na.rm = TRUE)                
                plot(itIndex, x[[1]]$scale$scale[j]*sqrt(x[[1]][[obj]][cbind(1:nrow(x[[1]][[obj]]), (order1[,k]-1)*LTp + jdiag[j])]),
                     type="l", xlab = xlab[i], ylab = ylab[i], col = col[1], main = main[i], ylim = YLIM, ...)
                lines(itIndex, x[[2]]$scale$scale[j]*sqrt(x[[2]][[obj]][cbind(1:nrow(x[[2]][[obj]]), (order2[,k]-1)*LTp + jdiag[j])]), col=col[2])
              }
            }        
          }  
        }  
      }  
    }

    else{

      ### Traceplot for the number of mixture components
      if (param == "K"){
        if (missing(ylab)) ylab <- "K"
        YLIM <- range(c(x[[1]][[obj]], x[[1]][[obj]]), na.rm = TRUE)
        plot(itIndex, x[[1]][[obj]], type="l", xlab=xlab[1], ylab=ylab[1], col=col[1], main=main[1], ylim = YLIM, ...)
        lines(itIndex, x[[2]][[obj]], col = col[2])
      }
      
      ### Traceplots of all other parameters
      else{
        if (missing(ylab)) ylab <- colnames(x[[1]][[obj]])
      
        for (i in 1:nparam){
          YLIM <- range(c(x[[1]][[obj]][, i], x[[2]][[obj]][, i]), na.rm = TRUE)            
          plot(itIndex, x[[1]][[obj]][, i], type="l", xlab=xlab[i], ylab=ylab[i], col=col[1], main=main[i], ylim = YLIM, ...)
          lines(itIndex, x[[2]][[obj]][, i], col=col[2])          
        }
      }      
    }  
  }    
  
  return(invisible(x))
}

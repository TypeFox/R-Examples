##
##  PURPOSE:   Drawing traceplots of selected parameters from fitted objects
##             * method for objects of class NMixMCMC
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   16/01/2010
##             01/04/2015:  a factor covariate on mixture weights allowed
##
##  FUNCTIONS: tracePlots.NMixMCMC (16/01/2010)
##
## ======================================================================

## *************************************************************
## tracePlots.NMixMCMC
## *************************************************************
tracePlots.NMixMCMC <- function(x, param = c("Emix", "SDmix", "Cormix", "K", "w", "mu", "sd", "gammaInv"),
                                relabel = FALSE, order,
                                auto.layout = TRUE, xlab = "Iteration", ylab, col = "slateblue", main = "", ...)
{
  param <- match.arg(param)

  ### Determine component in x where to look for required chains  
  if (param %in% c("Emix", "SDmix", "Cormix")) obj <- "mixture"
  else if (param == "sd") obj <- "Sigma"
       else               obj <- param

  if (param %in% c("w", "mu", "sd") & x$prior$priorK != "fixed") stop("Not implemented for this value of param.")
  
  ### Number of parameters to plot
  if (param %in% c("Emix", "SDmix")) nparam <- x$dim * x$nx_w
  else if (param == "Cormix") nparam <- ((x$dim * (x$dim + 1)) / 2 - x$dim) * x$nx_w
       else if (param == "K") nparam <- 1
            else if (param == "w") nparam <- x$K[1] * x$nx_w
                 else if (param %in% c("mu", "sd")) nparam <- x$K[1] * x$dim
                      else if (param == "gammaInv") nparam <- x$dim
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
  if (param %in% c("Emix", "SDmix", "Cormix")){
    if (param == "Emix"){
      if (x$nx_w == 1){                
        COLS <- paste("y.Mean.", 1:nparam, sep="")
      }else{
        nparam0 <- nparam / x$nx_w
        COLS <- paste(rep(paste("y.Mean.", 1:nparam0, sep=""), x$nx_w), "-", rep(x$lx_w, each = nparam0), sep = "")
      }    
    }else{
      if (param == "SDmix"){
        if (x$nx_w == 1){          
          COLS <- paste("y.SD.", 1:nparam, sep="")
        }else{
          nparam0 <- nparam / x$nx_w
          COLS <- paste(rep(paste("y.SD.", 1:nparam0, sep=""), x$nx_w), "-", rep(x$lx_w, each = nparam0), sep = "")
        }    
      }else{
        if (param == "Cormix"){
          Imat <- diag(x$dim)
          rowsI <- row(Imat)[lower.tri(row(Imat), diag=FALSE)]
          colsI <- col(Imat)[lower.tri(col(Imat), diag=FALSE)]
          if (x$nx_w == 1){
            COLS <- paste("y.Corr.", rowsI, ".", colsI, sep="")
          }else{
            nparam0 <- nparam / x$nx_w              
            COLS <- paste(rep(paste("y.Corr.", rowsI, ".", colsI, sep=""), x$nx_w), "-", rep(x$lx_w, each = nparam0), sep = "")
          }    
        }
      }
    }
    if (missing(ylab)) ylab <- COLS

    for (i in 1:nparam){
      plot(itIndex, x[[obj]][, COLS[i]], type = "l", xlab = xlab[i], ylab = ylab[i], col = col[i], main = main[i], ...)
    }  
  }

  else{

    ### Traceplots of mixture weights, means or standard deviations (possibly re-labeled)
    if (param %in% c("w", "mu", "sd")){###%%%
      if (!relabel) order <- matrix(rep(1:x$K[1], nrow(x[[obj]])), ncol = x$K[1], byrow = TRUE)
      if (relabel & missing(order)) order <- x$order
      if (nrow(order) != nrow(x[[obj]])) stop("order has incompatible number of rows")
      if (ncol(order) != x$K[1])         stop("order has incompatible number of columns")
      if (any(apply(order, 1, min) != 1))      stop("order must contain 1 in each row")
      if (any(apply(order, 1, max) != x$K[1])) stop("order must contain", x$K[1], "in each row")

      if (param == "w"){
        if (missing(ylab)) ylab <- colnames(x[[obj]])

        if (x$nx_w == 1){
          for (k in 1:x$K[1]){
            plot(itIndex, x[[obj]][cbind(1:nrow(x[[obj]]), order[,k])],
                 type = "l", xlab = xlab[k], ylab = ylab[k], col = col[k], main = main[k], ...)
          }
        }else{
          for (ixw in 1:x$nx_w){
            wixw <- x$w[, (ixw-1)*x$K[1] + (1:x$K[1])]
            for (k in 1:x$K[1]){
              plot(itIndex, wixw[cbind(1:nrow(wixw), order[,k])],
                   type = "l", xlab = xlab[(ixw - 1) * x$K[1] + k], ylab = ylab[(ixw - 1) * x$K[1] + k], col = col[(ixw - 1) * x$K[1] + k], main = main[(ixw - 1) * x$K[1] + k], ...)
            }            
          }
        }    
      }
      else{
        if (param == "mu"){     ### Draw traceplots of shifted and scaled mixture means
          if (missing(ylab)) ylab <- colnames(x[[obj]])

          for (k in 1:x$K[1]){
            for (j in 1:x$dim){
              i <- (k-1)*x$dim + j
              plot(itIndex, x$scale$shift[j] + x$scale$scale[j]*x[[obj]][cbind(1:nrow(x[[obj]]), (order[,k]-1)*x$dim + j)],
                   type = "l", xlab = xlab[i], ylab = ylab[i], col = col[i], main = main[i], ...)
            }
          }
        }
        else{
          if (param == "sd"){   ### Draw traceplots of scaled mixture standard deviations
            if (missing(ylab)) ylab <- paste("sd.", rep(1:x$K[1], each=x$dim), ".", rep(1:x$dim, x$K[1]), sep="")
            LTp <- (x$dim * (x$dim + 1)) / 2
            Idiag <- matrix(0, nrow=x$dim, ncol=x$dim)
            Idiag[lower.tri(Idiag, diag=TRUE)] <- 1:LTp
            jdiag <- diag(Idiag)
        
            for (k in 1:x$K[1]){
              for (j in 1:x$dim){
                i <- (k-1)*x$dim + j
                plot(itIndex, x$scale$scale[j]*sqrt(x[[obj]][cbind(1:nrow(x[[obj]]), (order[,k]-1)*LTp + jdiag[j])]),
                     type="l", xlab = xlab[i], ylab = ylab[i], col = col[i], main = main[i], ...)
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
        plot(itIndex, x[[obj]], type="l", xlab=xlab[1], ylab=ylab[1], col=col[1], main=main[1], ...)
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

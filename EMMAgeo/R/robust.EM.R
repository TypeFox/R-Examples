#' Function to extract robust end-members.
#' 
#' This function takes a matrix with end-member loadings and extracts those
#' whose modes fall into specified limits. The function returns a list with all
#' passing end-member loadings and scores, along with their respective
#' coumn-wise (variable-wise) measures of centrality and dispersion.
#' 
#' 
#' @param Vqsn Numeric matrix with m samples (rows) and n variables (columns).
#' @param limits Numeric matrix with two columns that contain the boundaries of
#' mode classes for each end-member. The first column contains the lower, the
#' second column the upper limit. If \code{classunits} are provided, the limits
#' are assumed to relate to these units, if omitted column-numbers of
#' \code{Vqsn} are used.
#' @param quantiles Optional numeric vector of length two with the quantiles to
#' be evaluated for the robust end-member loadings; default is \code{c(0.25,
#' 0.75)}.
#' @param Vqn Numeric matrix with optional normalised factor loadings. If
#' present, the same factor loadings as the respectively selected end-member
#' loadings are returned.
#' @param classunits Numeric vector, optional class units (e.g. phi classes or
#' micrometers) of the same length as columns of X.
#' @param ID Numeric or character vector, optional sample IDs of the same
#' length as columns of X.
#' @param plot Logical scalar, optional graphical output of the results,
#' default is FALSE. If set to TRUE, selected end-member loadings are plotted
#' in different colours, according to the specified classes. All end-member
#' loadings are plotted in pale colour, means and standard deviations are
#' plotted above in thicker lines. To plot median and quantile range instead of
#' mean and standard deviation, add \code{median = TRUE} as further plot
#' parameter.  See examples section for further advice.
#' @param legend Character scalar, specifing legend position (cf.
#' \code{\link{legend}}). If omitted, no legend will be plotted, default is no
#' legend.
#' @param \dots Additional arguments passed to the plot function. Use
#' \code{colour} instead of \code{col} to create different colours.
#' @param pm Logical scalar to enable pm.
#' @return A list object containing: \item{Vqsn.data}{A list with Vqsn values.}
#' \item{Vqsn.mean}{A matrix with Vqsn means.} \item{Vqsn.median}{A matrix with
#' Vqsn medians.} \item{Vqsn.sd}{A matrix with Vqsn standard deviations.}
#' \item{Vqsn.qt1}{A matrix with Vqsn quantiles 1.} \item{Vqsn.qt2}{A matrix
#' with Vqsn quantiles 2.} \item{Vqn.data}{A list with Vqn values.}
#' \item{Vqn.mean}{A matrix with Vqn means.} \item{Vqn.median}{A matrix with
#' Vqn medians.} \item{Vqn.sd}{A matrix with Vqn standard deviations.}
#' \item{Vqn.qt1}{A matrix with Vqn quantiles 1.} \item{Vqn.qt2}{A matrix with
#' Vqn quantiles 2.}
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{EMMA}}, \code{\link{test.robustness}},
#' \code{\link{define.limits}}
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.
#' @keywords EMMA
#' @examples
#' 
#' ## load example data, i.e. here TR
#' data(TR, envir = environment())
#' 
#' ## define end-member limits
#' limits = cbind(c(11, 31, 60, 78), 
#'                c(13, 33, 62, 80))
#' 
#' ## extract robust end-members with limits matrix
#' rEM <- robust.EM(Vqsn = TR$Vqsn, limits = limits,
#'                  plot = TRUE,
#'                  legend = "topleft", 
#'                  cex = 0.7, 
#'                  colour = c("orange", "navyblue", "springgreen4", "red4"),
#'                  median = TRUE)
#' 
#' @export robust.EM
robust.EM <- function(
  Vqsn,
  limits,
  quantiles,
  Vqn,
  classunits,
  ID,
  plot = FALSE,
  legend,
  ...,
  pm = FALSE
) {
  
  ## check/set class units vector and test for consistency
  if(missing(classunits) == TRUE) {
    
    classunits <- seq(from = 1,
                      to = ncol(Vqsn))
  }
  
  if(missing(quantiles) == TRUE) {
    
    quantiles <- c(0.25, 0.75)
  }
  
  if(ncol(Vqsn) != length(classunits)) {
    
    stop("Units vector is not of same length as variables.")
  }
    
  ## check/set ID vector and test for consistency
  if(missing(ID) == TRUE) {
    
    ID <- seq(from = 1, 
              to = nrow(Vqsn))
  }
  
  if(nrow(Vqsn) != length(ID)) {
    stop("ID vector is not of same length as variables.")
  }
  
  ## create modes vector
  modes <- rep(x = NA, 
               times = nrow(Vqsn))
  
  ## determine mode class for all end-member loadings
  for(i in 1:nrow(Vqsn)) {
    
    modes[i] <- classunits[Vqsn[i,1:ncol(Vqsn)] == max(Vqsn[i,1:ncol(Vqsn)])]
  }
  
  ## create dummy list structures
  EM.Vqsn.list <- list(matrix(nrow = 1, 
                              ncol = ncol(Vqsn)))
  EM.Vqn.list  <- list(matrix(nrow = 1, 
                              ncol = ncol(Vqsn)))
  
  ## select modes that fall into limits for all limit pairs
  for(i in 1:nrow(limits)) {

    ## assign valid loadings
    EM.Vqsn <- Vqsn[(modes >= limits[i,1] & modes <= limits[i,2]),]
    
    ## append loadings matrix to list
    EM.Vqsn.list[[length(EM.Vqsn.list) + 1]] <- EM.Vqsn
    
    ## test if Vqn data set is present and if so assign valid loadings to list
    if(missing(Vqn) != TRUE) {
      
      EM.Vqn <- Vqn[(modes >= limits[i,1] & modes <= limits[i,2]),]
      EM.Vqn.list[[length(EM.Vqn.list) + 1]] <- EM.Vqn
    }
  }
  ## remove dummy matrices from list
  EM.Vqsn.list[1] <- NULL
  EM.Vqn.list[1]  <- NULL
  
  ## infer empty classes
  empty.classes <- rep(NA, nrow(limits))
  for(i in 1:nrow(limits)) {
    empty.classes[i] <- ifelse(test = length(EM.Vqsn.list[[i]]) == 0, 
                               yes = TRUE, 
                               no = FALSE)
  }
  
  ## check for complete data sets
  if(sum(empty.classes) == 0) {
    
    ## CASE 1 - complete results list
  
    ## create dummy output matrices for existing data
    EM.Vqsn.mean   <- matrix(nrow = nrow(limits), 
                             ncol = ncol(Vqsn))
    EM.Vqsn.median <- matrix(nrow = nrow(limits), 
                             ncol = ncol(Vqsn))
    EM.Vqsn.sd     <- matrix(nrow = nrow(limits), 
                             ncol = ncol(Vqsn))
    EM.Vqsn.qt1    <- matrix(nrow = nrow(limits), 
                             ncol = ncol(Vqsn))
    EM.Vqsn.qt2    <- matrix(nrow = nrow(limits), 
                             ncol = ncol(Vqsn))
    EM.Vqn.mean    <- matrix(nrow = nrow(limits), 
                             ncol = ncol(Vqsn))
    EM.Vqn.median  <- matrix(nrow = nrow(limits), 
                             ncol = ncol(Vqsn))
    EM.Vqn.sd      <- matrix(nrow = nrow(limits), 
                             ncol = ncol(Vqsn))
    EM.Vqn.qt1     <- matrix(nrow = nrow(limits), 
                             ncol = ncol(Vqsn))
    EM.Vqn.qt2     <- matrix(nrow = nrow(limits), 
                             ncol = ncol(Vqsn))
    
    ## define quantiles function
    qts <- function(X, quantiles) {quantile(X, quantiles)}
    
    ## calculate mean, median, sd and quantiles for Vqsn and, if present, Vqn
    for (i in 1:nrow(limits)) {
      EM.Vqsn.mean[i,]    <- apply(X = EM.Vqsn.list[[i]], 
                                   MARGIN = 2, 
                                   FUN = mean)
      EM.Vqsn.median[i,]  <- apply(X = EM.Vqsn.list[[i]], 
                                   MARGIN = 2, 
                                   FUN = median)
      EM.Vqsn.sd[i,]      <- apply(X = EM.Vqsn.list[[i]], 
                                   MARGIN = 2, 
                                   FUN = sd)
      EM.Vqsn.qt1[i,]     <- apply(X = EM.Vqsn.list[[i]], 
                                   MARGIN = 2, 
                                   FUN = qts, quantiles[1])
      EM.Vqsn.qt2[i,]     <- apply(X = EM.Vqsn.list[[i]], 
                                   MARGIN = 2, 
                                   FUN = qts, quantiles[2])

      if(missing(Vqn) != TRUE) {
        
        EM.Vqn.mean[i,]   <- apply(X = EM.Vqn.list[[i]], 
                                   MARGIN = 2, 
                                   FUN = mean)
        EM.Vqn.median[i,]  <- apply(X = EM.Vqn.list[[i]], 
                                    MARGIN = 2, 
                                    FUN = median)
        EM.Vqn.sd[i,]     <- apply(X = EM.Vqn.list[[i]], 
                                   MARGIN = 2, 
                                   FUN = sd)
        EM.Vqn.qt1[i,]     <- apply(X = EM.Vqn.list[[i]], 
                                    MARGIN = 2, 
                                    FUN = qts, quantiles[1])
        EM.Vqn.qt2[i,]     <- apply(X = EM.Vqn.list[[i]], 
                                    MARGIN = 2, 
                                    FUN = qts, quantiles[2])
      }
    }
    
    ## optionally plot results
    if(plot == TRUE) {
      
      ## adjust plot margins
      par(oma = c(0, 1, 0, 0))
      
      ## read additional arguments list and check/set default values
      extraArgs <- list(...)
      main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
      {expression(paste("End-member loadings (", V[qsn], ")", sep = ""))}
      xlab <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab} else
      {"Class"}
      ylab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else
      {"Amount, relative"}
      ylim <- if("ylim" %in% names(extraArgs)) {extraArgs$ylim} else
      {c(0, max(Vqsn, na.rm = TRUE))}
      colour <- if("colour" %in% names(extraArgs)) {extraArgs$colour} else
      {seq(1, nrow(EM.Vqsn.mean))}
      if("legend" %in% names(extraArgs)) {legend.text <- extraArgs$legend} else
      {legend.text <- rep(NA, nrow(EM.Vqsn.mean))
       for(i in 1:nrow(EM.Vqsn.mean)) {
         mode <- classunits[EM.Vqsn.mean[i,] == max(EM.Vqsn.mean[i,])]
         legend.text[i] <- paste("EM ", i, " (", round(mode, 2), ")", 
                                 sep = "")}}
      legend.cex <- if("cex" %in% names(extraArgs)) {extraArgs$cex} else
      {1}
      legend.lty <- if("lty" %in% names(extraArgs)) {extraArgs$lty} else
      {1}
      median <- if("median" %in% names(extraArgs)) {extraArgs$median} else
      {FALSE}
      
      ## create plot colour vector in hsv-space
      plot.colour <- t(rgb2hsv(col2rgb(colour)))
      plot.colour.half <- plot.colour
      plot.colour.half[,2] <- 0.1
      if(sum(col2rgb(colour[1]) == c(0,0,0)) == 3) {
        plot.colour.half[1,] <- c(0, 0, 0.74)}
      
      ## plot all end-member loadings coloured by limit class    
      ## plot first curve of first end-member class
      plot(classunits, EM.Vqsn.list[[1]][1,],
           type = "l", 
           main = main,
           xlab = xlab,
           ylab = ylab,
           ylim = ylim,
           col  = plot.colour[,1])
      ## plot other curves of first end-member class
      if(nrow(EM.Vqsn.list[[1]] > 1)) {for(i in 1:nrow(EM.Vqsn.list[[1]])) {
        lines(classunits, EM.Vqsn.list[[1]][i,], 
              col = hsv(h = plot.colour.half[1,1],
                        s = plot.colour.half[1,2],
                        v = plot.colour.half[1,3]))}}
      ## plot all curves of all other end-member classes
      if(nrow(limits) > 1) {for (j in 2: nrow(limits)) {
        lines(classunits, EM.Vqsn.list[[j]][1,], 
              type = "l", 
              col = hsv(h = plot.colour.half[j,1],
                        s = plot.colour.half[j,2],
                        v = plot.colour.half[j,3]))
        if(nrow(EM.Vqsn.list[[j]] > 1)) {for (i in 1:nrow(EM.Vqsn.list[[j]])) {
          lines(classunits, EM.Vqsn.list[[j]][i,], 
                col = hsv(h = plot.colour.half[j,1],
                          s = plot.colour.half[j,2],
                          v = plot.colour.half[j,3]))
        }}
      }}
      ## optionally, plot mean and standard deviation curves
      if(median == FALSE){
        lines(classunits, EM.Vqsn.mean[1,] - EM.Vqsn.sd[1,],
              col = hsv(h = plot.colour[1,1],
                        s = plot.colour[1,2],
                        v = plot.colour[1,3]))
        lines(classunits, EM.Vqsn.mean[1,], lwd = 2,
              col = hsv(h = plot.colour[1,1],
                        s = plot.colour[1,2],
                        v = plot.colour[1,3]))
        lines(classunits, EM.Vqsn.mean[1,] + EM.Vqsn.sd[1,],
              col = hsv(h = plot.colour[1,1],
                        s = plot.colour[1,2],
                        v = plot.colour[1,3]))
        if(nrow(limits) > 1) for(i in 2:nrow(limits)) {
          lines(classunits, EM.Vqsn.mean[i,] - EM.Vqsn.sd[i,],
                col = hsv(h = plot.colour[i,1],
                          s = plot.colour[i,2],
                          v = plot.colour[i,3]))
          lines(classunits, EM.Vqsn.mean[i,], lwd = 2,
                col = hsv(h = plot.colour[i,1],
                          s = plot.colour[i,2],
                          v = plot.colour[i,3]))
          lines(classunits, EM.Vqsn.mean[i,] + EM.Vqsn.sd[i,],
                col = hsv(h = plot.colour[i,1],
                          s = plot.colour[i,2],
                          v = plot.colour[i,3]))
        }
      } else {
        ## alternatively plot median and quantile range
        lines(classunits, EM.Vqsn.qt1[1,],
              col = hsv(h = plot.colour[1,1],
                        s = plot.colour[1,2],
                        v = plot.colour[1,3]))
        lines(classunits, EM.Vqsn.median[1,], lwd = 2,
              col = hsv(h = plot.colour[1,1],
                        s = plot.colour[1,2],
                        v = plot.colour[1,3]))
        lines(classunits, EM.Vqsn.qt2[1,],
              col = hsv(h = plot.colour[1,1],
                        s = plot.colour[1,2],
                        v = plot.colour[1,3]))
        if(nrow(limits) > 1) for(i in 2:nrow(limits)) {
          lines(classunits, EM.Vqsn.qt1[i,],
                col = hsv(h = plot.colour[i,1],
                          s = plot.colour[i,2],
                          v = plot.colour[i,3]))
          lines(classunits, EM.Vqsn.median[i,], lwd = 2,
                col = hsv(h = plot.colour[i,1],
                          s = plot.colour[i,2],
                          v = plot.colour[i,3]))
          lines(classunits, EM.Vqsn.qt2[i,],
                col = hsv(h = plot.colour[i,1],
                          s = plot.colour[i,2],
                          v = plot.colour[i,3]))
        }
      }
      ## optionally add legend
      if(missing(legend) == FALSE) {legend.position <- legend
                                    legend(x = legend.position,
                                           legend = legend.text,
                                           col = hsv(h = plot.colour[,1],
                                                     s = plot.colour[,2],
                                                     v = plot.colour[,3]),
                                           cex = legend.cex,
                                           lty = legend.lty)}
    }
    
    ## optionally add pm
    if(pm == TRUE) {pm <- check.data(matrix(runif(4), ncol = 2), 
                                     5, 0.01, 100, invisible = FALSE)}
    
    ## readjust plot margins
    par(oma = c(0, 0, 0, 0))
    
  } else {
    ## CASE 2 - incomplete results list
    
    ## create dummy output matrices for empty data
    EM.Vqsn.mean   <- c()
    EM.Vqsn.median <- c()
    EM.Vqsn.sd     <- c()
    EM.Vqsn.qt1    <- c()
    EM.Vqsn.qt2    <- c()
    EM.Vqn.mean    <- c()
    EM.Vqn.median  <- c()
    EM.Vqn.sd      <- c()
    EM.Vqn.qt1     <- c()
    EM.Vqn.qt2     <- c()
    
    ## notify empty output
    print(paste("No end-members found that fit to the following limits: ",
          limits[empty.classes,1], " - ", limits[empty.classes,2],
          ". No output generated. Remove empty class.",
          sep = ""))
  }
  
  ## return output, either with or without Vqn data
  if(missing(Vqn) != TRUE) {
    ##value<< A list object containing:
    list(Vqsn.data   = EM.Vqsn.list,   ##<< A list with Vqsn values.
         Vqsn.mean   = EM.Vqsn.mean,   ##<< A matrix with Vqsn means.
         Vqsn.median = EM.Vqsn.median, ##<< A matrix with Vqsn medians.
         Vqsn.sd     = EM.Vqsn.sd, ##<< A matrix with Vqsn standard deviations.
         Vqsn.qt1    = EM.Vqsn.qt1,    ##<< A matrix with Vqsn quantiles 1.
         Vqsn.qt2    = EM.Vqsn.qt2,    ##<< A matrix with Vqsn quantiles 2.
         Vqn.data    = EM.Vqn.list,    ##<< A list with Vqn values.
         Vqn.mean    = EM.Vqn.mean,    ##<< A matrix with Vqn means.
         Vqn.median = EM.Vqn.median, ##<< A matrix with Vqn medians.
         Vqn.sd      = EM.Vqn.sd,   ##<< A matrix with Vqn standard deviations.
         Vqn.qt1    = EM.Vqn.qt1,    ##<< A matrix with Vqn quantiles 1.
         Vqn.qt2    = EM.Vqn.qt2)    ##<< A matrix with Vqn quantiles 2.
    ##end<<
  } else {
    list(Vqsn.data   = EM.Vqsn.list,
         Vqsn.mean   = EM.Vqsn.mean,
         Vqsn.median = EM.Vqsn.median,
         Vqsn.sd     = EM.Vqsn.sd,
         Vqsn.qt1    = EM.Vqsn.qt1,
         Vqsn.qt2    = EM.Vqsn.qt2)
  }
}
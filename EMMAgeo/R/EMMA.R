#' Function for end-member modelling analysis.
#' 
#' A multivariate data set (m samples composed of n variables) is decomposed by
#' eigenspace analysis and modelled with a given number of end-members (q).
#' Several steps of scaling, transformation, normalisation, eigen space
#' decomposition, factor rotation, data modelling and evaluation are performed.
#' 
#' The function values \code{$loadings} and \code{$scores} are redundant. They
#' are essentially the same as \code{$Vqsn} and \code{$Mqs}. However, they are
#' included for user convenience. \cr\cr We kindly thank Christoph Burow for 
#' his quick contribution to remove unnecessary loops.
#' 
#' @param X Numeric matrix with m samples (rows) and n variables (columns).
#' @param q Numeric scalar with number of end-members to be modelled.
#' @param l Numeric scalar with the weight tranformation limit, i.e.
#' quantiles, cf. Klovan & Imbrie (1971); default is 0.
#' @param c Numeric scalar specifying the constant sum scaling parameter, e.g.
#' 1, 100, 1000; default is 100.
#' @param Vqn Numeric matrix specifying optional unscaled user-defined
#' end-member loadings. If provided, these are used instead of model-derived
#' ones.
#' @param EM.ID Character vector with end-member names. If present, these will
#' be set as row-names of the output data set and used in the legend text.
#' @param classunits Numeric vector, optional class units (e.g. micrometers or
#' phi-units) of the same length as columns of X.
#' @param ID Numeric or character vector, optional sample IDs of the same
#' length as rows of X.
#' @param rotation Character scalar, rotation type, default is "Varimax" (cf.
#' Dietze et al., 2012). One out of the rotations provided in GPArotation is
#' possible (cf. \code{\link{rotations}}).
#' @param plot Logical scalar, optional graphical output of the results,
#' default is FALSE. If set to TRUE, end-member loadings and end-member scores
#' are plotted.
#' @param \dots Additional arguments passed to the plot function. Since the
#' function returns two plots some additional graphical parameters must be
#' specified as vector with the first element for the first plot and the second
#' element for the second plot. 
#' @param pm Logical scalar to enable pm.
#' @return A list with numeric matrix objects. \item{loadings}{Normalised
#' rescaled end-member loadings.} \item{scores}{Rescaled end-member scores.}
#' \item{Vqn}{Normalised end-member loadings.} \item{Vqsn}{Normalised rescaled
#' end-member loadings.} \item{Mqs}{Rescaled end-member scores.}
#' \item{Xm}{Modelled data.} \item{modes}{Mode class of end-member loadings.}
#' \item{Mqs.var}{Explained variance of end-members} \item{Em}{Absolute
#' row-wise model error.} \item{En}{Absolute column-wise model error.}
#' \item{Rm}{Row-wise (sample-wise) explained variance.} \item{Rn}{Column-wise
#' (variable-wise) explained variance.} \item{ol}{Number of overlapping
#' end-members.}
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{test.parameters}}, \code{\link{rotations}},
#' \code{\link{eigen}}, \code{\link{nnls}}
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.\cr
#' Klovan JE, Imbrie J. 1971. An Algorithm and FORTRAN-IV Program for
#' Large-Scale Q-Mode Factor Analysis and Calculation of Factor Scores.
#' Mathematical Geology 3: 61-77. Miesch AT. 1976. Q-Mode factor analysis of
#' geochemical and petrologic data matrices with constant row sums. U.S.
#' Geological Survey Professsional Papers 574.
#' @keywords EMMA
#' @examples
#' 
#' ## load example data and set phi-vector
#' data(X, envir = environment())
#' phi <- seq(from = 1, to = 10, length.out = ncol(X))
#' 
#' ## perform EMMA with 5 end-members
#' EM <- EMMA(X = X, q = 5, l = 0.05, c = 100, plot = TRUE)
#' 
#' ## perform EMMA with 4 end-members and more graphical settings
#' EM <- EMMA(X = X, q = 4, l = 0.05, c = 100, 
#'            plot = TRUE,
#'            EM.ID = c("EM 1", "EM 2", "EM 3", "EM 4"),
#'            classunits = phi,
#'            xlab = c(expression(paste("Class [", phi, "]")), "Sample ID"),
#'            cex = 0.7,
#'            col = colors()[c(441, 496, 499, 506)])
#' 
#' @export EMMA
EMMA <- function(
 X,
 q, 
 l,
 c,
 Vqn,
 EM.ID,
 classunits,
 ID,
 rotation = "Varimax",
 plot = FALSE,
 ...,
 pm = FALSE
) {
  
  ## check for l vs. lw
  if("lw" %in% names(list(...))) {
    stop('Parameter "lw" is depreciated. Use "l" instead.')
  }
   
  ## check/set default values
  if(missing(l) == TRUE) {
    l <- 0
  }
  
  if(missing(c) == TRUE) {
    c <- 100
  }
  
  ## check/set class units vector and test for consistency
  if(missing(classunits) == TRUE) {
    classunits <- 1:ncol(X)
  }
  
  if(ncol(X) != length(classunits)) {
    stop("Units vector is not of same length as variables.")
  }
  
  ## check/set ID vector and test for consistency
  if(missing(ID) == TRUE) {
    ID <- 1:nrow(X)
  }
  
  if(nrow(X) != length(ID)) {
    stop("ID vector is not of same length as variables.")
  }
  
  ## End-member modelling
  ## rescale X to constant sum c
  X  <- X / apply(X = X, MARGIN = 1, FUN = sum) * c

  ## calculate weight limit quantiles column-wise
  ls <- apply(X <- X, 
              MARGIN = 2, 
              FUN = quantile, 
              probs = c(l, 1 - l), 
              type = 5)
  
  ## perform weight-transformation
  W <- t((t(X) - ls[1,]) / (ls[2,] - ls[1,]))

  ## create similarity matrix as outer product
  A <- t(W) %*% W

  ## perform eigenspace decomposition
  EIG <- eigen(A)

  ## assign raw eigenvectors V
  V <- EIG$vectors[,order(seq(ncol(A), 1, -1))]
  Vf <- V[,order(seq(ncol(A), 1, -1))]

  ## rotate eigenvectors and assign factor loadings
  Vr <- do.call(what = rotation, 
                args = list(Vf[,1:q]))
  Vq <- Vr$loadings[,seq(q, 1, -1)]

  ## rescale factor loadings and transpose matrix
  Vqr <- t(Vq) / apply(X = Vq, MARGIN = 2, FUN = sum) * c

  ## optionally calculate normalised factor loadings
  if(missing(Vqn) == TRUE) {
    Vqn <- ((Vqr - apply(X = Vqr, MARGIN = 1, FUN = min)) / (
      apply(X = Vqr, MARGIN = 1, FUN = max) - 
        apply(X = Vqr, MARGIN = 1, FUN = min)))
  } else {
    ## or use respective values from input data  
    Vqn <- Vqn
  }

  ## get eigenvector scores as non-negative least squares estimate
  Mq <- t(apply(W, MARGIN = 1, FUN = function(row) {
    limSolve::nnls(t(Vqn), as.vector(t(row)))$X
  }))
  
  ## modelled values matrix
  Wm <- Mq %*% Vqn

  ## calculate scaling vector after Miesch (1976)
  ls <- apply(X <- X, 
              MARGIN = 2, 
              FUN = quantile, 
              probs = c(l, 1 - l), 
              type = 5)
  
  s <- c - sum(ls[1,]) / apply(X = Vqn * (ls[2,] - ls[1,]), 
                               MARGIN = 1, 
                               FUN = sum)

  ## rescale end-member loadings after Miesch (1976)
  Vqs <- t(apply(cbind(t(t(s)), Vqn), MARGIN = 1, FUN = function(row) {
    row[1] * row[2:length(row)] * (ls[2, ] - ls[1, ]) + ls[1, ]
  }))

  ## normalise end-member loadings
  Vqsn <- Vqs / apply(Vqs, 1, sum) * c

  ## rescale end-member scores
  Mqs <- t(t(Mq) / s) / apply(t(t(Mq) / s), 1, sum)
  
  ## Model evaluation
  ## get number of overlapping end-member loadings
  ol <- sum(apply(Vqsn, MARGIN = 1, FUN = function(row) {
    if (row[row == max(row)] < max(row[row == max(row)])) { 1 } else { 0 }
  }))
  
  ## calculate modelled data set
  Xm <- Mqs %*% Vqs

  ## normalise modelled data set
  Xm  <- Xm / apply(X = Xm, MARGIN = 1, FUN = sum) * c
  
  ## evaluate absolute error and explained variances
  Mqs.var <- diag(var(Mqs)) / 
    sum(diag(var(Mqs))) * 100             # explained scores variance
  Em <- as.vector(apply(abs(X - Xm), 1, mean)) # absolute row-wise model error
  En <- as.vector(apply(abs(X - Xm), 2, mean)) # absolute column-wise model error
  Rm <- diag(cor(t(X), t(Xm))^2)          # row-wise explained variance
  Rn <- diag(cor(X, Xm)^2)                # column-wise explained variance
  RMSEm <- sqrt(mean(x = Em^2, na.rm = TRUE))
  RMSEn <- sqrt(mean(x = En^2, na.rm = TRUE))
  mRm <- mean(Rm, na.rm = TRUE)
  mRn <- mean(Rn, na.rm = TRUE)
  mRt <- mean(c(mRm, mRn), na.rm = TRUE)

  ## Sort output by mode position for consistent output
  ## create auxiliary x-unit and index variable, calculate index
  x.unit <- 1:ncol(X)
  ind <- numeric(q)
  for (i in 1:q) {ind[i] <- x.unit[Vqsn[i,] == max(Vqsn[i,])]}
  ind.sort <- seq(1, nrow(Vqsn))[order(ind)]
  
  Vqsn <- Vqsn[ind.sort,]
  Vqn <- Vqn[ind.sort,]
  Mqs <- Mqs[,ind.sort]
  Mqs.var <- Mqs.var[ind.sort]
  
  if(missing(EM.ID) == TRUE) {
    EM.ID <- paste("EM", 1:q, sep = "")
  }
  
  ## determine mode class for all end-member loadings
  modes <- numeric(q)
  for(i in 1:q) {modes[i] <- classunits[Vqsn[i,1:ncol(
    Vqsn)] == max(Vqsn[i,1:ncol(Vqsn)])]}
 
  ## optionally, plot end-member loadings and scores
  if(plot == TRUE) {
    
    ## read additional arguments list and check/set default values
    extraArgs <- list(...)
    
    if("main" %in% names(extraArgs)) {
      main <- extraArgs$main
    } else {
      main <- c("End-member loadings",
                "End-member scores")
    }
    
    if("xlab" %in% names(extraArgs)) {
      xlab <- extraArgs$xlab
    } else {
      xlab <- c("Class",
                "Samples")
    }
    
    if("ylab" %in% names(extraArgs)) {
      ylab <- extraArgs$ylab
    } else {
      ylab <- c("Relative amounts",
                "Relative amounts")
    }
    
    if("ylim" %in% names(extraArgs)) {
      ylim <- extraArgs$ylim
    } else {
      ylim <- rbind(c(0, max(Vqsn, na.rm = TRUE)),
                    c(0, 1))
    }
    
    if("log" %in% names(extraArgs)) {
      log <- extraArgs$log
    } else {
      log <- ""
    }

    if("col" %in% names(extraArgs)) {
      col <- extraArgs$col
    } else {
      col <- seq(1, q)
    }
    
    if("cex" %in% names(extraArgs)) {
      cex <- extraArgs$cex
    } else {
      cex <- 1
    }
    
    if("lty" %in% names(extraArgs)) {
      lty <- extraArgs$lty
    } else {
      lty <- 1
    }

    if("lwd" %in% names(extraArgs)) {
      lwd <- extraArgs$lwd
    } else {
      lwd <- 1
    }

    ## setup plot area
    ## save origial parameters
    par.old <- op <- par(no.readonly = TRUE)
    
    ## adjust margins
    par(mar = c(3.9, 4.5, 3.1, 1.1))
    
    ## define layout
    layout(rbind(c(1, 1, 2, 2), 
                 c(1, 1, 2, 2), 
                 c(3, 3, 4, 4), 
                 c(3, 3, 4, 4),
                 c(3, 3, 4, 4), 
                 c(3, 3, 4, 4),
                 c(5, 5, 5, 5))
           )

    ## plot col-wise explained variance
    plot(x = classunits, 
         y = Rn, 
         type = "b", 
         main = paste("Class-wise explained variance (", 
                      signif(x = mean(Rn * 100), digits = 2), " %)", sep = ""),
         xlab = xlab[1],
         ylab = expression(R^2),
         log = log)
    
    ## plot row-wise explained variance
    plot(x = 1:nrow(X), 
         y = Rm, 
         type = "b", 
         main = paste("Sample-wise explained variance (", 
                      signif(x = mean(Rm * 100), digits = 2), " %)", sep = ""),
         xlab = xlab[2],
         ylab = expression(R^2))
    
    ## plot end-member loadings
    plot(classunits, Vqsn[1,], type = "l", 
         main = main[1],
         xlab = xlab[1],
         ylab = ylab[1],
         ylim = as.vector(ylim[1,]),
         log = log,
         col = col[1])
    
    if(nrow(Vqsn) >= 2) {
      for(i in 2:nrow(Vqsn)) {
        lines(x = classunits, 
              y = Vqsn[i,], 
              col = col[i])
      }
    }
    
    ## plot end-member scores
    barplot(t(Mqs), names.arg = ID,
            main = main[2],
            xlab = xlab[2], 
            ylab = ylab[2],
            ylim = as.vector(ylim[2,]),
            col = col,
            horiz = FALSE)
    
    ## plot legend-like information
    par(mar = c(1, 1, 1, 1))
    
    plot(NA, 
         xlim = c(0, 1), 
         ylim = c(0, 1), 
         axes = FALSE, 
         ann = FALSE, 
         frame.plot = TRUE)
    
    mtext(line = -2, 
          text = "End-member ID (mode position | explained variance)", 
          cex = 0.9 * cex)
    
    legend(x = "bottom", 
           legend = paste(EM.ID, " (", signif(x = modes, digits = 2), " | ", 
                          signif(x = Mqs.var, digits = 2), " %)", sep = ""), 
           col = col, 
           lty = lty, 
           lwd = lwd, 
           horiz = TRUE, 
           box.lty = 0)
    
    ## reset plot parameters
    par(par.old)
  }

  ## optionally, assign EM.IDs to Vqsn and Vqn matrices
  if(missing(EM.ID) == FALSE) {
    rownames(Vqn) <- EM.ID
    rownames(Vqsn) <- EM.ID
  }
  
  ## Output list with numeric matrix objects.
  list(loadings = Vqsn,    ## Normalised rescaled end-member loadings.
       scores   = Mqs,     ## Rescaled end-member scores.
       Vqn      = Vqn,     ## Normalised end-member loadings.
       Vqsn     = Vqsn,    ## Normalised rescaled end-member loadings.
       Mqs      = Mqs,     ## Rescaled end-member scores.
       Xm       = Xm,      ## Modelled data.
       modes    = modes,   ## Mode class of end-member loadings.
       Mqs.var  = Mqs.var, ## Explained variance of end-members
       Em       = Em,      ## Absolute row-wise model error.
       En       = En,      ## Absolute column-wise model error.
       RMSEm    = RMSEm,   ## Sample-wise RMSE.
       RMSEn    = RMSEn,   ## Class-wise RMSE.
       Rm       = Rm,      ## Row-wise (sample-wise) explained variance.
       Rn       = Rn,      ## Column-wise (variable-wise) explained variance.
       mRm      = mRm,      ## Average class-wise explained variance.
       mRn      = mRn,      ## Average sample-wise average explained variance.
       mRt      = mRt,      ## Total average explained variance.
       ol       = ol)      ## Number of overlapping end-members.
}
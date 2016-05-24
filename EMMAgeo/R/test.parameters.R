#' Function to evaluate influence of model parameters.
#' 
#' All possible combinations of number of end-members and weight transformation
#' limits are used to perform EMMA. The function returns matrices of absolute
#' and relative measures of individual model performance.
#' 
#' The mean total explained variance mRt may be used to define a maximum number
#' of meaningful end-members for subsequent modelling, e.g. as the number of
#' end-members, which reaches the first local mRt maximum.\cr\cr Overlapping is
#' defined as one end-member having its mode within the "area" of any other
#' end-member, which is genetically not explainable.\cr\cr Special 
#' acknowledgements go to Christoph Burow for his efforts to implement the 
#' multicore functionality to this function.
#' 
#' @param X Numeric matrix with m samples (rows) and n variables (columns).
#' @param q Numeric vector of length two, specifying the minimum and maximum
#' number of end-members to be modelled.
#' @param l Numeric vector specifying the weight tranformation limit, i.e.
#' quantile; default is 0.
#' @param c Numeric scalar specifying the constant sum scaling parameter, e.g.
#' 1, 100, 1000; default is 0.
#' @param rotation Character scalar, rotation type, default is "Varimax" (cf.
#' Dietze et al., 2012). One out of the rotations provided in GPArotation is
#' possible (cf. \code{\link{rotations}}).
#' @param plot Character scalar, optional graphical output of the results.
#' Specify which tested parameter will be plotted: "mEm" (mean absolute
#' row-wise error), "mEn" (mean absolute column-wise error), "mRm" (mean
#' relative row-wise error), "mRn" (mean relative column-wise error), "mRt"
#' (mean relative total error), "ol" (number of overlapping end-members). All
#' plots except "ol" are colour-coded bitmaps of q, l and the specified test
#' parameter and line-plots the specified parameter vs. q.
#' @param legend Character scalar, specifying legend position (cf.
#' \code{\link{legend}}).  If omitted, no legend will be plotted, default is no
#' legend.
#' @param progressbar Logical scalar, optionally show a progress bar, default
#' is \code{FALSE}. Only available if option \code{multicore} is not used.
#' @param multicore Logical scalar, optionally ditribute calculations to all 
#' available cores of the computer, default is \code{TRUE}.
#' @param \dots Additional arguments passed to the plot function. Since the
#' function returns two plots (except for plot option "ol"), additional
#' graphical parameters must be specified as vector with the first element for
#' the first plot and the second element for the second plot. If graphical
#' parameters are natively vectors (e.g. a sequence of colours), they must be
#' specified as matrices with each vector as a row. A legend can only be added
#' to the second plot. Colours only apply to the second plot as well. If
#' colours are specified, \code{colour} should be used instead of \code{col}.
#' See example section for further advice.
#' @param pm Logical scalar to enable pm.
#' @return A list with result objects \item{mEm}{Absolute row-wise model
#' error.} \item{mEn}{Absolute column-wise model error.} \item{mRm}{Mean
#' row-wise explained variance.} \item{mRn}{Mean column-wise explained
#' variance.} \item{mRt}{Mean total explained variance.} \item{ol}{Number of
#' overlapping end-member loadings.} \item{q.max}{Maximum number of meaningful
#' end-members.}
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{EMMA}}
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set
#' data(X, envir = environment())
#' 
#' ## truncate the data set for faster computation
#' X.trunc <- X[1:20,]
#' 
#' ## define test parameters
#' q <- 2:8 # number of end-members
#' l <- seq(from = 0, to = 0.3, by = 0.1)
#' 
#' ## test parameter influence and plot mean total explained variance
#' TP <- test.parameters(X = X.trunc, q = q, l = l, plot = "mRt",
#'                       legend = "bottomright", cex = 0.7,
#'                       multicore = FALSE,
#'                       colour = rgb((1:7) / 7, 0.9, 0.2, 1))
#' 
#' ## show maximum number of end-members
#' TP$q.max
#' 
#' @export test.parameters
test.parameters <-function(
  X,
  q, 
  l = 0,
  c = 100,
  rotation = "Varimax",
  plot = FALSE,
  legend,
  progressbar = FALSE,
  multicore = FALSE,
  ...,
  pm = FALSE
){
  
  ## check for l vs. lw
  if("lw" %in% names(list(...))) {
    stop('Parameter "lw" is depreciated. Use "l" instead.')
  }
  
  ## definition of target vectors and matrices
  q.t <- rep(q, length(l))
  l.t <- rep(seq(min(l), max(l), length.out = length(l)), 
              each = q[length(q)] - 1)
  mEm <- rep(NA, length(q.t))
  mEn <- mEm
  mRm <- mEm
  mRn <- mEm
  mRt <- mEm
  ol  <- mEm
  i.pb <- 1
  
  ## optionally setup initial progress bar
  if(progressbar == TRUE & multicore == FALSE) {
    pb <- txtProgressBar(min = 0, max = length(l.t), char = "=", style = 3)
  }
  
  if(multicore == TRUE) {
    
    ## detect cores
    cores <- parallel::detectCores()
    
    ## initiate cluster
    cl <- parallel::makeCluster(getOption("mc.cores", cores))
    
    ## run EMMA in parallel
    EM <- parallel::clusterMap(cl, function(qi, li, X, c, rotation) {
      tryCatch(expr = {
        library(EMMAgeo)
        tmp <- EMMA(X = X, 
                    q = qi, 
                    l = li, 
                    c = c, 
                    rotation = rotation)
        c(mean(tmp$Em), 
          mean(tmp$En), 
          mean(tmp$Rm), 
          mean(tmp$Rn), 
          mean(c(tmp$Rm, tmp$Rn)), 
          tmp$ol)
      },
      error = function(e) { NA })
    },
    q.t, 
    l.t, 
    SIMPLIFY = TRUE, 
    MoreArgs = list(X = X, 
                    c = c, 
                    rotation = rotation))
    
    ## stop cluster
    parallel::stopCluster(cl) 
    
    # assign putput to target variables
    mEm <- EM[1, ]
    mEn <- EM[2, ]
    mRm <- EM[3, ]
    mRn <- EM[4, ]
    mRt <- EM[5, ]
    ol <- EM[6, ]
  } else {
    ## loop through all test q-vector elements
    for(i in 1:length(q.t)) {
      ## perform EMMA with respective parameter combinations
      EM <- try(EMMA(X, 
                     q = q.t[i], 
                     l = l.t[i], 
                     c = c,
                     rotation = rotation), silent = TRUE)
      
      ## test if EMMA resulted an error
      test.EM <- as.logical(grepl("Error", EM[1]))
      if (test.EM == TRUE) {
        ## if EMMA was erroneous, set result values to NA
        mEm[i] <- NA
        mEn[i] <- NA
        mRm[i] <- NA
        mRn[i] <- NA
        mRt[i] <- NA
        ol[i]  <- NA
        
        ## update progress bar
        if(progressbar == TRUE) {
          setTxtProgressBar(pb, i)
        }
        
      } else {
        ## if EMMA was successful, assign result values
        mEm[i] <- mean(abs(EM$Em)) # absolute row-wise model error vector
        mEn[i] <- mean(abs(EM$En)) # absolute column-wise model error vector
        mRm[i] <- mean(EM$Rm) # mean row-wise explained variance vector
        mRn[i] <- mean(EM$Rn) # mean row-wise explained variance vector
        mRt[i] <- mean(c(EM$Rm, EM$Rn)) # mean total explained variance vector
        ol[i]  <- EM$ol # mean total explained variance vector is NA
        
        if(progressbar == TRUE & multicore == FALSE) {
          setTxtProgressBar(pb, i) # update progress bar
        }
      }
    }
    
    ## close progress bar
    if(progressbar == TRUE & multicore == FALSE) {
      close(pb)
    }
  }
  
  ## convert vectors to matrices
  dim(mEm) = c(q[length(q)] - 1, length(l))
  dim(mEn) = c(q[length(q)] - 1, length(l))
  dim(mRm) = c(q[length(q)] - 1, length(l))
  dim(mRn) = c(q[length(q)] - 1, length(l))
  dim(mRt) = c(q[length(q)] - 1, length(l))
  dim(ol) =  c(q[length(q)] - 1, length(l))
  
  ## assign row- and col-names
  rownames(mEm) <- q
  colnames(mEm) <- l
  
  rownames(mEn) <- q
  colnames(mEn) <- l

  rownames(mRm) <- q
  colnames(mRm) <- l

  rownames(mRn) <- q
  colnames(mRn) <- l

  rownames(mRt) <- q
  colnames(mRt) <- l

  rownames(ol) <- q
  colnames(ol) <- l

  ## define output vector
  q.max <- rep(NA, length(l))
  
  ## loop through all l-values
  for (i in 1:length(l)) {
    ## assign maximum number of q, i.e. last q before explained variance drops
    mRti <- mRt[1:(nrow(mRt) - 1),i]
    mRtj <- mRt[2:nrow(mRt),i]
    dmRt <- mRtj - mRti
    q.max[i] <- q[dmRt < 0][1]
  }
  
  ## plot outputs
  ## adjust plot margins
  par(oma = c(0, 1, 0, 0))
  ## read plot-independent plot parameters and check/set default values
  extraArgs <- list(...)
  colour <- if("colour" %in% names(extraArgs)) {extraArgs$colour} else
  {seq(1, length(l))}
  legend.cex <- if("cex" %in% names(extraArgs)) {extraArgs$cex} else
  {1}
  legend.lty <- if("lty" %in% names(extraArgs)) {extraArgs$lty} else
  {1}
  legend.text <- if("legend" %in% names(extraArgs)) {extraArgs$legend} else
  {round(l, 3)}
  legend.title <- if("title" %in% names(extraArgs)) {extraArgs$title} else
  {"l"}
  
  ## optional plot mEm
  if(plot == "mEm") {
    ## read additional arguments list and check/set default values
    extraArgs <- list(...)
    main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
    {c(expression(paste("Model error (", bar(E)[m], ")", sep = "")), 
       expression(paste("Model error (", bar(E)[m], ")", sep = "")))}
    xlab <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab} else
    {c("Number of end-members (q)", 
       "Number of end-members (q)")}
    ylab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else
    {c(expression(paste("Weight limit (", l[w], ")", sep = "")),
       expression(paste("Error (", bar(E)[m], ")", sep = "")))}
    
    ## setup plot area
    par(mfcol = c(1, 2))
    
    ## plot matrix q-l as image
    image(q, l, mEm,
          col = rainbow(n = 250, start = 0.2, end = 1),
          main = main[1],
          xlab = xlab[1], 
          ylab = ylab[1],
          zlim = c(min(mEm, na.rm = TRUE), max(mEm, na.rm = TRUE)))
    if(length(l) > 1) contour(q, l, mEm, add = TRUE,
                               zlim = c(min(mEm, na.rm = TRUE), max(mEm, na.rm = TRUE)))
    
    ## plot mEm vs. q
    plot(seq(min(q), max(q)), mEm[,1], type = "l", 
         main = main[2],
         xlab = xlab[2],
         ylab = ylab[2],
         ylim = range(mEm, na.rm = TRUE),
         col = colour[1])
    if(nrow(mEm > 1)) {
     
      for(i in 2:length(l)) {
        
        lines(mEm[,i], col = colour[i])
      } 
    }
    
    ## optionally add legend
    if(missing(legend) == FALSE) {
      legend.position <- legend
      legend(x = legend.position,
             legend = legend.text,
             col = colour,
             cex = legend.cex,
             lty = legend.lty,
             title = legend.title)
    }
    
    ## reset plot area
    par(mfcol = c(1,1))
  }
  
  ## optional plot mEn
  if(plot == "mEn") {
    ## read additional arguments list and check/set default values
    extraArgs <- list(...)
    main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
    {c(expression(paste("Model error (", bar(E)[n], ")", sep = "")), 
       expression(paste("Model error (", bar(E)[n], ")", sep = "")))}
    xlab <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab} else
    {c("Number of end-members (q)", 
       "Number of end-members (q)")}
    ylab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else
    {c(expression(paste("Weight limit (", l[w], ")", sep = "")),
       expression(paste("Error (", bar(E)[n], ")", sep = "")))}
    
    ## setup plot area
    par(mfcol = c(1,2))
    
    ## plot matrix q-l as image
    image(q, l, mEn, 
          col = rainbow(n = 250, s = 1, v = 1, start = 0.2, end = 1),
          main = main[1],
          xlab = xlab[1], 
          ylab = ylab[1],
          zlim = c(min(mEn, na.rm = TRUE), max(mEn, na.rm = TRUE)))
    if(length(l) > 1) contour(q, l, mEn, add = TRUE,
                               zlim = c(min(mEn, na.rm = TRUE), max(mEn, na.rm = TRUE)))
    
    ## plot mEn vs. q
    plot(seq(min(q), max(q)), mEn[,1], type = "l", 
         main = main[2],
         xlab = xlab[2],
         ylab = ylab[2],
         ylim = range(mEn, na.rm = TRUE),
         col = colour[1])
    if(nrow(mEm > 1)) for (i in 2:length(l)) {lines(mEn[,i], col = colour[i])}
    
    ## optionally add legend
    if(missing(legend) == FALSE) {
      legend.position <- legend
      legend(x = legend.position,
             legend = legend.text,
             col = colour,
             cex = legend.cex,
             lty = legend.lty,
             title = legend.title)
    }
    
    ## reset plot area
    par(mfcol = c(1, 1))
  }
  
  ## optional plot mRm
  if(plot == "mRm") {
    ## read additional arguments list and check/set default values
    extraArgs <- list(...)
    main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
    {c(expression(paste("Explained variance (", bar(R^2)[m], ")", sep = "")), 
       expression(paste("Explained variance (", bar(R^2)[m], ")", sep = "")))}
    xlab <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab} else
    {c("Number of end-members (q)", 
       "Number of end-members (q)")}
    ylab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else
    {c(expression(paste("Weight limit (", l[w], ")", sep = "")),
       expression(paste("Explained variance (", bar(R^2)[m], ")", sep = "")))}
    
    ## setup plot area
    par(mfcol = c(1,2))
    
    ## plot matrix q-l as image
    image(q, l, mRm, 
          col = rainbow(n = 250, s = 1, v = 1, start = 0.2, end = 1),
          main = main[1],
          xlab = xlab[1], 
          ylab = ylab[1],
          zlim = c(min(mRm, na.rm = TRUE), max(mRm, na.rm = TRUE)))
    if(length(l) > 1) contour(q, l, mRm, add = TRUE,
                               zlim = c(min(mRm, na.rm = TRUE), max(mRm, na.rm = TRUE)))
    
    ## plot mRm vs. q
    plot(seq(min(q), max(q)), mRm[,1], type = "l", 
         main = main[2],
         xlab = xlab[2],
         ylab = ylab[2],
         ylim = range(mRm, na.rm = TRUE),
         col = colour[1])
    if(nrow(mRm > 1)) for (i in 2:length(l)) {lines(mRm[,i], col = colour[i])}
    
    ## optionally add legend
    if(missing(legend) == FALSE) {
      legend.position <- legend
      legend(x = legend.position,
             legend = legend.text,
             col = colour,
             cex = legend.cex,
             lty = legend.lty,
             title = legend.title)
    }
    
    ## reset plot area
    par(mfcol = c(1,1))
  }
  
  ## optional plot mRn
  if(plot == "mRn") {
    ## read additional arguments list and check/set default values
    extraArgs <- list(...)
    main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
    {c(expression(paste("Explained variance (", bar(R^2)[n], ")", sep = "")), 
       expression(paste("Explained variance (", bar(R^2)[n], ")", sep = "")))}
    xlab <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab} else
    {c("Number of end-members (q)", 
       "Number of end-members (q)")}
    ylab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else
    {c(expression(paste("Weight limit (", l[w], ")", sep = "")),
       expression(paste("Explained variance (", bar(R^2)[n], ")", sep = "")))}
    
    ## setup plot area
    par(mfcol = c(1,2))
    
    ## plot matrix q-l as image
    image(q, l, mRn, 
          col = rainbow(n = 250, s = 1, v = 1, start = 0.2, end = 1),
          main = main[1],
          xlab = xlab[1], 
          ylab = ylab[1],
          zlim = c(min(mRn, na.rm = TRUE), max(mRn, na.rm = TRUE)))
    if(length(l) > 1) contour(q, l, mRn, add = TRUE,
                               zlim = c(min(mRn, na.rm = TRUE), max(mRn, na.rm = TRUE)))
    
    ## plot mRn vs. q
    plot(seq(min(q), max(q)), mRn[,1], type = "l", 
         main = main[2],
         xlab = xlab[2],
         ylab = ylab[2],
         ylim = range(mRn, na.rm = TRUE),
         col = colour[1])
    if(nrow(mRn > 1)) for (i in 2:length(l)) {lines(mRn[,i], col = colour[i])}
    
    ## optionally add legend
    if(missing(legend) == FALSE) {
      legend.position <- legend
      legend(x = legend.position,
             legend = legend.text,
             col = colour,
             cex = legend.cex,
             lty = legend.lty,
             title = legend.title)
    }
    
    ## reset plot area
    par(mfcol = c(1,1))
  }
  
  ## optional plot mRt
  if(plot == "mRt") {
    ## read additional arguments list and check/set default values
    extraArgs <- list(...)
    main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
    {c(expression(paste("Explained variance (", bar(R^2)[t], ")", sep = "")), 
       expression(paste("Explained variance (", bar(R^2)[t], ")", sep = "")))}
    xlab <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab} else
    {c("Number of end-members (q)", 
       "Number of end-members (q)")}
    ylab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else
    {c(expression(paste("Weight limit (", l[w], ")", sep = "")),
       expression(paste("Explained variance (", bar(R^2)[t], ")", sep = "")))}
    
    ## setup plot area
    par(mfcol = c(1,2))
    
    ## plot matrix q-l as image
    image(q, l, mRt, 
          col = rainbow(n = 250, s = 1, v = 1, start = 0.2, end = 1),
          main = main[1],
          xlab = xlab[1], 
          ylab = ylab[1],
          zlim = c(min(mRt, na.rm = TRUE), max(mRt, na.rm = TRUE)))
    if(length(l) > 1) contour(q, l, mRt, add = TRUE,
                               zlim = c(min(mRt, na.rm = TRUE), max(mRt, na.rm = TRUE)))
    
    ## plot mRt vs. q
    plot(seq(min(q), max(q)), mRt[,1], type = "l", 
         main = main[2],
         xlab = xlab[2],
         ylab = ylab[2],
         ylim = range(mRt, na.rm = TRUE),
         col = colour[1])
    if(nrow(mRt > 1)) for (i in 2:length(l)) {lines(mRt[,i], col = colour[i])}
    
    ## optionally add legend
    if(missing(legend) == FALSE) {
      legend.position <- legend
      legend(x = legend.position,
             legend = legend.text,
             col = colour,
             cex = legend.cex,
             lty = legend.lty,
             title = legend.title)
    }
    
    ## reset plot area
    par(mfcol = c(1,1))
  }  
  
  ## optional plot ol
  if(plot == "ol") {
    ## read additional arguments list and check/set default values
    extraArgs <- list(...)
    main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
    {"Overlapping end-member loadings"}
    xlab <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab} else
    {"Number of end-members (q)"}
    ylab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else
    {expression(paste("Weight limit (", l[w], ")", sep = ""))}
    cex <- if("cex" %in% names(extraArgs)) {extraArgs$cex} else
    {0.6}
    
    ## plot matrix q-l as image
    image(q, l, ol, 
          col = rev(rainbow(n = 250, s = 1, v = 1, start = 0.2, end = 1)),
          main = main, 
          xlab = xlab, 
          ylab = ylab, 
          zlim = c(min(ol, na.rm = TRUE), max(ol, na.rm = TRUE)))
    text(q.t, l.t, ol, cex = cex)
  }
  
  ## optionally add pm
  if(pm == TRUE) {pm <- check.data(matrix(runif(4), ncol = 2),
                                   5, 0.01, 100, invisible = FALSE)}
  
  ## readjust plot margins
  par(oma = c(0, 0, 0, 0))
  
  ## return result
  list(mEm   = mEm,
       mEn   = mEn,
       mRm   = mRm,
       mRn   = mRn,
       mRt   = mRt,
       ol    = ol, 
       q.max = q.max)
}
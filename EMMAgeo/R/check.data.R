#' Function to check data consistency.
#' 
#' The input data matrix (X), number of end-members (q), weight transformation
#' limits (l) and constant sum scaling parameter (c) are checked for
#' consistency. This includes checking for absence of missing values, columns
#' containing only zero-values and for numeric data type of variables.
#' Furthermore, a check tests if l is below the maximum possible value,
#' preventing numerical instability prior to factor rotation.
#' 
#' 
#' @param X Numeric matrix with m samples (rows) and n variables (columns).
#' @param q Numeric scalar with number of end-members to be modelled.
#' @param l Numeric scalar or vector specifying the weight transformation
#' limit, i.e.  quantile.
#' @param c Numeric scalar specifying the constant sum scaling parameter, e.g.
#' 1, 100, 1000.
#' @param invisible Logical scalar setting visibility option.
#' @param \dots Further arguments passed to the function.
#' @return Character vector with test results.
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
#' ## perform data set check
#' check.data(X = X, 
#'            q = 6, 
#'            l = seq(from = 0, 
#'                    to = 0.2, 
#'                    by = 0.01), 
#'            c = 1)
#' 
#' @export check.data
check.data <- function(
  X, 
  q, 
  l,
  c,
  invisible = TRUE,
  ...
){
  
  ## create result vector
  result <- NA
  
  ## check for l vs. lw
  if("lw" %in% names(list(...))) {
    stop('Parameter "lw" is depreciated. Use "l" instead.')
  }
  
  ## test data type
  result[length(result) + 1] <- ifelse(is.numeric(X) == FALSE, 
    "    Warning: data matrix X contains non-numeric values.",
    "Data matrix passed test... OK")
  
  result[length(result) + 1] <- ifelse(is.numeric(q) == FALSE, 
    "    Warning: end-member vector q has non-numeric values.",
    "End-member vector passed test... OK")
 
  result[length(result) + 1] <- ifelse(is.numeric(l) == FALSE, 
    paste("    Warning: weight transformation limit vector l", 
          "contains non-numeric values."),
    "Weight transformation limit vector passed test... OK")
  
  result[length(result) + 1] <- ifelse(is.numeric(c) == FALSE,
    "    Warning: constant sum scaling parameter c is non-numeric.",
    "Scaling parameter passed test... OK")
  
  ## check for samples with NA-values
  result[length(result) + 1] <- ifelse(sum(complete.cases(X)) < nrow(X),
    paste("    Warning: The following samples comprise NA-values: ",
          seq(1, nrow(X))[!complete.cases(X)],
          ".",
          sep = ""), "NA-test passed... OK")
  if(sum(complete.cases(X)) < nrow(X)) {X <- X[complete.cases(X),]}
  
  ## test if columns contain only 0 values
  X.0 <- apply(X, 2, sum, na.rm = TRUE)
  X.unmet <- seq(1, ncol(X))[X.0 == 0]
  m.unmet <- paste(X.unmet, collapse = ", ")
  if(length(X.unmet) > 0) {
    paste("    Warning: the following columns contain only zero values: ",
      m.unmet, ".",
      sep = "")} else {
        result[length(result) + 1] <- "Test for zero-only values passed... OK"
      }
   
  ## test range of vector l
  l.max <- test.l(X, l)$l.max
  if(length(l.max) == 0) {
    result[length(result) + 1] <- paste("    Warning: weight transformation",
          "limit is out of range.")
  } else {
  result[length(result) + 1] <- ifelse(max(l) > l.max,
    paste("    Note: weight transformation limit(s) are out",
          "of range. Maximum value is", l.max),
    "Maximum weight transformation limit value passed test... OK")
  }
  
  ## test if all samples sum up to constant value
  X.c <- round(apply(X, 1, sum) - rep(c, nrow(X)), 5)
  X.unmet <- seq(1, nrow(X))[X.c != 0]
  m.unmet <- paste(X.unmet, collapse = ", ")
  result[length(result) + 1] <- ifelse(length(X.unmet) >= 1,
      paste("    Note: the following rows do not sum up to the specified c: ",
            m.unmet, ".", sep = ""),
      "All samples sum up to constant sum... OK")
  
  if(invisible == FALSE) {
    
    par(new = TRUE)
    plot(NA, xlim = c(0, 10), ylim = c(0, 10), 
         main = "", xlab = "", ylab = "", 
         axes = FALSE, frame.plot = FALSE)
    
    n_frames      <- 25
    t_animation   <- 2.5
    
    dt            <- t_animation / n_frames
    x1            <- seq(0, 10, length.out = n_frames)
    y1            <- rep(1.5, n_frames)
    r1            <- 0.5
    
    x2            <- seq(0, 16, length.out = n_frames)
    y2            <- rep(8.5, n_frames)
    r2            <- 0.5
    
    x4            <- seq(11, 0, length.out = n_frames)
    y4            <- rep(5, n_frames)
    r4            <- 0.5
    
    # set angles for each step of mouth opening
    angles_mouth <- rep(c(0.01, 0.25, 0.5, 0.25), length.out = n_frames)
    
    for(i in 1:n_frames){
      # define circles
      shape::filledcircle(r1 = r1, r2 = 0.00001, mid = c(x1[i], y1[i]), 
                   from = angles_mouth[i], to = 2 * pi - angles_mouth[i], 
                   col = "yellow")
      shape::filledcircle(r1 = r2, r2 = 0.00001, mid = c(x2[i], y2[i]), 
                   from = angles_mouth[i], to = 2 * pi - angles_mouth[i], 
                   col = "yellow")
      shape::filledcircle(r1 = r4, r2 = 0.00001, mid = c(x4[i], y4[i]), 
                   from = angles_mouth[i] + 3, to = 2 * pi - angles_mouth[i] + 3, 
                   col = "yellow")
      
      # define eyes
      points(x1[i] + 0.2, y1[i] + 0.75, pch = 21, bg = 1, cex = 0.7)
      points(x2[i] + 0.2, y2[i] + 0.75, pch = 21, bg = 1, cex = 0.7)
      points(x4[i] - 0.05, y4[i] + 0.75, pch = 21, bg = 1, cex = 0.7)
      
      Sys.sleep(dt)
      
      # define white background
      shape::plotcircle(r = 1.1 * r1, mid = c(x1[i], y1[i]), 
                 col = "white", lcol = "white")
      shape::plotcircle(r = 1.1 * r2, mid = c(x2[i], y2[i]), 
                 col = "white", lcol = "white")
      shape::plotcircle(r = 1.1 * r4, mid = c(x4[i], y4[i]), 
                 col = "white", lcol = "white")
    }
  }
  
  ## return result
  return(result[2:length(result)])
}
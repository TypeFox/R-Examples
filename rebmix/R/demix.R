demix <- function(x = NULL, 
  Preprocessing = NULL, 
  pdf = NULL,
  k = NULL, 
  xmin = NULL, 
  xmax = NULL, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  if (is.null(x)) {
    stop(sQuote("x"), " must not be NULL!", call. = FALSE)
  }

  if ((!is.numeric(x)) && (!is.data.frame(x))) {
    stop(sQuote("x"), " numeric or data frame is requested!", call. = FALSE)
  }
  
  x <- as.matrix(x)
  
  d <- ncol(x)
  n <- nrow(x) 
  
  if (is.null(Preprocessing)) {
    stop(sQuote("Preprocessing"), " must not be NULL!", call. = FALSE)
  }

  if (!is.character(Preprocessing)) {
    stop(sQuote("Preprocessing"), " character is requested!", call. = FALSE)
  }

  Preprocessing <- match.arg(Preprocessing, .rebmix$Preprocessing, several.ok = FALSE)    
  
  if (is.null(pdf)) {
    stop(sQuote("pdf"), " must not be NULL!", call. = FALSE)
  }

  if (!is.character(pdf)) {
    stop(sQuote("pdf"), " character vector is requested!", call. = FALSE)
  }

  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)
  
  Variables <- NULL
    
  for (i in 1:length(.rebmix$pdf)) {
    Variables[which(pdf == .rebmix$pdf[i])] <- .rebmix$pdf.Variables[i]
  }  
  
  if (is.null(k)) {
    stop(sQuote("k"), " must not be NULL!", call. = FALSE)
  }

  if (!is.wholenumber(k)) {
    stop(sQuote("k"), " integer is requested!", call. = FALSE)
  }

  if (!(k > 0)) {
    stop(sQuote("k"), " must be greater than 0!", call. = FALSE)
  }
  
  if (is.null(xmin)) {
    xmin <- apply(x, 2, min)
  }
  else {
    xmin <- xmin
  }
  
  if (is.null(xmax)) {
    xmax <- apply(x, 2, max)
  }
  else {
    xmax <- xmax
  }

  if (Preprocessing == .rebmix$Preprocessing[1]) {
    h <- array(data = 0.0, dim = d, dimnames = NULL)  
    y0 <- array(data = 0.0, dim = d, dimnames = NULL)  
    
    for (i in 1:d) {
      if (Variables[i] == .rebmix$Variables[1]) {
        h[i] = (xmax[i] - xmin[i]) / k; y0[i] = xmin[i] + 0.5 * h[i]
      }
      else 
      if (Variables[i] == .rebmix$Variables[2]) {
        h[i] = 1.0; y0[i] = xmin[i]
      }
    }    

    output <- .C("RPreprocessingH",
      h = as.double(h),
      y0 = as.double(y0),
      length.pdf = as.integer(d),
      pdf = as.character(pdf),
      k = as.integer(k),
      n = as.integer(n),
      d = as.integer(d),
      x = as.double(unlist(x)),
      y = double(n * (d + 1)),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in preprocessing!", call. = FALSE); return(NA)
    }

    length(output$y) <- output$k * (output$d + 1); dim(output$y) <- c(output$k, output$d + 1)

    output$y[, d + 1] <- output$y[, d + 1] / prod(h) / n

    output <- as.data.frame(output$y, stringsAsFactors = FALSE)

    colnames(output) <- c(paste("x", if (d > 1) 1:d else "", sep = ""), "f")
  } 
  else 
  if (Preprocessing == .rebmix$Preprocessing[2]) {
    h <- array(data = 0.0, dim = d, dimnames = NULL)  
    
    for (i in 1:d) {
      if (Variables[i] == .rebmix$Variables[1]) {
        h[i] = (xmax[i] - xmin[i]) / k
      }
      else 
      if (Variables[i] == .rebmix$Variables[2]) {
        h[i] = 1.0
      }
    }   
      
    output <- .C("RPreprocessingPW",
      h = as.double(h),
      n = as.integer(n),
      d = as.integer(d),
      x = as.double(unlist(x)),
      y = double(n * (d + 2)),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in preprocessing!", call. = FALSE); return(NA)
    }

    dim(output$y) <- c(n, d + 2)

    output$y[, d + 2] <- output$y[, d + 2] / prod(h) / n

    output <- as.data.frame(output$y[, -(d + 1)], stringsAsFactors = FALSE)

    colnames(output) <- c(paste("x", if (d > 1) 1:d else "", sep = ""), "f")
  } 
  else
  if (Preprocessing == .rebmix$Preprocessing[3]) {
    h <- array(data = 0.0, dim = d, dimnames = NULL)  
    
    for (i in 1:d) {
      h[i] = xmax[i] - xmin[i]
    } 

    output <- .C("RPreprocessingKNN",
      k = as.integer(k),
      h = as.double(h),
      n = as.integer(n),
      d = as.integer(d),
      x = as.double(unlist(x)),
      y = double(n * (d + 3)),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in preprocessing!", call. = FALSE); return(NA)
    }

    dim(output$y) <- c(n, d + 3)

    output$y[, d + 2] <- k / output$y[, d + 2] / n

    output <- as.data.frame(output$y[, c(-(d + 1), -(d + 3))], stringsAsFactors = FALSE)

    colnames(output) <- c(paste("x", if (d > 1) 1:d else "", sep = ""), "f")
  }
  
  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  invisible(output)
} ## demix

# Class RNGMIX

setClass(Class = "RNGMIX",
slots = c(Dataset.name = "character",
  rseed = "numeric",
  n = "numeric",
  Theta = "list",
  Dataset = "list",
  Zt = "factor",
  w = "numeric",
  Variables = "character",
  ymin = "numeric",
  ymax = "numeric"),
prototype = list(rseed = -1))

setMethod("initialize", "RNGMIX", 
function(.Object, ...,
  Dataset.name,
  rseed,
  n,
  Theta)
{
  # Dataset.name.
  
  if (missing(Dataset.name) || (length(Dataset.name) == 0)) {
    stop(sQuote("Dataset.name"), " must not be empty!", call. = FALSE)
  }

  if (!is.character(Dataset.name)) {
    stop(sQuote("Dataset.name"), " character vector is requested!", call. = FALSE)
  }

  # rseed.

  if (missing(rseed) || (length(rseed) == 0)) rseed <- .Object@rseed
  
  if (!is.wholenumber(rseed)) {
    stop(sQuote("rseed"), " integer is requested!", call. = FALSE)
  }

  length(rseed) <- 1
  
  if (rseed > -1) {
    stop(sQuote("rseed"), " must be less than 0!", call. = FALSE)
  }

  # n.
  
  if (missing(n) || (length(n) == 0)) {
    stop(sQuote("n"), " must not be empty!", call. = FALSE)
  }

  if (!is.wholenumber(n)) {
    stop(sQuote("n"), " integer vector is requested!", call. = FALSE)
  }
  
  if (!all(n > 0)) {
    stop("all ", sQuote("n"), " must be greater than 0!", call. = FALSE)
  }
  
  c <- length(n)
  
  # Theta.  
  
  if (missing(Theta) || (length(Theta) == 0)) {
    stop(sQuote("Theta"), " must not be empty!", call. = FALSE)
  }
  
  if (!is.list(Theta)) {
    stop(sQuote("Theta"), " list is requested!!", call. = FALSE)
  }  
  
  Names <- names(Theta)
  
  if (length(grep("pdf", Names)) == 0) {
    stop(sQuote("pdf1"), " in " , sQuote("Theta"), " character vector is requested!", call. = FALSE)
  }
  
  if (length(grep("theta1", Names)) == 0) {
    stop(sQuote("theta1.1"), " in " , sQuote("Theta"), " numeric vector is requested!", call. = FALSE)
  }
  
  if (length(grep("theta2", Names)) == 0) {
    stop(sQuote("theta2.1"), " in " , sQuote("Theta"), " numeric vector is requested!", call. = FALSE)
  }     
  
  length(grep("pdf", Names))  
  
  j <- 0; length.pdf <- length(Theta[[1]])
  
  for (i in grep("pdf", Names)) {  
    pdf <- as.character(Theta[[i]])
  
    pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)
    
    if (length(pdf) != length.pdf) {
      stop("lengths of ", sQuote("pdfi"), " in " , sQuote("Theta"), " must be equal!", call. = FALSE)
    }    
  
    Theta[[i]] <- pdf; j <- j + 1
  }

  if ((length.pdf > 1) && (j != c)) {
    stop(sQuote("pdfi"), " in " , sQuote("Theta"), " and ", sQuote("n"), " must match!", call. = FALSE)
  } 
  
  j <- 0; length.theta1 <- length(Theta[[2]])
  
  for (i in grep("theta1", Names)) {  
    theta1 <- as.numeric(Theta[[i]])
   
    if (length(theta1) != length.theta1) {
      stop("lengths of ", sQuote("theta1.i"), " in " , sQuote("Theta"), " must be equal!", call. = FALSE)
    }    

    j <- j + 1
  }

  if ((length.pdf > 1) && (j != c)) {
    stop(sQuote("theta1.i"), " in " , sQuote("Theta"), " and ", sQuote("n"), " must match!", call. = FALSE)
  } 
  
  j <- 0; length.theta2 <- length(Theta[[3]])
  
  for (i in grep("theta2", Names)) {  
    theta2 <- as.numeric(Theta[[i]])
   
    if (length(theta2) != length.theta2) {
      stop("lengths of ", sQuote("theta2.i"), " in " , sQuote("Theta"), " must be equal!", call. = FALSE)
    }

    j <- j + 1    
  }

  if ((length.pdf > 1) && (j != c)) {
    stop(sQuote("theta2.i"), " in " , sQuote("Theta"), " and ", sQuote("n"), " must match!", call. = FALSE)
  } 
 
  # Variables.

  for (i in 1:length(.rebmix$pdf)) {
    .Object@Variables[which(pdf == .rebmix$pdf[i])] <- .rebmix$pdf.Variables[i]
  }
  
  callNextMethod(.Object, ...,
    Dataset.name = Dataset.name,
    rseed = rseed,
    n = n,
    Theta = Theta)
}) ## initialize

setMethod("show",
          signature(object = "RNGMIX"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RNGMIX is requested!", call. = FALSE)
  }
  
  cat("An object of class ", "\"", class(object), "\"", "\n", sep = "")
  
  cat("Slot \"w\":", "\n", sep = "")

  print(object@w, quote = FALSE)

  cat("Slot \"ymin\":", "\n", sep = "")

  print(object@ymin, quote = FALSE)

  cat("Slot \"ymax\":", "\n", sep = "")

  print(object@ymax, quote = FALSE)

  rm(list = ls())
}) ## show

# Class RNGMVNORM
                          
setClass("RNGMVNORM", contains = "RNGMIX")

setMethod("show",
          signature(object = "RNGMVNORM"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RNGMVNORM is requested!", call. = FALSE)
  }
  
  cat("An object of class ", "\"", class(object), "\"", "\n", sep = "")
  
  cat("Slot \"w\":", "\n", sep = "")

  print(object@w, quote = FALSE)

  cat("Slot \"ymin\":", "\n", sep = "")

  print(object@ymin, quote = FALSE)

  cat("Slot \"ymax\":", "\n", sep = "")

  print(object@ymax, quote = FALSE)

  rm(list = ls())
}) ## show

# Class REBMIX

setClass(Class = "REBMIX",
slots = c(Dataset = "list",
  Preprocessing = "character",
  cmax = "numeric",
  Criterion = "character",
  Variables = "character",
  pdf = "character",
  theta1 = "numeric",
  theta2 = "numeric",
  K = "ANY",
  y0 = "numeric",
  ymin = "numeric",
  ymax = "numeric",                                      
  ar = "numeric",   
  Restraints = "character",
  w = "list",
  Theta = "list",
  summary = "ANY",
  pos = "numeric",
  opt.c = "list",
  opt.IC = "list",
  opt.logL = "list",
  opt.D = "list",
  all.K = "list",
  all.IC = "list"),
prototype = list(cmax = 15,
  Criterion = "AIC",
  ar = 0.1,
  Restraints = "loose",
  pos = 1)) 

setMethod("initialize", "REBMIX", 
function(.Object, ...,
  Dataset,
  Preprocessing,
  cmax,
  Criterion,
  pdf,
  theta1,
  theta2,
  K,
  y0,
  ymin,
  ymax,
  ar,
  Restraints)
{
  # Dataset.
  
  if (missing(Dataset) || (length(Dataset) == 0)) {
    stop(sQuote("Dataset"), " must not be empty!", call. = FALSE)
  }
  
  if (!is.list(Dataset)) {
    stop(sQuote("Dataset"), " list is requested!", call. = FALSE)
  }   
  
  if (is.null(names(Dataset))) {
    names(Dataset) <- paste("dataset", 1:length(Dataset), sep = "")
  }

  if (!all(unlist(lapply(Dataset, is.data.frame)))) {
    stop(sQuote("Dataset"), " list of data frames is requested!", call. = FALSE)
  }

  d <- unique(unlist(lapply(Dataset, ncol)))

  if (length(d) != 1) {
    stop(sQuote("Dataset"), " numbers of columns in data frames must be equal!", call. = FALSE)
  }

  if (!all(unlist(lapply(Dataset, ncol)) > 0)) {
    stop(sQuote("Dataset"), " numbers of columns in data frames must be greater than 0!", call. = FALSE)
  }

  if (!all(unlist(lapply(Dataset, nrow)) > 1)) {
    stop(sQuote("Dataset"), " numbers of rows in data frames must be greater than 1!", call. = FALSE)
  }

  # Preprocessing.
  
  if (missing(Preprocessing) || (length(Preprocessing) == 0)) {
    stop(sQuote("Preprocessing"), " must not be empty!", call. = FALSE)
  }

  if (!is.character(Preprocessing)) {
    stop(sQuote("Preprocessing"), " character vector is requested!", call. = FALSE)
  }   

  Preprocessing <- match.arg(Preprocessing, .rebmix$Preprocessing, several.ok = TRUE)
  
  # cmax.
  
  if (missing(cmax) || (length(cmax) == 0)) cmax <- .Object@cmax
  
  if (!is.wholenumber(cmax)) {
    stop(sQuote("cmax"), " integer is requested!", call. = FALSE)
  }

  if (cmax < 1) {
    stop(sQuote("cmax"), " must be greater than 0!", call. = FALSE)
  }  

  # Criterion.

  if (missing(Criterion) || (length(Criterion) == 0)) Criterion <- .Object@Criterion

  if (!is.character(Criterion)) {
    stop(sQuote("Criterion"), " character vector is requested!", call. = FALSE)
  } 

  Criterion <- match.arg(Criterion, .rebmix$Criterion, several.ok = TRUE)
  
  # pdf.
  
  if (missing(pdf) || (length(pdf) == 0)) {
    stop(sQuote("pdf"), " must not be empty!", call. = FALSE)
  }

  if (!is.character(pdf)) {
    stop(sQuote("pdf"), " character vector is requested!", call. = FALSE)
  } 

  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)

  if (length(pdf) != d) {
    stop("lengths of ", sQuote("pdf"), " and ", sQuote("d"), " must match!", call. = FALSE)
  }
  
  # theta1.
  
  if (missing(theta1) || (length(theta1) == 0)) {
    theta1 <- .Object@theta1
  }
  else {
    if (length(theta1) != d) {
      stop("lengths of ", sQuote("theta1"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
    
    class(theta1) <- "numeric"
  }
  
  # theta2.
  
  if (missing(theta2) || (length(theta2) == 0)) {
    theta2 <- .Object@theta2
  }
  else {
    if (length(theta2) != d) {
      stop("lengths of ", sQuote("theta2"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
    
    class(theta2) <- "numeric"
  }   
  
  # K. 

  if (missing(K) || (length(K) == 0)) {
    stop(sQuote("K"), " must not be empty!", call. = FALSE)
  }

  if (is.list(K)) {
    if (!all(unlist(lapply(K, is.wholenumber) == TRUE))) {
      stop(sQuote("K"), " list of integer vectors is requested!", call. = FALSE)
    }

    if (!all(unlist(lapply(K, function(x) all(x > 0)))) == TRUE) {
      stop("all ", sQuote("K"), " must be greater than 0!", call. = FALSE)
    }

    if (length(K) != length(Preprocessing)) {
      stop("lengths of ", sQuote("Preprocessing"), " and ", sQuote("K"), " must match!", call. = FALSE)
    }
  }
  else 
  if (is.numeric(K)) {
    if (!is.wholenumber(K)) {
      stop(sQuote("K"), " integer vector is requested!", call. = FALSE)
    }

    if (!all(K > 0)) {
      stop("all ", sQuote("K"), " must be greater than 0!", call. = FALSE)
    }
  }
  else {
    stop(sQuote("K"), " list of integer vectors or integer vector is requested!", call. = FALSE)
  }
  
  # y0.
  
  if (missing(y0) || (length(y0) == 0)) {
    y0 <- .Object@y0  
  }
  else {
    if (length(y0) != d) {
      stop("lengths of ", sQuote("y0"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }      
  }
  
  # ymin.
  
  if (missing(ymin) || (length(ymin) == 0)) {
    ymin <- .Object@ymin  
  }
  else {
    if (length(ymin) != d) {
      stop("lengths of ", sQuote("ymin"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }      
  }
  
  # ymax.
  
  if (missing(ymax) || (length(ymax) == 0)) {
    ymax <- .Object@ymax  
  }
  else {
    if (length(ymax) != d) {
      stop("lengths of ", sQuote("ymax"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }      
  }
  
  # ar.
  
  if (missing(ar) || (length(ar) == 0)) ar <- .Object@ar
  
  if (!is.numeric(ar)) {
    stop(sQuote("ar"), " numeric is requested!", call. = FALSE)
  }
  
  length(ar) <- 1

  if ((ar <= 0.0) || (ar > 1.0)) {
    stop(sQuote("ar"), " must be greater than 0.0 and less or equal than 1.0!", call. = FALSE)
  }     

  # Restraints.

  if (missing(Restraints) || (length(Restraints) == 0)) Restraints <- .Object@Restraints
  
  if (!is.character(Restraints)) {
    stop(sQuote("Restraints"), " character is requested!", call. = FALSE)
  }
  
  Restraints <- match.arg(Restraints, .rebmix$Restraints, several.ok = FALSE)
  
  # Variables.

  for (i in 1:length(.rebmix$pdf)) {
    .Object@Variables[which(pdf == .rebmix$pdf[i])] <- .rebmix$pdf.Variables[i]
  }
  
  callNextMethod(.Object, ...,
    Dataset = Dataset,
    Preprocessing = Preprocessing,
    cmax = cmax,
    Criterion = Criterion,
    pdf = pdf,
    theta1 = theta1,
    theta2 = theta2,
    K = K,
    y0 = y0,
    ymin = ymin,
    ymax = ymax,
    ar = ar,
    Restraints = Restraints)
}) ## initialize
                    
setMethod("show",
          signature(object = "REBMIX"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  if (!is.wholenumber(object@pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(object@pos) <- 1

  if ((object@pos < 1) || (object@pos > nrow(object@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(object@summary), "!", call. = FALSE)
  }  
  
  cat("An object of class ", "\"", class(object), "\"", "\n", sep = "")
  
  cat("Slot \"w\":", "\n", sep = "")

  print(object@w[[object@pos]], quote = FALSE)

  cat("Slot \"Theta\":", "\n", sep = "")

  print(object@Theta[[object@pos]], quote = FALSE)

  cat("Slot \"summary\":", "\n", sep = "")
  
  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(object@summary), nomatch = 0)
  
  DF <- object@summary[object@pos, p]

  is.num <- sapply(DF, is.number); DF[is.num] <- lapply(DF[is.num], as.number)

  print(DF, quote = FALSE) 

  rm(list = ls())
}) ## show

## Class REBMVNORM

setClass("REBMVNORM", contains = "REBMIX")

setMethod("show",
          signature(object = "REBMVNORM"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMVNORM is requested!", call. = FALSE)
  }
  
  if (!is.wholenumber(object@pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(object@pos) <- 1

  if ((object@pos < 1) || (object@pos > nrow(object@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(object@summary), "!", call. = FALSE)
  }  
  
  cat("An object of class ", "\"", class(object), "\"", "\n", sep = "")
  
  cat("Slot \"w\":", "\n", sep = "")

  print(object@w[[object@pos]], quote = FALSE)

  cat("Slot \"Theta\":", "\n", sep = "")

  print(object@Theta[[object@pos]], quote = FALSE)

  cat("Slot \"summary\":", "\n", sep = "")
  
  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(object@summary), nomatch = 0)
  
  DF <- object@summary[object@pos, p]

  is.num <- sapply(DF, is.number); DF[is.num] <- lapply(DF[is.num], as.number)

  print(DF, quote = FALSE)

  rm(list = ls())
}) ## show

## Class REBMIX.boot

setClass("REBMIX.boot",
slots = c(x = "ANY",
  pos = "numeric",
  Bootstrap = "character",
  B = "numeric", 
  n = "numeric",
  replace = "logical", 
  prob = "numeric",
  c = "numeric",
  c.se = "numeric",
  c.cv = "numeric",
  c.mode = "numeric",
  c.prob = "numeric",
  w = "matrix",
  w.se = "numeric",
  w.cv = "numeric",
  Theta = "list",
  Theta.se = "list",
  Theta.cv = "list"),
prototype = list(pos = 1,
  Bootstrap = "parametric",
  B = 100,
  replace = TRUE))

setMethod("initialize", "REBMIX.boot", 
function(.Object, ...,
  x,
  pos,
  Bootstrap,
  B,
  n,
  replace,
  prob)
{
  model <- gsub("\\.boot", "", .Object@class[1])
  
  # x.

  if (missing(x) || (length(x) == 0)) {
    stop(sQuote("x"), " must not be empty!", call. = FALSE)
  }

  if (class(x) != model) {
    stop(sQuote("x"), " object of class ", model, " is requested!", call. = FALSE)
  }

  # pos.

  if (missing(pos) || (length(pos) == 0)) pos <- .Object@pos
  
  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }
  
  # Bootstrap.

  if (missing(Bootstrap) || (length(Bootstrap) == 0)) Bootstrap <- .Object@Bootstrap

  if (!is.character(Bootstrap)) {
    stop(sQuote("Bootstrap"), " character is requested!", call. = FALSE)
  } 
  
  Bootstrap <- match.arg(Bootstrap, .rebmix.boot$Bootstrap, several.ok = FALSE) 

  # B.

  if (missing(B) || (length(B) == 0)) B <- .Object@B
  
  if (!is.wholenumber(B)) {
    stop(sQuote("B"), " integer is requested!", call. = FALSE)
  }
  
  length(B) <- 1

  if (B < 1) {
    stop(sQuote("B"), " must be greater than 0!", call. = FALSE)
  }
  
  # n.

  if (missing(n) || (length(n) == 0)) n <- .Object@n

  nmax <- nrow(as.matrix(x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]]))
  
  if (length(n) == 0) {
    n <- nmax
  }
  else {
    if (!is.wholenumber(n)) {
      stop(sQuote("n"), " integer is requested!", call. = FALSE)
    }
  
    if ((n < 1) || (n > nmax)) {
      stop(sQuote("n"), " must be greater than 0 and less or equal than ", nmax, "!", call. = FALSE)
    }
  }

  # replace.

  if (missing(replace) || (length(replace) == 0)) replace <- .Object@replace
  
  if (!is.logical(replace)) {
    stop(sQuote("replace"), " logical is requested!", call. = FALSE)
  }

  # prob.

  if (missing(prob) || (length(prob) == 0)) {
    prob <- .Object@prob
  }
  else {
    if (!is.numeric(prob)) {
      stop(sQuote("prob"), " numeric vector is requested!", call. = FALSE)
    }

    if (length(prob) != length(n)) {
      stop("lengths of ", sQuote("prob"), " and ", sQuote("n"), " must match!", call. = FALSE)
    }
  }  
  
  callNextMethod(.Object, ...,
    x = x,
    pos = pos,
    Bootstrap = Bootstrap,
    B = B,
    n = n,
    replace = replace,
    prob = prob)
}) ## initialize

setMethod("show",
          signature(object = "REBMIX.boot"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMIX.boot is requested!", call. = FALSE)
  }
  
  cat("An object of class ", "\"", class(object), "\"", "\n", sep = "")  
  
  cat("Slot \"c\":", "\n", sep = "")

  print(object@c, quote = FALSE)
  
  cat("Slot \"c.se\":", "\n", sep = "")

  print(object@c.se, quote = FALSE)  

  cat("Slot \"c.cv\":", "\n", sep = "")

  print(object@c.cv, quote = FALSE)
  
  cat("Slot \"c.mode\":", "\n", sep = "")

  print(object@c.mode, quote = FALSE) 
  
  cat("Slot \"c.prob\":", "\n", sep = "")

  print(object@c.prob, quote = FALSE)        

  rm(list = ls())
}) ## show

## Class REBMVNORM.boot

setClass("REBMVNORM.boot", contains = "REBMIX.boot")

setMethod("show",
          signature(object = "REBMVNORM.boot"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMVNORM.boot is requested!", call. = FALSE)
  }
  
  cat("An object of class ", "\"", class(object), "\"", "\n", sep = "")  
  
  cat("Slot \"c\":", "\n", sep = "")

  print(object@c, quote = FALSE)
  
  cat("Slot \"c.se\":", "\n", sep = "")

  print(object@c.se, quote = FALSE)  

  cat("Slot \"c.cv\":", "\n", sep = "")

  print(object@c.cv, quote = FALSE)
  
  cat("Slot \"c.mode\":", "\n", sep = "")

  print(object@c.mode, quote = FALSE) 
  
  cat("Slot \"c.prob\":", "\n", sep = "")

  print(object@c.prob, quote = FALSE)        

  rm(list = ls())
}) ## show

## Class RCLRMIX

setClass("RCLRMIX",
slots = c(x = "ANY",
  pos = "numeric",
  Zt = "factor",
  Zp = "factor"),
prototype = list(pos = 1))

setMethod("initialize", "RCLRMIX", 
function(.Object, ...,
  x,
  pos,
  Zt)
{
  model <- gsub("RCLR", "REB", .Object@class[1])

  # x.

  if (missing(x) || (length(x) == 0)) {
    stop(sQuote("x"), " must not be empty!", call. = FALSE)
  }

  if (class(x) != model) {
    stop(sQuote("x"), " object of class ", model, " is requested!", call. = FALSE)
  }  

  # pos.

  if (missing(pos) || (length(pos) == 0)) pos <- .Object@pos
  
  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }
  
  # Zt.
  
  if (missing(Zt) || (length(Zt) == 0)) Zt <- .Object@Zt
  
  if (!is.factor(Zt)) {
    stop(sQuote("Zt"), " factor is requested!", call. = FALSE)
  }
  
  levels(Zt) <- 1:length(levels(Zt))  
  
  callNextMethod(.Object, ...,
    x = x,
    pos = pos,
    Zt = Zt)
}) ## initialize

setMethod("show",
          signature(object = "RCLRMIX"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLRMIX is requested!", call. = FALSE)
  }
  
  cat("An object of class ", "\"", class(object), "\"", "\n", sep = "")  
  
  cat("Slot \"Zp\":", sep = "")

  print(object@Zp, quote = FALSE)

  rm(list = ls())
}) ## show

## Class RCLRMVNORM

setClass("RCLRMVNORM", contains = "RCLRMIX")

setMethod("show",
          signature(object = "RCLRMVNORM"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLRMVNORM is requested!", call. = FALSE)
  }
  
  cat("An object of class ", "\"", class(object), "\"", "\n", sep = "")  
  
  cat("Slot \"Zp\":", sep = "")

  print(object@Zp, quote = FALSE)

  rm(list = ls())
}) ## show

## Class RCLSMIX

setClass("RCLSMIX",
slots = c(x = "list",
  Dataset = "data.frame",
  o = "numeric",
  s = "numeric",
  ntrain = "numeric",
  P = "numeric",
  ntest = "numeric",
  Zt = "factor",
  Zp = "factor",
  CM = "table",
  Accuracy = "numeric",
  Error = "numeric",
  Precission = "numeric",
  Sensitivity = "numeric",
  Specificity = "numeric"),
prototype = list(CM = table(0)))

setMethod("initialize", "RCLSMIX", 
function(.Object, ...,
  x,
  Dataset,
  Zt)
{
  model <- gsub("RCLS", "REB", .Object@class[1])

  # x.

  if (missing(x) || (length(x) == 0)) {
    stop(sQuote("x"), " must not be empty!", call. = FALSE)
  } 
  
  if (!all(unlist(lapply(x, class)) == model)) {
    stop(sQuote("x"), " list of ", model, " objects is requested!", call. = FALSE)
  }
  
  # o.
  
  .Object@o <- length(x)
  
  # s.
  
  s <- unique(unlist(lapply(x, function(x) length(x@Dataset))))
  
  if (length(s) != 1) {
    stop("lengths of ", sQuote("Dataset"), " in ", sQuote("x"), " must be equal!", call. = FALSE)
  }
  
  if (s == 1) {
    stop(sQuote("s"), " must be greater than 1!", call. = FALSE)
  }  
  
  .Object@s <- s
  
  # ntrain.
  
  ntrain <- matrix(unlist(lapply(x, function(x) lapply(x@Dataset, nrow))), ncol = s, byrow = TRUE)
  
  ntrain <- ntrain[!duplicated(ntrain),] 
  
  if (length(ntrain) != s) {
    stop(sQuote("Dataset"), " in ", sQuote("x"), " numbers of rows in data frames must be equal!", call. = FALSE)
  }
  
  .Object@ntrain <- ntrain
  
  # P.
  
  .Object@P <- ntrain / sum(ntrain)
  
  # Dataset.
  
  if (missing(Dataset) || (length(Dataset) == 0)) {
    stop(sQuote("Dataset"), " must not be empty!", call. = FALSE)
  }
  
  if (!is.data.frame(Dataset)) {
    stop(sQuote("Dataset"), " data frame is requested!", call. = FALSE)  
  }
  
  # ntest.
  
  .Object@ntest <- nrow(Dataset)
  
  # Zt.
  
  if (missing(Zt) || (length(Zt) == 0)) Zt <- .Object@Zt
  
  if (!is.factor(Zt)) {
    stop(sQuote("Zt"), " factor is requested!", call. = FALSE)
  }
  
  levels(Zt) <- 1:length(levels(Zt))
  
  callNextMethod(.Object, ...,
    x = x,
    Dataset = Dataset,
    Zt = Zt)
}) ## initialize

setMethod("show",
          signature(object = "RCLSMIX"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLSMIX is requested!", call. = FALSE)
  }
  
  cat("An object of class ", "\"", class(object), "\"", "\n", sep = "")  
  
  cat("Slot \"CM\":", sep = "")

  print(object@CM, quote = FALSE)
  
  cat("Slot \"Error\":", "\n", sep = "")

  print(object@Error, quote = FALSE)
  
  cat("Slot \"Precission\":", "\n", sep = "")
  
  names(object@Precission) <- NULL

  print(object@Precission, quote = FALSE) 
  
  cat("Slot \"Sensitivity\":", "\n", sep = "")
  
  names(object@Sensitivity) <- NULL

  print(object@Sensitivity, quote = FALSE)   
  
  cat("Slot \"Specificity\":", "\n", sep = "")
  
  names(object@Specificity) <- NULL

  print(object@Specificity, quote = FALSE)      

  rm(list = ls())
}) ## show

## Class RCLSMVNORM

setClass("RCLSMVNORM", contains = "RCLSMIX")

setMethod("show",
          signature(object = "RCLSMVNORM"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLSMVNORM is requested!", call. = FALSE)
  }
  
  cat("An object of class ", "\"", class(object), "\"", "\n", sep = "")  
  
  cat("Slot \"CM\":", sep = "")

  print(object@CM, quote = FALSE)
  
  cat("Slot \"Error\":", "\n", sep = "")

  print(object@Error, quote = FALSE)
  
  cat("Slot \"Precission\":", "\n", sep = "")
  
  names(object@Precission) <- NULL

  print(object@Precission, quote = FALSE) 
  
  cat("Slot \"Sensitivity\":", "\n", sep = "")
  
  names(object@Sensitivity) <- NULL

  print(object@Sensitivity, quote = FALSE)   
  
  cat("Slot \"Specificity\":", "\n", sep = "")
  
  names(object@Specificity) <- NULL

  print(object@Specificity, quote = FALSE)      

  rm(list = ls())
}) ## show

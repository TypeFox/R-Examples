########## R function: aspmFormRead ##########

# For reading in a formula for asp() and
# constructing a list of all information
# required for the fit.

# Last changed: 16 JUN 2006

"aspmFormRead" <-
  function (form, omit.missing = FALSE) 
{
  char.vec <- as.character(form)
  resp.name <- char.vec[2]
  resp.val <- eval(parse(text = resp.name))
  rhs <- paste(break.string(char.vec[3]), collapse = "")
  rhs <- rm.char(rhs, "\n")
  rhs <- rm.char(rhs, "\t")
  rhs <- break.string(rhs, "+")
  lin <- list()
  pen <- list()
  krige <- list()
  off.set <- NULL
  lin$names <- NULL
  lin$x <- NULL
  pen$name <- NULL
  pen$x <- NULL
  pen$adf <- list()
  pen$degree <- NULL
  pen$knots <- pen$var.knot <- list()
  pen$spar <- list()
  pen$basis <- NULL
  pen$adap <- NULL
  krige$name <- NULL
  krige$x <- NULL
  krige$adf <- NULL
  krige$knots <- krige$var.knot <- list()
  krige$spar <- NULL
  krige$bdry <- NULL
  krige$degree <- NULL
  krige$adap <- NULL
  krige$num.knots <- NULL
  if (rhs[1] == "1") {
    lin <- NULL
    pen <- NULL
    krige <- NULL
  }
  if (rhs[1] != "1") {
    for (i in 1:length(rhs)) {
      term <- rhs[i]
      type <- spmTermType(term)
      if (type == "offset") 
        off.set <- eval(parse(text = term))
      if (type == "lin") {
        lin$name <- c(lin$name, term)
        lin$x <- cbind(lin$x, eval(parse(text = term)))
      }
      if (type == "pen") {
        out <- aspmPenRead(term)
        pen$name <- c(pen$name, out$name)
        pen$x <- cbind(pen$x, out$var)
        pen$adf <- c(pen$adf, list(out$adf))
        pen$degree <- c(pen$degree, out$degree)
        pen$knots <- c(pen$knots, list(out$knots))
        pen$var.knot <- c(pen$var.knot, list(out$var.knot))
        pen$spar <- c(pen$spar, list(out$spar))
        pen$basis <- c(pen$basis, out$basis)
        pen$adap <- c(pen$adap, out$adap)
      }
      if (!is.null(pen$name)) {
        trunc.poly.basis <- sum(pen$basis == "trunc.poly")
        if (trunc.poly.basis > 0) 
          pen$basis <- "trunc.poly"
        else pen$basis <- "tps"
      }
      if (type == "krige") {
        out <- aspmKrigeRead(term)
        krige$name <- out$name
        krige$x <- out$var
        krige$adf <- out$adf
        krige$knots <- out$knots
        krige$spar <- out$spar
        krige$bdry <- out$bdry
        krige$degree <- out$degree
        krige$num.knots <- out$num.knots
        krige$var.knot <- out$var.knot
        krige$adap <- out$adap
      }
    }
  }
  if ((is.null(lin$x)) & (is.null(krige$x))) {
    if (length(pen$name) == 1) {
      if (pen$adf[[1]] == "miss") {
        term <- rhs[1]
        arg.list <- substring(term, 3, (nchar(term) - 
                                        1))
        out <- arg.search(arg.list, "df=")
        present <- out$present
        arg <- out$arg
        if (present == FALSE) 
          pen$adf[[1]] <- "miss"
        if (present == TRUE) 
          pen$adf[[1]] <- spmArgRead(arg)$val - 1
      }
    }
  }
  if (is.null(lin$x)) 
    lin <- NULL
  if (is.null(pen$x)) 
    pen <- NULL
  if (is.null(krige$x)) 
    krige <- NULL
  spm.info <- list(formula = form, y = resp.val, intercept = TRUE, 
                   lin = lin, pen = pen, krige = krige, off.set = off.set)
  datmat <- spm.info$y
  if (!is.null(spm.info$lin)) 
    datmat <- cbind(datmat, spm.info$lin$x)
  if (!is.null(spm.info$pen)) 
    datmat <- cbind(datmat, spm.info$pen$x)
  if (!is.null(spm.info$krige)) 
    datmat <- cbind(datmat, spm.info$krige$x)
  if (is.null(omit.missing)) 
    if (sum(is.na(datmat)) > 0) 
      stop("Missing data present and omit.missing not true")
  if (!is.null(omit.missing)) {
    if (omit.missing == TRUE) {
      indNA <- NULL
      for (j in 1:ncol(datmat)) indNA <- union(indNA, (1:nrow(datmat))[is.na(datmat[, 
                                                                                    j] == TRUE)])
      if (length(indNA) > 0) {
        spm.info$y <- spm.info$y[-indNA]
        if (!is.null(spm.info$lin)) 
          spm.info$lin$x <- spm.info$lin$x[-indNA, ]
        if (!is.null(spm.info$pen)) 
          spm.info$pen$x <- spm.info$pen$x[-indNA, ]
        if (!is.null(spm.info$krige)) 
          spm.info$krige$x <- spm.info$krige$x[-indNA, 
                                               ]
      }
    }
  }
  return(spm.info)
}

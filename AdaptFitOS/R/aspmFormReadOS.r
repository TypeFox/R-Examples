########## R function: aspmFormRead ##########

# For reading in a formula for asp2() and
# constructing a list of all information
# required for the fit.

"aspmFormReadOS" <-
  function (form, omit.missing = FALSE,constrasts=NULL) 
{
  char.vec <- as.character(form)
  resp.name <- char.vec[2]
  resp.val <- eval(parse(text = resp.name))
  rhs <- paste(SemiPar::break.string(char.vec[3]), collapse = "")
  rhs <- rm.char(rhs, "\n")
  rhs <- rm.char(rhs, "\t")
  rhs <- SemiPar::break.string(rhs, "+")
  lin <- list()
  constrasts <-list()
  pen <- list()
  krige <- list()
  off.set <- NULL
  lin$name <- NULL
  lin$name.orig <- NULL
  lin$x <- NULL
  pen$name <- NULL
  pen$x <- NULL
  pen$adf <- list()
  pen$degree <- list()
  pen$knots <- pen$var.knots <-pen$var.basis <-pen$var.degree <- list()
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
        if (is.factor(eval(parse(text = term)))){
          temp=model.matrix(as.formula(paste("~",term)),constrasts=contrasts)
          lin$contrasts <- c(lin$contrasts, attr(temp, "contrasts")[[term]])
          names(lin$contrasts)[[length( lin$contrasts)]]=term
          temp=temp[,-1,drop=F]
          lin$name <- c(lin$name, dimnames(temp)[[2]])
          lin$name.orig <- c(lin$name.orig, term)
          temp=as.matrix(temp)
          dimnames(temp)[[2]]=NULL
        }
        else {
          lin$name <- c(lin$name, term)
          lin$name.orig <- c(lin$name.orig, term)
          temp=as.matrix(as.data.frame(eval(parse(text = term))))
          dimnames(temp)[[2]] =NULL
        }
        lin$x <- cbind(lin$x, temp)
        rm(temp)
      }
      if (type == "pen") {
        out <- aspmPenReadOS(term)
        pen$name <- c(pen$name, out$name)
        pen$x <- cbind(pen$x, out$var)
        pen$adf <- c(pen$adf, list(out$adf))
        pen$degree <- c(pen$degree, list(out$degree))
        pen$knots <- c(pen$knots, list(out$knots))
        pen$var.knots <- c(pen$var.knots, list(out$var.knots))
        pen$var.basis <- c(pen$var.basis, list(out$var.basis))
        pen$var.degree <- c(pen$var.degree, list(out$var.degree))
        pen$spar <- c(pen$spar, list(out$spar))
        pen$basis <- c(pen$basis, out$basis)
        pen$adap <- c(pen$adap, out$adap)
      }
      if (!is.null(pen$name)) {
        trunc.poly.basis <- sum(pen$basis == "trunc.poly")
        if (trunc.poly.basis > 0) 
          pen$basis <- "trunc.poly"
#         else pen$basis <- "tps"
       else pen$basis <- pen$basis[1]
      }
      if (type == "krige") {
        out <- aspmKrigeReadOS(term)
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
        arg.list <- substring(term, 3, (nchar(term) -  1))
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
                   lin = lin, pen = pen, krige = krige, off.set = off.set, contrasts=contrasts)
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
      for (j in 1:ncol(datmat)) indNA <- union(indNA, (1:nrow(datmat))[is.na(datmat[,j] == TRUE)])
      if (length(indNA) > 0) {
        spm.info$y <- spm.info$y[-indNA]
        if (!is.null(spm.info$lin)) 
          spm.info$lin$x <- spm.info$lin$x[-indNA, ]
        if (!is.null(spm.info$pen)) 
          spm.info$pen$x <- spm.info$pen$x[-indNA, ]
        if (!is.null(spm.info$krige)) 
          spm.info$krige$x <- spm.info$krige$x[-indNA, ]
      }
    }
  }


if (!is.null(spm.info$pen$basis))
for (j in 1:length(spm.info$pen$degree))
    if (spm.info$pen$basis=="os" & length(spm.info$pen$degree[[j]])==1){
      if (((spm.info$pen$degree[[j]]+1)/2)%%1 != 0) {warning("If only degree of B-spline basis is given, degree must be chosen such that q=(p+1)/2 is an integer. Set to its default p=3, q=2."); spm.info$pen$degree[[j]] =c(3,2)}
      else spm.info$pen$degree[[j]] <- c(spm.info$pen$degree[[j]], (spm.info$pen$degree[[j]]+1)/2)
    }


  return(spm.info)
}

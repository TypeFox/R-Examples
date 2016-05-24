
# Generic function condition
condition <- function(x, ...)
  UseMethod("condition")

# default method
condition.default <- function(x, tt, ...){
  if (inherits(tt, "truthTab")) {}
  else {tt <- truthTab(tt)}
  stt <- deparse(substitute(tt))
  tt <- if (inherits(tt, "")) tt else (tt)
  x <- gsub(" ", "", as.character(x))
  x <- sub("<->", "->", x)
  cond0 <- sub("->", " -> ", x)
  cond0 <- gsub("\\+", " + ", cond0)
  splitcond <- strsplit(x, "->")
  cond.nms <- allConds <- unlist(splitcond)
  for (nm in names(tt)){
    allConds <- gsub(tolower(nm), paste("(!", nm, ")", sep = ""), allConds)
    allConds <- gsub(toupper(nm), nm, allConds)
    }
  condTbls <- lapply(allConds, function(co) sign(eval(parse(text = co), tt)))
  names(condTbls) <- cond.nms
  out <- split(condTbls, rep(seq_along(splitcond), sapply(splitcond, length)))
  out <- lapply(out, as.data.frame, optional = TRUE)
  names(out) <- lapply(out, function(condtbl) paste(names(condtbl), collapse = " -> "))
  out <- lapply(out, data.frame, n = as.vector(attr(tt, "n")),
                cases = as.vector(attr(tt, "cases")), check.names = FALSE)
  out <- lapply(out, "class<-", c("cond", "data.frame"))
  structure(lapply(out, appendInfo, tt), class = "listof")
  }

# method for class condTbl (here, x is a conTbl)
condition.condTbl <- function(x, tt, ...)
  condition.default(x[["condition"]], tt, ...)

# print method for class cond
print.cond <- function(x, digits = 3, print.table = TRUE, 
                          row.names = FALSE, ...){
  if (print.table) print.data.frame(x, row.names = row.names, ...)
  info <- attr(x, "info")
  if (length(setdiff(names(x), c("n", "cases"))) == 2){
    cat("Consistency: ", formatC(attr(x, "consistency"), format = "f", digits = digits),
        " (", info["szy"], "/", info["sy"], ")\n",
        "Coverage:    ", formatC(attr(x, "coverage"), format = "f", digits = digits),
        " (", info["szy"], "/", info["sz"], ")\n", 
        "Total no. of cases: ", info["sumf"], "\n", sep = "")
    if (!is.null(nu <- attr(x, "n.unique")) && !is.null(uc <- attr(x, "unique.coverage")))
      cat(paste(format(rep(c("Unique Coverages: ", ""), c(1, length(nu) - 1))),
                format(names(uc)), " : ", formatC(uc, format = "f", digits = digits), 
                " (", nu, "/", info["sz"], ")", sep = ""),
          sep = "\n")
    }
  else
    cat("Frequency:   ", formatC(attr(x, "freq"), format = "f", digits = digits),
        " (", info["sy"], "/", info["sumf"], ")\n", sep = "")
  }


# auxiliary function appendInfo
# extracts consistency, coverage and other features from a "cond"-objects and a truth table
# and stores them in several attributes that are assigned to this object
appendInfo <- function(x, tt){
  y <- as.logical(x[[1]])
  f <- x$n
  sy <- sum(f[y])
  sumf <- sum(f)
  if (length(setdiff(names(x), c("n", "cases"))) == 2){
    z <- as.logical(x[[2]])
    sz <- sum(f[z])
    szy <- sum(f[z & y])
    attr(x, "info") <- c(sy = sy, sz = sz, szy = szy, sumf = sumf)
    attr(x, "consistency") <- szy/sy
    attr(x, "coverage") <- szy/sz
    # unique coverages
    sConds0 <- sConds <- unlist(strsplit(names(x)[1], "\\+"))
    if (length(sConds) > 1){
      for (nm in names(tt)){
        sConds <- gsub(tolower(nm), paste("(!", nm, ")", sep = ""), sConds)
        sConds <- gsub(toupper(nm), nm, sConds)
        }
      ctbl <- lapply(sConds, function(co) sign(eval(parse(text = co), tt)))
      ctbl <- as.data.frame(ctbl)
      names(ctbl) <- sConds0
      nu <- integer(length(ctbl))
      uc <- numeric(length(ctbl))
      names(uc) <- sConds0
      for (i in seq_along(ctbl)){
        nu[i] <- sum(ctbl[i] * (rowSums(ctbl[-i]) == 0) * z * f)
        uc[i] <- nu[i] / sz
        }
      attr(x, "n.unique") <- nu  
      attr(x, "unique.coverage") <- uc
      }
    }
  else {
    attr(x, "info") <- c(sy = sy, sumf = sumf)
    attr(x, "freq") <- sy/sumf
    }
  return(x)
  }


# function group.by.outcome
group.by.outcome <- function(condlst, cases = TRUE){
  stopifnot(is.list(condlst), all(sapply(condlst, inherits, "cond")))
  outc <- sapply(strsplit(names(condlst), "-> "), "[", 2)
  outc[is.na(outc)] <- "(No outcome)"
  out <- lapply(split(condlst, outc), grbyout1, cases = cases)
  structure(out, class = "listof")
  }  
grbyout1 <- function(x, cases){
  secondCols <- do.call(cbind, lapply(x, "[", 2))
  outName <- names(secondCols)[1]
  if (!all(apply(secondCols, 1, function(col) length(unique(col)) == 1)))
    stop("Response ", outName, "is not identical in all condition tables.")
  ncases <- do.call(cbind, lapply(x, "[", "n"))
  if (!all(apply(ncases, 1, function(col) length(unique(col)) == 1)))
    stop("n is not identical in all condition tables.")
  if (cases){  
    Cases <- do.call(cbind, lapply(x, "[", "cases"))
    if (!all(apply(Cases, 1, function(col) length(unique(col)) == 1)))
      stop("Cases are not identical in all condition tables.")
    }
  do.call(cbind, c(lapply(x, "[", 1), 
                   x[[1]][c(outName, "n", if (cases) "cases")]))
  }  

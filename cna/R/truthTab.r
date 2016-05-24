# truthTab
truthTab <- function(x, frequency = NULL,
                     switch = FALSE, case.cutoff = 0){
  nm.x <- deparse(substitute(x))
  stopifnot(is.data.frame(x) || is.matrix(x),
            length(case.cutoff) == 1)
  if (inherits(x, "truthTab")){
    nofcases <- attr(x, "n")
    if (!is.null(frequency) && !is.null(nofcases) &&
        !all(frequency == nofcases))
      stop("Case Frequencies in the truth table ", nm.x, " and argument frequency are contradictory!")
    if (!is.null(nofcases))
      frequency <- nofcases
    }
  if (is.matrix(x)){
    if (is.null(colnames(x)))
      colnames(x) <- LETTERS[seq_len(ncol(x))]
    x <- as.data.frame(x)
    }
  for (i in seq_along(x)) mode(x[[i]]) <- "integer"
  if (!all(unlist(x) %in% 0:1))
    stop("x must contain only 0 and 1.")
  splitInput <- split(x, cx <- combine(x))
  tt <- do.call(rbind, lapply(splitInput, "[", 1, ))
  rownames(tt) <- NULL
  tt[switch] <- lapply(tt[switch], function(x) 1-x)
  f <- if (is.null(frequency)) sapply(splitInput, nrow)
    else as.vector(tapply(frequency, cx, sum))
  if (length(f) != nrow(tt) || any(f < 0) || any(is.na(f)))
    stop("Inadmissible frequency argument")
  # delete "rare" cases
  del.cases <- f <= case.cutoff
  if (any(del.cases)){
    message("Note: ", sum(f[del.cases]), " of ", sum(f),
            " cases are removed due to case.cutoff = " , case.cutoff, ".")
    tt <- tt[!del.cases, , drop = FALSE]
    f <- f[!del.cases]
    }
  # Identify constant factors
  is.constant <- sapply(tt, function(col) length(unique(col)) == 1)
  if (any(is.constant)) message("Note: The following factors are constant: ",
                                paste(names(tt)[is.constant], collapse = ", "))
  # Identify equivalent factors
  colpairs <- combn(length(tt), 2, FUN = function(x) tt[x], simplify = FALSE)
  pidcomp <- sapply(colpairs, function(X) all(X[[1]] == X[[2]]) || all(rowSums(X) == 1))
  if (any(pidcomp)){
    pairs <- t(combn(length(tt), 2, FUN = function(x) names(tt)[x]))[pidcomp, , drop = FALSE]
    pairnms <- apply(pairs, 1, function(p) paste(p, collapse = " and "))
    message("Note: The following pairs of factors are identical or complementary: ",
            paste(pairnms, collapse = ", "))
    }
  class(tt) <- c("truthTab", "data.frame")
  attr(tt, "n") <- f
  attr(tt, "cases") <-
    as.vector(sapply(splitInput[!del.cases],
                     function(x) paste(rownames(x), collapse = ",")))
  tt
  }
  
# print method for class truthTab
print.truthTab <- function(x, row.names = FALSE, show.cases = FALSE, ...){
  if (is.null(attr(x, "n")))
    warning("Attribute \"n\" is missing")
  if (is.null(attr(x, "cases")))
    warning("Attribute \"cases\" is missing")
  df.args <- list(x, n = attr(x, "n"))
  if (show.cases) df.args$cases <- attr(x, "cases")
  prntx <- do.call(data.frame, df.args)
  print(prntx, row.names = row.names, ...)
  cat("Total no.of.cases:", sum(attr(x, "n")), "\n")
  invisible(prntx)
  }

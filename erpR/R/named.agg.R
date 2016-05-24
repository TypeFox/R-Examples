named.agg <- function (formula, data, FUN, ..., subset, na.action = na.omit) 
{
  if (missing(formula) || !inherits(formula, "formula")) 
    stop("'formula' missing or incorrect")
  if (length(formula) != 3L) 
    stop("'formula' must have both left and right hand sides")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) 
    m$data <- as.data.frame(data)
  m$... <- m$FUN <- NULL
  m[[1L]] <- as.name("model.frame")
  if (formula[[2L]] == ".") {
    rhs <- unlist(strsplit(deparse(formula[[3L]]), " *[:+] *"))
    lhs <- sprintf("cbind(%s)", paste(setdiff(names(data), 
                                              rhs), collapse = ","))
    lhs
    m[[2L]][[2L]] <- parse(text = lhs)[[1L]]
  }
  mf <- eval(m, parent.frame())
  if (is.matrix(mf[[1L]])) {
    lhs <- as.data.frame(mf[[1L]])
    names(lhs) <- as.character(m[[2L]][[2L]])[-1L]
    myOut <- aggregate.data.frame(lhs, mf[-1L], FUN = FUN, ...)
    colnames(myOut) <- c(names(mf[-1L]), 
                        names(lhs))
                        # paste(names(lhs), deparse(substitute(FUN)), sep = "."))
  }
  else {
    myOut <- aggregate.data.frame(mf[1L], mf[-1L], FUN = FUN, ...)
    colnames(myOut) <- c(names(mf[-1L]), 
    					 strsplit(gsub("cbind\\(|\\)|\\s", "", names(mf[1L])), ",")[[1]])
                         #                    names(mf[1L])), ",")[[1]]
                         #paste(strsplit(gsub("cbind\\(|\\)|\\s", "", 
                         #                    names(mf[1L])), ",")[[1]],
                         #     deparse(substitute(FUN)), sep = "."))
  } 

  myOut
}

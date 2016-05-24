mutate_wrapper <- function(...) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1L]] <- quote(mutate)
  mc[[length(mc) + 1]] <- quote(dat)
  
  retFun <- function(dat) {
    eval(mc) %>% as.data.frame
  }
  
  preserve_attributes(retFun)
}

preserve_attributes <- function(fun) {
  force(fun)
  function(dat) {
    attOfX <- attributes(dat)
    res <- fun(dat)
    attOfRes <- attributes(res)
    attToPreserve <- names(attOfX)[!(names(attOfX) %in% names(attOfRes))]
    attributes(res) <- c(attributes(res), attributes(dat)[attToPreserve])
    res
  }
}

apply_by <- function(by, fun) {
  # by: variable names used for split
  force(fun)
  force(by)
  retFun <- function(dat) {
    stopifnot(all(by %in% names(dat)))
    out <- split(dat, dat[by]) %>% lapply(fun) %>% rbind_all
    as.data.frame(out)
  }
  preserve_attributes(retFun)
}

mapply_by <- function(by, funs) {
  # by: variable names used for split
  force(funs)
  force(by)
  retFun <- function(dat) {
    stopifnot(all(by %in% names(dat)))
    out <- mapply(function(dat, fun) fun(dat), dat = split(dat, dat[by]), fun = funs, SIMPLIFY = FALSE) %>% rbind_all
    as.data.frame(out)
  }
  preserve_attributes(retFun)
}

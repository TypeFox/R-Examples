#' @rdname spv
#' @method spv data.frame
#' @export
spv.data.frame <- function(n, design, type = c("spherical", "cuboidal", "lhs", "mlhs", "slhs", "rslhs", "custom"), 
                           formula, at = FALSE, keepfun, sample, unscaled = FALSE, ...){
  cll <- match.call()
  type <- tolower(type)
  type <- match.arg(type)
  if (missing(sample)){
    sample <- sampler(n = n, design = design, type = type, at = at, ...)
    if (!missing(keepfun)) {
      repeat{
        keep <- keepfun(sample)
        cnt <- sum(keep)
        sample <- sample[keep, ]
        if (cnt >= n) break
        rate <- cnt/n
        cat("Retained samples:",  round(cnt, digits = 2), 
            "-- Adding some more...\n")
        addsample <- sampler(n = max(ceiling((n - cnt)/rate), ceiling(n/10)), 
                             design = design, type = type, at = at, ...)
        sample <- rbind(sample, addsample)
      }
      cat("Final sample of size", nrow(sample))
    }
  }
  ndes <- nrow(design)
  n <- nrow(sample)
  m <- ncol(design)
  if (is(formula, "formula")){
    formula <- as.formula(paste("~", paste(attr(terms(formula, data = sample), "term.labels"), 
                                           collapse = " + ")))
    mat <- model.matrix(formula, data = as.data.frame(sample))
    mod.mat <- model.matrix(formula, data = design)
    p <- ncol(mod.mat)
    FtF.inv <- solve(crossprod(mod.mat))
    tmp <- .Fortran("fds", as.integer(p), as.integer(n), as.integer(ndes), 
                    as.double(FtF.inv), as.double(mat), double(n), 
                    PACKAGE = "vdg")
    spv <- tmp[[6]]
    if (unscaled) spv <- spv / ndes
    out <- list(spv = spv, sample = sample, type = type, call = cll, 
                formula = formula, at = at, FtF.inv = FtF.inv, ndes = ndes, 
                unscaled = unscaled)
    class(out) <- c("spv", "list")
    return(out)
  }
  if (is.list(formula)){
    formula <- lapply(formula, function(x) 
      as.formula(paste("~", paste(attr(terms(x, data = design[[1]]), "term.labels"), collapse = " + "))))
    nr <- length(formula)
    nms <- names(formula)
    if (is.null(nms)) nms <- paste0("formula", seq_along(formula))
    if (length(unique(nms)) != nr) stop("Formula names must be unique.")
    spvformula <- function(formula, design, sample, call, unscaled){
      ndes <- nrow(design)
      n <- nrow(sample)
      mat <- model.matrix(formula, data = as.data.frame(sample))
      m <- ncol(design)
      mod.mat <- model.matrix(formula, data = as.data.frame(design))
      p <- ncol(mod.mat)
      FtF.inv <- solve(crossprod(mod.mat))
      tmp <- .Fortran("fds", as.integer(p), as.integer(n), as.integer(ndes), 
                      as.double(FtF.inv), as.double(mat), double(n), 
                      PACKAGE = "vdg")
      spv <- tmp[[6]]
      if (unscaled) spv <- spv / ndes
      out <- list(spv = spv, sample = sample, type = type, call = call, 
                  formula = formula, at = at, FtF.inv = FtF.inv, ndes = ndes, 
                  unscaled = unscaled)
      class(out) <- c("spv", "list")
      out
    }
    cl <- makeCluster(getOption("cl.cores", min(detectCores() - 1, nr)))  
    clusterEvalQ(cl, library(vdg))
    out <- parLapply(cl, formula, spvformula, design = design, sample = sample, 
                     call = cll, unscaled = unscaled)
    stopCluster(cl)
    
    names(out) <- nms
    class(out) <- c("spvforlist", "list")
    return(out)
  }
}
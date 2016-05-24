#' @rdname spv
#' @method spv list
#' @export
spv.list <- function(n, design, type = c("spherical", "cuboidal", "lhs", "mlhs", "slhs", "rslhs", "custom"), 
                     formula, at = FALSE, keepfun, sample, unscaled = FALSE, ...){
  cll <- match.call()
  type <- tolower(type)
  type <- match.arg(type)
  nr <- length(design)
  desnms <- names(design)
  if (is.null(desnms)) desnms <- paste0("design", seq_along(design))
  if (length(unique(desnms)) != nr) stop("Design names must be unique.")
  if (missing(sample)){
    sample <- sampler(n = n, design = design[[1]], type = type, at = at, ...)
    if (!missing(keepfun)) {
      repeat{
        keep <- keepfun(sample)
        cnt <- sum(keep)
        sample <- sample[keep, ]
        if (cnt >= n) break
        rate <- cnt/n
        message("Retained samples:",  round(cnt, digits = 2), 
            "-- Adding some more...")
        addsample <- sampler(n = max(ceiling((n - cnt)/rate), ceiling(n/10)), 
                             design = design[[1]], type = type, at = at, ...)
        sample <- rbind(sample, addsample)
      }
      message("Final sample of size", nrow(sample))
    }
  }
  spvdesign <- function(design, sample, formula, call, unscaled){
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
    out <- list(spv = spv, sample = sample, type = type, call = call, at = at, 
                FtF.inv = FtF.inv, formula = formula, ndes = ndes, 
                unscaled = unscaled)
    class(out) <- c("spv", "list")
    out
  }
  cl <- makeCluster(getOption("cl.cores", min(detectCores() - 1, nr)))
  on.exit(stopCluster(cl))
  if (is(formula, "formula")){
    formula <- as.formula(paste("~", paste(attr(terms(formula, data = design[[1]]), "term.labels"), 
                                           collapse = " + ")))
    out <- parLapply(cl, design, spvdesign, sample = sample, formula = formula, 
                     call = cll, unscaled = unscaled)
    names(out) <- desnms
    class(out) <- c("spvlist", "list")
    return(out)
  }
  if (is.list(formula)){
    formula <- lapply(formula, function(x) 
      as.formula(paste("~", paste(attr(terms(x, data = design[[1]]), "term.labels"), collapse = " + "))))
    nf <- length(formula)
    fornms <- names(formula)
    if (is.null(fornms)) fornms <- paste0("formula", seq_along(formula))
    if (length(unique(fornms)) != nf) stop("Formula names must be unique.")
    out <- lapply(formula, function(y) {
      out <- parLapply(cl, design, spvdesign, sample = sample, formula = y, 
                       call = cll, unscaled = unscaled)
      names(out) <- desnms
      class(out) <- c("spvlist", "list")
      return(out)
      })
    names(out) <- fornms
    class(out) <- c("spvlistforlist", "list")
    return(out)
  }
}
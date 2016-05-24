planor2design <- function(x, ...){
      if (!"planordesign" %in% class(x)) stop("Function planor2design requires a planordesign object")
      aus <- x@design
      fn <- x@factors
      nruns <- x@nunits
      dk <- x@designkey
      if (!nruns==nrow(aus)) stop("this planordesign object cannot be forced to class design")
      attr(aus, "desnum") <- model.matrix(~.,aus)[,-1]
      ro <- 1:nruns
      attr(aus, "run.order") <- data.frame(run.no.in.std.order=ro, run.no=ro, run.no.std.rp=ro)
      attr(aus, "design.info") <- list(type="planor",
        nruns = nruns, nfactors = length(fn@levels),
        factor.names = fn@levels, nlevels = fn@fact.info$nlev,
        replications = 1,
        repeat.only = FALSE,
        randomize = FALSE,
        seed = NULL,
        response.names = NULL,
        generators = dk,
        creator = "planor")
      class(aus) <- c("design", class(aus))
      aus
}
data2design <- function(x, quantitative=rep(FALSE, ncol(x)), ...){
      xnam <- deparse(substitute(x))
      if (!("data.frame" %in% class(x) || "matrix" %in% class(x)))
           stop("function data2design requires a data.frame or a matrix")
      if ("design" %in% class(x)){
         message(xnam, " already had class design, nothing changed")
         return(x)
      }
      if (is.matrix(x)) x <- as.data.frame(x)
      nc <- ncol(x)
      nr <- nrow(x)
      fnam <- colnames(x)
      fn <- lapply(fnam, function(obj) sort(unique(x[[obj]])))      ## bug fix 10/08/15
      names(fn) <- fnam
      nlevels <- sapply(fn, length)
      if (!nc == length(quantitative)) stop("quantitative has wrong length")
      nlevels[quantitative] <- NA
      if (any(nlevels[!quantitative]>15))
          warning("Qualitative factor(s) with more than 15 levels? \nForgot to use quantitative option?")
      fn[quantitative] <- lapply(fn[quantitative], function(obj) range(obj))
      for (i in 1:nc){
         if (is.numeric(x[[i]]) && !quantitative[i]) x[[i]] <- as.factor(x[[i]])
         }
      aus <- x
      attr(aus, "desnum") <- model.matrix(~.,aus)[,-1]
      ro <- 1:nr
      attr(aus, "run.order") <- data.frame(run.no.in.std.order=ro, run.no=ro, run.no.std.rp=ro)
      attr(aus, "design.info") <- list(type="external",
        nruns = nr, nfactors = nc,
        factor.names = fn, nlevels = nlevels,
        replications = 1,
        repeat.only = FALSE,
        randomize = FALSE,
        seed = NULL,
        response.names = NULL,
        creator = "external")
      class(aus) <- c("design", class(aus))
      aus
}
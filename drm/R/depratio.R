`depratio` <-
  function (formula, data, subset, ord = 2, boot.ci = FALSE, n.boot = NULL, 
            ci.width = 0.95) 
{
  call <- match.call()
  m <- match.call(expand = FALSE)
  m$ord <- m$boot.ci <- m$n.boot <- NULL
  m$na.action <- "na.include"
  Terms <- if (missing(data)) 
    terms(formula, c("cluster", "Time"))
  else terms(formula, c("cluster", "Time"), data = data)
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  cluster <- attr(Terms, "specials")$cluster
  Time <- attr(Terms, "specials")$Time
  if (!length(cluster)) 
    stop("Cluster not specified. Use e.g. y ~ cluster(id) + Time(time)")
  if (length(cluster) > 1) 
    stop("Only one cluster allowed. Use e.g. y ~ cluster(id) + Time(time)")
  if (!length(Time)) 
    stop("Time not specified. Use e.g. y ~ cluster(id) + Time(time)")
  if (length(Time) > 1) 
    stop("Only one Time-dimension allowed. Use e.g. y ~ cluster(id) + Time(time)")
  tempc <- untangle.specials(Terms, "cluster", 1:10)
  tempt <- untangle.specials(Terms, "Time", 1:10)
  cl <- strata(m[, tempc$vars], shortlabel = TRUE)
  ords <- attr(Terms, "ord")[c(tempc$terms, tempt$terms)]
  if (any(ords > 1)) 
    stop("cluster() or Time() can not be used in an interaction")
  ords <- c(cluster, Time, attr(Terms, "response"))
  mord <- eval(parse(text = paste("order(", paste("m[,", ords, 
                       "]", collapse = ","), ")")))
  m <- m[mord, ]
  cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
  ordlabel <- unique(as.character(strata(m[, tempt$vars], shortlabel = TRUE)))
  nrep <- table(table(cluster))
  if (length(nrep) > 1) 
    stop(paste("Can handle only balanced data frame.", "For unbalanced data frames, fill NA's to obtain equal cluster size", 
               sep = "\n"))
  nrep <- as.numeric(names(nrep))
  y <- model.extract(m, "response")
  nclass <- length(levels(factor(y)))
  ymat <- matrix(as.numeric(as.factor(y)), ncol = nrep, byrow = TRUE)
  dr.est <- function(ymat, nrep, ord) {
    path <- apply(ymat, 1, paste, collapse = "")
    ylist <- lapply(seq(ncol(ymat)), function(i, ymat) ymat[, 
                                                            i], ymat = ymat)
    num <- table(rep(1:nrep, rep(nrow(ymat), nrep)), ymat)
    path <- apply(ymat, 1, paste, sep = "")
    pa <- lapply(seq(ncol(ymat) - (ord - 1)), function(i, 
                                                       ymat, ord) apply(ymat[, i:(i + (ord - 1))], 1, paste, 
                                                                        collapse = ""), ymat = ymat, ord = ord)
    n.num <- rep(seq(length(pa)), rep(ord, length(pa)))
    n.denom <- c(sapply(seq(nrep - (ord - 1)), function(i, 
                                                        ord) i:(i + (ord - 1)), ord = ord))
    margn <- lapply(seq(length(n.denom)), function(i, n.num, 
                                                   n.denom, pa, ylist) {
      if (length(grep("NA", pa[[n.num[i]]])) > 0) 
        table(ylist[[n.denom[i]]][-grep("NA", pa[[n.num[i]]])])
      else table(ylist[[n.denom[i]]])
    }, n.num = n.num, n.denom = n.denom, pa = pa, ylist = ylist)
    pa <- lapply(pa, table)
    pa <- lapply(pa, function(i) {
      if (length(grep("NA", names(i))) > 0) 
        i[-grep("NA", names(i))]
      else i
    })
    nn <- sapply(pa, sum)
    pos <- as.matrix(expand.grid(rep(list(2:nclass), ord)))
    aa <- apply(pos, 1, paste, collapse = "")
    ntmp <- c(1, 1 + cumsum(rep(ord, length(margn)/ord)))
    ntmp <- ntmp[-length(ntmp)]
    marg <- matrix(sapply(1:nrow(pos), function(j, pos, ord, ntmp)
                          sapply(ntmp, function(i, j, pos, ord)
                                 eval(parse(text = paste(paste("margn[[", i:(i + (ord - 1)),
                                              "]][(pos[", j, ",", seq(ncol(pos)),
                                              "])]", collapse = "*"), collapse = "*"))), j = j, 
                                 pos = pos, ord = ord), pos = pos, ord = ord, ntmp = ntmp),nrow=(nrep - (ord - 1)))
    joint <- matrix(sapply(pa, function(i, aa)
                           i[match(aa, names(i))], aa = aa), ncol = (nrep - (ord - 1)))
    tau <- (nn^(ord - 1) * joint)/t(marg)
    tnam <- apply(sapply(1:(dim(tau)[2]), function(i, ord, 
                                                   ordlabel) ordlabel[i:(i + ord - 1)], ordlabel = ordlabel, 
                         ord = ord), 2, paste, collapse = ",")
    rnam <- apply(pos - 1, 1, paste, collapse = "")
    if (length(unique(rnam)) == 1) 
      rnam <- ""
    dimnames(tau) <- list(paste("tau", rnam, sep = ""), tnam)
    dimnames(joint) <- list(paste("freq", rnam, sep = ""), 
                            tnam)
    list(tau = tau, freq = joint)
  }
  tau <- dr.est(ymat, nrep, ord)
  if (boot.ci) {
    if (is.null(n.boot)) 
      stop("number of bootstrap samples (n.boot) is missing")
    cat(paste("Sampling for bootstrap CI's. (n=", n.boot, 
              ")\n", sep = ""))
    drsamp <- lapply(1:n.boot, function(i, ymat, nrep, ord) {
      if (i %in% c(100 * c(1:(n.boot/100)))) 
        cat("[", i, "]\n")
      samp <- ymat[sample(1:nrow(ymat), nrow(ymat), TRUE), 
                   ]
      out <- dr.est(samp, nrep, ord)$tau
      out[is.na(out)] <- 0
      out
    }, ymat = ymat, nrep = nrep, ord = ord)
    a <- array(unlist(drsamp), dim = c(nrow(drsamp[[1]]), 
                                 ncol(drsamp[[1]]), n.boot))
    boot <- list(lcl = apply(a, c(1, 2), function(i) quantile(i, 
                   probs = (1 - ci.width)/2)), ucl = apply(a, c(1, 2), 
                                                 function(i) quantile(i, probs = 1 - (1 - ci.width)/2)))
    dimnames(boot[[1]]) <- dimnames(boot[[2]]) <- dimnames(tau[[1]])
  }
  else (boot <- NULL)
  structure(list(tau = tau[[1]], freq = tau[[2]], boot = boot, 
                 call = call), class = "depratio")
}

create.data <-
    function (nvars = c(100, 100, 100, 100, 600),
              cors = c(0.8, 0, 0.8, 0, 0),
              associations = c(0.5, 0.5, 0.3, 0.3, 0),
              firstonly = c(TRUE,FALSE, TRUE, FALSE, FALSE),
              nsamples = 100,
              censoring = "none",
              labelswapprob = 0,
              response = "timetoevent",
              basehaz = 0.2,
              logisticintercept = 0)
{
    if (labelswapprob < 0)
        stop("labelswapprob cannot be negative")
    if (labelswapprob > 0 & response == "timetoevent")
        stop("labelswapprob is only implemented for binary response")
    if (!class(nvars) %in% c("numeric","integer"))
        stop("nvars must be a numeric vector")
    if (!class(cors) %in% c("numeric","integer"))
        stop("cors must be a numeric vector")
    if (class(firstonly) != "logical")
        stop("firstonly must be a logical vector")
    if (!class(associations) %in% c("numeric","integer"))
        stop("associations must be a numeric vector")
    if (length(nvars) != length(cors) | length(nvars) != length(firstonly) |
        length(nvars) != length(associations))
        stop("nvars, cors, firstonly, and associations must all have the same length.")
    x.out <- matrix(0, ncol = sum(nvars), nrow = nsamples)
    definecors <- data.frame(start = c(1, cumsum(nvars[-length(nvars)]) +
                             1), end = cumsum(nvars), cors = cors, associations = associations,
                             num = nvars, firstonly = firstonly, row.names = letters[1:length(nvars)])
    Sigma <- matrix(0, ncol = sum(nvars), nrow = sum(nvars))
    wts <- rep(0, sum(nvars))
    for (i in 1:nrow(definecors)) {
        thisrange <- definecors[i, "start"]:definecors[i, "end"]
        Sigma[thisrange, thisrange] <- definecors[i, "cors"]
        diag(Sigma) <- 1
        x.out[, thisrange] <- mvrnorm(n = nsamples, mu = rep(0,
                                                    nvars[i]), Sigma = Sigma[thisrange, thisrange])
        if (definecors[i, "firstonly"]) {
            wts[definecors[i, "start"]] <- definecors[i, "associations"]
        }
        else {
            wts[definecors[i, "start"]:definecors[i, "end"]] <- definecors[i,
                                                                           "associations"]
        }
        varnames <- paste(letters[i], 1:nvars[i], sep = ".")
        names(wts)[definecors[i, "start"]:definecors[i, "end"]] <- varnames
    }
    names(wts) <- make.unique(names(wts))
    dimnames(Sigma) <- list(colnames = names(wts), rownames = names(wts))
    colnames(x.out) <- names(wts)
    betaX <- x.out %*% wts
    x.out <- data.frame(x.out)
    if (identical(response, "timetoevent")) {
        h = basehaz * exp(betaX[, 1])
        x.out$time <- rexp(length(h), h)
        x.out$cens <- 1
        if(class(censoring)=="numeric" | class(censoring)=="integer"){
          if(length(censoring)==2){
            censtimes <- runif(length(h),min=censoring[1],max=censoring[2])
          }else if(length(censoring)==1){
            censtimes <- rep(censoring,length(h))
          }
          x.out$cens[x.out$time>censtimes] <- 0
          x.out$time[x.out$time>censtimes] <- censtimes[x.out$time>censtimes]
        }
    }
    else if (identical(response, "binary")) {
        p <- 1/(1 + exp(-(betaX + logisticintercept)))
        x.out$outcome <- rbinom(length(p), 1, p)
        if(labelswapprob > 0){
            do.swap <- runif(length(p)) < labelswapprob
            new.outcome <- x.out$outcome
            new.outcome[x.out$outcome==1 & do.swap] <- 0
            new.outcome[x.out$outcome==0 & do.swap] <- 1
            x.out$outcome <- new.outcome
        }
        x.out$outcome <- factor(x.out$outcome)
    }
    else stop("response must be either timetoevent or binary")
    return(list(summary = definecors, associations = wts, covariance = Sigma,
                data = x.out))
}

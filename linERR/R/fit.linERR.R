fit.linERR <-
  function(formula, beta=NULL, data, ages, lag=0)
  {
    Call <- match.call()
    if (missing(data)) 
      data <- environment(formula)
    
    mf <- match.call()
    mm <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, mm)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- suppressWarnings(eval(mf, parent.frame()))
    mt <- attr(mf, "terms")
    attr(mt,"predvars") <- gsub(" ", "", attr(mt,"predvars"), fixed = TRUE)
    loglin.part <- strsplit(strsplit(x=as.character(attr(mt,"predvars")[3]), split="|", fixed = TRUE)[[1]][1], split="+", fixed=TRUE)
    lin.part <- strsplit(strsplit(x=as.character(attr(mt,"predvars")[3]), split="|", fixed = TRUE)[[1]], split="+", fixed=TRUE)
    for (i in 1:length(lin.part))
    {
      if (length(grep("strata", lin.part[[i]])) != 0) lin.part[[i]] <- lin.part[[i]][-grep("strata", lin.part[[i]], fixed=TRUE)]
      lin.part[[i]] <- gsub("\n", "", lin.part[[i]])
    }
    lin.part1 <- lin.part[[2]]
    max.exp <- length(lin.part1)
    doses <- matrix(nrow=dim(data)[1], ncol=max.exp)
    doses[, 1] <- cbind(eval(parse(text=paste0(Call$data, "$", lin.part1[1]))))
    for (i in 2:max.exp)
    {
      doses[, i] <- cbind(eval(parse(text=paste0(Call$data, "$", lin.part1[i]))))
    }
    data$cumd <- apply(doses, 1, sum)
    
    end <- round(eval(parse(text=paste0(Call$data, "$", strsplit(x=attr(mt,"predvars")[[2]], ",", fixed=TRUE)[[1]][2]))), 8)
    status <- eval(parse(text=paste0(Call$data, "$", gsub(")", "", x=strsplit(x=attr(mt, "predvars")[[2]], ",", fixed=TRUE)[[1]][3]))))
    data[, ncol(data)+1] <- end
    data[, ncol(data)+1] <- status
    data <- data[, c(1:4, ncol(data)-1, ncol(data), (5:ncol(data))[-c(ncol(data)-1, ncol(data))])]
    data[, ncol(data)] <- NULL
    data[, ncol(data)] <- NULL
    colnames(data)[5]  <- "end"
    colnames(data)[6]  <- "status"
    
    if (length(grep("strata", formula))==0) model <- coxph(eval(parse(text=attr(mt,"predvars")[[2]])) ~ cumd, data = data)
    if (length(grep("strata", formula))>0) model <- coxph(eval(parse(text=attr(mt,"predvars")[[2]])) ~ cumd + strata(eval(parse(text=formula[[3]][[3]][[3]][[2]]))), data = data)
    rsets <- as.data.frame(coxph.detail(model, riskmat=TRUE)$riskmat)
    rsets$id <- seq(1, dim(data)[1], 1)
    id <- rsets$id
    
    ages  <- as.data.frame(cbind(id, ages))
    doses <- as.data.frame(cbind(id, doses))
    
    data_2   <- as.list(data)
    data_2   <- lapply(data_2, as.numeric)
    rsets_2  <- as.list(rsets)
    rsets_2  <- lapply(rsets_2, as.numeric)
    ages_2   <- as.list(ages)
    ages_2   <- lapply(ages_2, as.numeric)
    doses_2  <- as.list(doses)
    doses_2  <- lapply(doses_2, as.numeric)
    covariates1 <- NULL
    suppressWarnings(if (is.na(as.integer(loglin.part[[1]])))
    {
      covariates1 <- sapply(paste0(Call$data, "$", loglin.part[[1]]), function(x) eval(parse(text=x)))
    })
    covars1_2 <- ifelse(loglin.part[[1]]==1, NA, apply(covariates1, 2, as.list))
    if (all(is.na(covars1_2))) covars1_2 <- NULL
    covars1_2 <- lapply(covars1_2, as.numeric)
    
    colnames(rsets)[1:(length(rsets_2)-1)] <- round(as.numeric(colnames(rsets)[1:(length(rsets_2)-1)]), 8)
    names(rsets_2)[1:(length(rsets_2)-1)] <- round(as.numeric(colnames(rsets)[1:(length(rsets_2)-1)]), 8)
    failtimes <- as.list(colnames(rsets)[1:(length(rsets_2)-1)])
    failtimes <- lapply(failtimes, as.numeric)
    suppressWarnings(if (is.null(beta)) {
      if (!is.na(as.integer(loglin.part[[1]])) & all(as.integer(loglin.part[[1]]))==1 & length(lin.part)<=2) {
        beta <- 0.1
      } 
      if (!is.na(as.integer(loglin.part[[1]])) & all(as.integer(loglin.part[[1]]))==1 & length(lin.part)>2) {
        beta <- rep(0.1, length(lin.part)-1)
      }
      if (is.na(as.integer(loglin.part[[1]])) & length(lin.part)==2){
        if (length(grep("strata", formula))==0) beta <- c(0.1, coxph(eval(parse(text=attr(mt,"predvars")[[2]])) ~ 
                                                                       as.matrix(covariates1), data = data)$coef)
        if (length(grep("strata", formula))!=0) beta <- c(0.1, coxph(eval(parse(text=attr(mt,"predvars")[[2]])) ~ 
                                                                       as.matrix(covariates1) + strata(eval(parse(text=formula[[3]][[3]][[3]][[2]]))), data = data)$coef)
      }
      if (is.na(as.integer(loglin.part[[1]])) & length(lin.part)>2){
        if (length(grep("strata", formula))==0) beta <- c(0.1, coxph(eval(parse(text=attr(mt,"predvars")[[2]])) ~ as.matrix(covariates1), data = data)$coef,
                                                          rep(0.1, length(lin.part)-2))
        if (length(grep("strata", formula))!=0) beta <- c(0.1, coxph(eval(parse(text=attr(mt,"predvars")[[2]])) ~ as.matrix(covariates1) + strata(eval(parse(text=formula[[3]][[3]][[3]][[2]]))), data = data)$coef,
                                                          rep(0.1, length(lin.part)-2))
      }
    })
    
    add.arguments <- function(f,n){
      # adds n arguments to a function f; returns that new function 
      t = paste("arg <- alist(",
                paste(sapply(1:n, function(i) paste("x",i, "=",sep="")), collapse=","),
                ")", sep="")
      formals(f) <- eval(parse(text=t))
      f
    }
    p.est <- function()
    {
      beta3  <- vector()
      for (i in 1:length(beta))
      {
        beta3[i] <- eval(parse(text = paste0("x", i)))
      }
      beta_2   <- as.list(beta3)
      beta_2   <- lapply(beta_2, as.numeric)
      suppressWarnings(if (loglin.part[[1]]==1)
      {
        ncovs1 <- 0
      }else{
        ncovs1 <- length(loglin.part[[1]])
      })
      res1 <- .Call("llhood", beta_2, data_2, rsets_2, ages_2, doses_2, covars1_2, as.integer(length(doses_2[[1]])),
                    as.integer(length(rsets_2)-1), as.integer(length(data_2)), failtimes, as.integer(ncovs1), 
                    as.integer(max.exp), as.numeric(lag))
      return(-sum(res1*is.finite(res1), na.rm=T))
    }
    p.est <- add.arguments(p.est, length(beta))
    beta <- as.list(beta)
    names(beta) <- paste0("x", seq(1, length(beta)))
    llim <- c(-1/max(data$cumd), rep(-Inf, sum(length(covars1_2), length(lin.part)-2)))
    res <- mle(p.est, start=beta, method="L-BFGS-B", lower=llim, upper=Inf)
    vcov <- vcov(res)
    aic  <- AIC(res)
    attr(res, "lowb") <- llim
    attr(res, "beta") <- beta
    attr(res, "max.exp") <- max.exp
    attr(res, "covariates1") <- covariates1
    attr(res, "data_2") <- data_2
    attr(res, "rsets_2") <- rsets_2
    attr(res, "doses_2") <- doses_2
    attr(res, "ages_2") <- ages_2
    attr(res, "call")$start <- beta
    attr(res, "call")$lower <- llim
    attr(res, "vcov") <- vcov
    attr(res, "aic") <- aic
    attr(res, "Call") <- Call
    attr(res, "llike") <- -attr(res, "details")$value
    attr(res, "deviance") <- 2*attr(res, "details")$value
    class(res) <- c("fit.linERR")
    return(res)
  }

flim <-
function(formula, data, id, obstime,
                  t.values=NULL, method="locf", lambda=NULL, art.cens=NULL) {
  data <- as.data.frame(data)
  info <- FlimFormula(formula)
  covariates <- info$covariates
  responses <- info$responses
  Cdata <- as.data.frame(data[, c(id, obstime, covariates, responses, art.cens)])
  Cdata <- Cdata[order(Cdata[, 1], Cdata[, 2]),]
  first.obs <- match(unique(Cdata[, 1]), Cdata[, 1])
  if(any(is.na(Cdata[first.obs, ]))) {
    indx <- apply(Cdata[first.obs, ], 1, function(x) any(is.na(x)))
    who <- Cdata[first.obs[which(indx)], 1]
    print(c("Check the following IDs: ", who))
    stop("There are missing values in the first observation set for some IDs")
  }
  dat.cf <- numeric()
  if(!is.null(art.cens)) {
    dat.split <- split(Cdata, Cdata[, id])
    for(j in 1:length(dat.split)) {
      datj <- dat.split[[j]]
      nj <- dim(datj)[1]
      datj$art.cens <- rep(0, nj)
      if(any(datj[, art.cens]==1)) {
        start.cens <- min(which(datj[, art.cens]==1))
        if(start.cens < nj) {
          datj[(start.cens + 1):nj, "art.cens"] <- 1
        }  
      }
      dat.split[[j]] <- datj
    }
    dat.test <- do.call("rbind", dat.split)
    dat.test[dat.test$art.cens==1, responses] <- NA
    Cdata <- dat.test
  }
  
  Cdata$obs.type <- rep(1, dim(Cdata)[1])
  times <- SetImputationTimes(Cdata, t.values)
  Hdata <- FillInGaps(Cdata, times)
  Hdata <- ImputateGaps(Hdata, times, responses, method)
  Hdata <- ExpandDataset(Hdata, times)  
  Hdata <- SetStaticVariables(Hdata, covariates)
  Hdata <- SetIncrements(Hdata, times, responses) 
  flim.fit <- ImputateMissingValues(Hdata, times, covariates,
                                    method, lambda, info) 
  model.fits <- flim.fit$model.fits
  Hdata <- flim.fit$dataset
  Hdata <- Hdata[, c(id, obstime, responses, covariates, "obs.type")]
  rownames(Hdata) <- seq(1, dim(Hdata)[1])
  
  if(!is.null(art.cens)) {
    no.obs <- tapply(Cdata[, 1], Cdata[, 1], length)
    dat.org.split <- split(Cdata, Cdata[, id])
    dat.flim.split <- split(Hdata, Hdata[, id])
    for(j in 1:length(no.obs)) {
      datj <- dat.org.split[[j]]
      dat.flim.j <- dat.flim.split[[j]]
      datj[, responses] <- dat.flim.j[1:no.obs[j], responses]
      dat.org.split[[j]] <- datj
    }
    dat.cf <- do.call("rbind", dat.org.split)
    Hdata <- data
    name.1 <- names(Hdata)
    name.2 <- paste0(responses, ".cf")
    Hdata <- cbind(Hdata, dat.cf[, responses])
    names(Hdata) <- c(name.1, name.2)
  }
  
  returnme <- list(df = Hdata,
                   dataset = Hdata,
                   fit = model.fits,
                   times = times,
                   t.values = t.values,
                   clean = Cdata,
                   covariates = covariates,
                   responses = responses,
                   method = method,
                   info = info,
                   formula = info$reg.fmlas,
                   call = match.call(),
                   lambda = lambda,
                   org = data,
                   dat.cf = dat.cf)
  class(returnme) <- "flim"
  returnme
}

ImputateMissingValues <-
function(dat, times, covariates,
                                  method, lambda, info) {
  long.res <- info$responses 
  long.inc <- paste0(info$responses, ".inc")
  reg.fmlas <- info$reg.fmlas
  if(is.null(lambda)) {
    RegFun <- function(fmla) lm(fmla, data = dataset.i, na.action = na.omit)
  } else RegFun <- function(fmla) lm.ridge(fmla, data = dataset.i,
                                           na.action = na.omit, lambda=lambda)
  savefitshere <- list()
  for(i in 1:(length(times)-1)) {
    dataset.i <- dat[dat[, 2] == times[i], ]
    if(i>=2 & method=="recursiveX") {
      dataset.i[, long.inc] <- dat[dat[, 2] == times[i+1], long.res] -
        dataset.i[, long.res]
    }
    fit.list.i <- lapply(reg.fmlas, RegFun)
    savefitshere[[i]] <- fit.list.i
    if(sum(is.na(dataset.i)) == 0) {next}  
    if(is.null(lambda)) {
      for(k in 1:length(long.res)) {
        is.na.k <- which(is.na(dataset.i[, long.inc[k]]))
        dataset.i[is.na.k, long.inc[k]] <- 
          predict(fit.list.i[[k]], newdata = dataset.i[is.na.k, ])
        Z <- dataset.i[is.na.k, ]
        if(method == "recursive") {
          id.na.next <- dat[dat[, 2] == times[i+1] & dat[, 1] %in% Z[, 1] &
                           is.na(dat[, long.res[k]]), 1]
          dat[dat[ ,2] == times[i+1] & dat[, 1] %in% id.na.next, long.res[k]] <-
            Z[Z[, 1] %in% id.na.next, long.inc[k]] +
            Z[Z[, 1] %in% id.na.next, long.res[k]]
        } else {
          dat[dat[ ,2] == times[i+1] & dat[, 1] %in% Z[, 1], long.res[k]] <-
            Z[, long.inc[k]] + Z[, long.res[k]]
        }
      }
    } else {
      for(k in 1:length(long.res)) {
        is.na.k <- which(is.na(dataset.i[, long.inc[k]]))
        Z <- dataset.i[is.na.k, ]
        Z[, long.inc[k]] <- 1
        design.mat <- model.matrix(reg.fmlas[[k]], Z)
        Z[, long.inc[k]] <- design.mat %*% coef(fit.list.i[[k]])
        
        if(method == "recursive") {
          id.na.next <- dat[dat[, 2] == times[i+1] & dat[, 1] %in% Z[, 1] &
                              is.na(dat[, long.res[k]]), 1]
          dat[dat[ ,2] == times[i+1] & dat[, 1] %in% id.na.next, long.res[k]] <-
            Z[Z[, 1] %in% id.na.next, long.inc[k]] +
            Z[Z[, 1] %in% id.na.next, long.res[k]]
        } else {
          dat[dat[ ,2] == times[i+1] & dat[, 1] %in% Z[, 1], long.res[k]] <-
            Z[, long.inc[k]] + Z[, long.res[k]]
        }
      }
    }    
  }
  returnthis <- list(dataset = dat,
                     model.fits = savefitshere)
  returnthis
}

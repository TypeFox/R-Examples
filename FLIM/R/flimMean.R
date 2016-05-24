flimMean <-
function(flimobject, response=NULL, grouping=NULL) {
  if(is.null(response)) response <- flimobject$info$responses[1]
  data <- flimobject$df
  data.obs <- flimobject$clean
  data.obs <- data.obs[data.obs[, 2] <=  max(flimobject$times), ]
  if(is.null(grouping)) {
    fid <- as.matrix(tapply(data[, response], data[, 2], mean))
    obs <- as.matrix(tapply(data.obs[, response], data.obs[, 2], mean, na.rm=TRUE))
    retme <- as.data.frame(cbind(fid, obs))
    names(retme) <- c("hypothetical", "observed")
    retme
  }
  else {
    fid <- tapply(data[, response], list(data[, 2], data[, grouping]), mean)
    obs <- tapply(data.obs[, response],
                  list(data.obs[, 2], data.obs[, grouping]), mean, na.rm=TRUE)
    retme <- as.data.frame(cbind(fid, obs))
    nl <- length(unique(data[, grouping]))
    #names(retme)[1:nl] <- paste(names(retme)[1:nl], "hyp")
    names(retme)[(nl+1):(nl*2)] <- paste(names(retme)[(nl+1):(nl*2)], "obs")
    retme
  }
}

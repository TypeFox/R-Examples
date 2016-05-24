expand.data <-
function(fu, formula, breaks, data)
{
  Nsubject <- nrow(fu)

  ## The model matrix to be merged up to the expanded data
  if (!missing(formula)) {
    covariate.frame <- data.frame(model.matrix(formula,data)[,-1])
    cov.names <- names(covariate.frame)
    covariate.frame$id <- 1:Nsubject
  }
  else {
    covariate.frame <- NULL
    cov.names <- character(0)
  }
    
  ## 1-responses for the survival intervals
  Nbreak <- length(breaks)
  ILENGTH <- function(x) {
    pmax(pmin(x[2], breaks[-1]) - pmax(x[1], breaks[-Nbreak]), 0)
  }
  well.mat <- -matrix(apply(fu[,c(1,2), drop=FALSE], 1, ILENGTH),
                      nrow=Nbreak-1)
  
  ## For the sake of the stability of the fitting procedure, each record
  ## with a 1-response is further split into separate records for each
  ## follow-up interval:  
  id.vec <- c(diag(Nbreak-1))
  well.mat <- matrix(apply(well.mat, 2, "*", id.vec), nrow=Nbreak-1)
  well.id <- rep(1:Nsubject, each=Nbreak-1)
  valid.cols <- apply(well.mat!=0, 2, any)
  well.mat <- well.mat[,valid.cols, drop=FALSE] #Remove cols that are all zero
  well.id <- well.id[valid.cols]
  
  ## 0-responses for the event intervals
  is.case <- !is.na(fu[,3])             # Observed to become ill
  fu.cases <- subset(fu, is.case)
  ill.mat <- -matrix(apply(fu.cases[,c(2,3), drop=FALSE], 1, ILENGTH),
                     nrow=Nbreak-1)
  ill.id <- which(is.case)
  
  ## The dataframe for analysis is one observation per survival interval
  ## (well.mat) and one per event interval (ill.mat):
  rates.frame <- as.data.frame(t(cbind(well.mat, ill.mat)))
  rates.names <- paste("(", breaks[-Nbreak],",",breaks[-1],")", sep="")
  names(rates.frame) <- rates.names
  rates.frame$y <- rep(c(1,0), c(ncol(well.mat), ncol(ill.mat)))
  rates.frame$id <- c(well.id, ill.id)
  
  ## Merge the covariates on to the model matrix for the baseline
  if (!is.null(covariate.frame)) {
    model.frame <- merge(rates.frame, covariate.frame, by="id")
  }
  else {
    model.frame <- rates.frame
  }
  
  model.frame[["id"]] <- NULL #Remove id variable
  return( list( rates.frame = model.frame[,rates.names, drop=FALSE],
                  cov.frame = model.frame[,  cov.names, drop=FALSE],
                          y = model.frame[,        "y", drop=TRUE ] ) )
}

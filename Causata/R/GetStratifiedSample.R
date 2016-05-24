GetStratifiedSample <- function(connect, query, stratification.variable, stratification.variable.name, stratification.value=0) {  
  
  # calculate the sample rates and the size of sample to query for
  sample.probabilities <- SamplingProbabilities(stratification.variable, stratification.value)
  limit <- Limit(query)
  # the limit must be defined, can't be null.
  stopifnot(!is.null(limit))
  sample.size.A <- as.integer(limit * sample.probabilities$A / (sample.probabilities$A + sample.probabilities$B))
  sample.size.B <- limit - sample.size.A
  
  # get Causata.Data for each strata  
  data.A <- GetSampledDataFrame(connect, query, paste("`", stratification.variable.name, "`", " = " , stratification.value, sep=""), sample.probabilities$A, sample.size.A)
  data.B <- GetSampledDataFrame(connect, query, paste("`", stratification.variable.name, "`", " <> ", stratification.value, sep=""), sample.probabilities$B, sample.size.B)
  # create a new temporary column in the data frames to store the weights.  This will be removed later.
  ssvar <- "stratified.sample.weights.98961" 
  data.A[[ssvar]] <- 1 / sample.probabilities$A
  data.B[[ssvar]] <- 1 / sample.probabilities$B
  
  # merge together into one Causata.Data object
  merged.data.frame <- rbind(data.A$df, data.B$df)
  weights <- merged.data.frame[[ssvar]] # copy weights from data frame
  merged.data.frame[[ssvar]] <- NULL # remove weights from data frame
  return(list(df=merged.data.frame, weights=weights))
}


GetSampledDataFrame <- function(connect, query, stratification.sql, sample.probability, limit) {
  # Function to get a data frame of causata data of the given size and the weights reflecting the sample.probability
  # query: query to fetch causata data
  # stratification.sql: sql clause to filter the results for this stratum
  # sample.probability: the probability with which this record is being sampled
  # limit: the number of records to fetch
  stratum.query <- query
  stratum.query$filters <- append(stratum.query$filters, stratification.sql)
  Limit(stratum.query) <- limit
  print(as.character(stratum.query))
  stratum.data <- GetData(connect, stratum.query)
  if(nrow(stratum.data) != limit){
    stop(paste("Error -- stratified sample required", limit, "records, but", nrow(stratum.data), "were returned."))
  }
  return(stratum.data)
}


SamplingProbabilities <- function (stratification.variable, stratification.value=0) {
  # Function to find the probability of sampling each stata for the given data frame
  # stratification.variable: values of this variable will determine which stratum the record belongs to
  # stratification.value: records for which the value of the stratification.variable equal this value will be
  #                       assigned to stratum A, otherwise the record will be assigned to stratum B.
    
  stratum.A.size <- sum(stratification.variable == stratification.value, na.rm=TRUE)
  stratum.B.size <- sum(stratification.variable != stratification.value, na.rm=TRUE)
  
  if (stratum.A.size == 0) {
    warning(paste("There are no records where stratification.variable == ", stratification.value, " and is not mising", sep=""))
    stratum.A.size <- 1
  }
  if (stratum.B.size == 0) {
    warning(paste("There are no records where stratification.variable != ", stratification.value, " and is not missing", sep=""))
    stratum.B.size <- 1
  }
  # calculate sample rates assuming stratum.A.size > stratum.B.size
  stratum.A.sample.rate <- sqrt(stratum.B.size / stratum.A.size)
  stratum.B.sample.rate <- 1.0
  
  # if stratum.B.size is bigger then reverse
  if (stratum.A.sample.rate > 1.0) {
    stratum.B.sample.rate <- 1.0 / stratum.A.sample.rate
    stratum.A.sample.rate <- 1.0
  }
  
  return(list(A=stratum.A.sample.rate, B=stratum.B.sample.rate))
}
ddalpha.test <- function(learn, test, ...){
  ops <- options(warn = -1)
  on.exit(options(ops))
  
  tryCatch({
    time <- system.time(
      ddalpha <- ddalpha.train(learn, ...)
    )
    C = ddalpha$dimension+1
    cc = ddalpha.classify(objects = test[,-C],ddalpha = ddalpha)
    if (is.numeric(test[,C])){
      if(is.factor(cc[[1]]) || is.character(cc[[1]])){
        cc <- unlist(lapply(cc, as.character))
        cc[cc == "Ignored"] <- NA
      }
      equal = (cc == test[,C])
    } else {
      cc <- unlist(lapply(cc, as.numeric)) 
      equal = (cc == as.numeric(test[,C]))      
    }
    
    if(!(T %in% equal) && !(F %in% equal)) 
    { return(NA)}
    error = sum(!equal,na.rm = T)/(sum(!equal,na.rm = T)+sum(equal,na.rm = T))
    return(list(error = error, correct = sum(equal,na.rm = T), incorrect = sum(!equal,na.rm = T), 
                total = length(cc)-sum(is.na(equal)), ignored = sum(is.na(equal)), n = length(cc),
                time = time[1])) 
  }
  #  tryCatch({}
  , error = function(e) {
    print ("ERROR T")
    print (e)
  }, finally = {          
  })
  return (NA)
}


ddalpha.getErrorRateCV <- function(data, numchunks = 10,  ...){
  n = nrow(data)
  numchunks = min(n, numchunks)
  chunksize = ceiling(n/numchunks)
  
  sample = seq(from = 1, by = numchunks, length.out = chunksize)   
  
  errors = 0
  total = 0
  times = c()
  
  for (i in 1:numchunks){
    sample = sample[sample<=n]
    learn = data[-sample,]
    test  = data[sample,]
    
    el = ddalpha.test(learn, test, ...)
    if(is.list(el)){
      errors = errors + el$incorrect
      total = total + el$total
      times = c(times,el$time)
    }
    
    sample = sample+1
  }
  
  return (list(errors = errors/total, time = mean(times), time_sd = sd(times)))  
}


ddalpha.getErrorRatePart <- function(data, size = 0.3, times = 10,  ...){
  
  if (!is.numeric(size) || size <=0 || size >= nrow(data)) stop("Wrong size of excluded sequences")
  
  if(size < 1)
    size = max(1, size*nrow(data)) # at least 1 point
  
  size = as.integer(size)
  
  indexes = 1:nrow(data)
  
  errors = c()
  total = 0
  time = c()
  
  for (i in 1:times){
    samp = sample(indexes, size)
    learn = data[-samp,]
    test  = data[samp,]
    
    el = ddalpha.test(learn, test, ...)
    if(is.list(el)){
      errors = c(errors,el$incorrect/el$total)
      time = c(time,el$time)
    }
  }
  
  return (list(errors = mean(errors), errors_sd = sd(errors), time = mean(time), time_sd = sd(times)))  
}
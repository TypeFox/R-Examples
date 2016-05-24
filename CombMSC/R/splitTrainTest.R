`splitTrainTest` <-
function(dat, numTrain = length(dat)-10){
  if (!inherits(dat, "ts")) stop("dat must be a time series object.")
  props = attributes(dat)$tsp
  train = ts(dat[1:numTrain], start = props[1], freq= props[3])
  test = ts(dat[(numTrain+1):length(dat)], start = (props[1]+numTrain/props[3]), freq=props[3])
   return(list(train= train, test=test))
}


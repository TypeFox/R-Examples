plus.freqdom = function(e1,e2){
  R = e1
  lags = freqdom.lags(e1)
  for (i in 1:length(lags))
    R$operators[i,,] = e1$operators[i,,] + e2$operators[i,,]
  R
}

plus.timedom = function(e1,e2){
  R = e1
  lags = freqdom.lags(e1)
  for (i in 1:length(lags))
    R$operators[i,,] = e1$operators[i,,] + e2$operators[i,,]
  R
}

#' @S3method '+' freqdom
"+.freqdom" = function (e1,e2) plus.freqdom(e1,e2)

#' @S3method '+' timedom
"+.timedom" = function (e1,e2) plus.freqdom(e1,e2)

minus.freqdom = function(e1,e2){
  R = e1
  R$lags = union(e1$lags,e2$lags)
  R$operators = array(0,c(length(R$lags),dim(R$operators)[2:3]))
  
  for (i in 1:length(R$lags)){
    lag = R$lags[i]
    i1 = which(e1$lags == lag)
    i2 = which(e2$lags == lag)
    if (sum(i1)==0)
      R$operators[i,,] = - e2$operators[i2,,]
    else if (sum(i2)==0)
      R$operators[i,,] = e1$operators[i1,,]
    else
      R$operators[i,,] = e1$operators[i1,,] - e2$operators[i2,,]
  }
  R
}

#' @S3method '-' freqdom
"-.freqdom" = function (e1,e2) minus.freqdom(e1,e2)

#' @S3method '-' timedom
"-.timedom" = function (e1,e2) minus.freqdom(e1,e2)

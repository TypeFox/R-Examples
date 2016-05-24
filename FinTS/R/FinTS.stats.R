FinTS.stats <- function(x){
#  {
#   Don't need as zoo is in depends
#    if(require('zoo'))
      Start <- min(index(x))
#    else
#      warning('requires(zoo) for Start; ',
#              ' returning 1')
#  }
  N <- sum(!is.na(x))
  Mean <- mean(x, na.rm=TRUE)
  Sd <- sd(x, na.rm=TRUE)
  {
    if(require('e1071')){
      sk <- skewness(x, na.rm=TRUE)
      kurt <- kurtosis(x, na.rm=TRUE)
    }
    else{
      warning('requires(e1071) for skewness and kurtosis;',
              ' returning NAs')
      sk <- NA
      kurt <- NA
    }
  }
  Min <- min(x, na.rm=TRUE)
  Max <- max(x,na.rm=TRUE)
#
  data.frame(Start=Start, Size=N, Mean=Mean, Standard.Deviation=Sd,
    Skewness=sk, Excess.Kurtosis=kurt, Minimum=Min,
    Maximum=Max)
}

#<<BEGIN>>
mcdatanocontrol <- function(data, type=c("V","U","VU","0"), nsv=ndvar(), nsu=ndunc(), nvariates=1, outm="each")
#ISALIAS mcdata
#DETAILS
#\samp{mcdatanocontrol} is a dangerous version of \samp{mcnode} which forces the dimension
#of data to be \samp{(nsv x nsu x nvariates)} and gives the atributes and the class
#without any control. This function is useful when your model is tested since
#it is much more quicker.
#
#--------------------------------------------
{
  dim(data) <- NULL
  data[1:(nsv*nsu*nvariates)] <- data
  dim(data) <- c(nsv,nsu,nvariates)
  class(data) <- "mcnode"
  attr(data,which="type") <- type
  attr(data,which="outm") <- outm
  }

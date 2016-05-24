# Function for standardizing regression predictors by dividing by 2 sd'
rescale <- function (x, binary.inputs="center"){
# function to rescale by subtracting the mean and dividing by 2 sd's
  if (!is.numeric(x)){
    x <- as.numeric(factor(x))
    x.obs <- x[!is.na(x)] 
  }
  x.obs <- x[!is.na(x)]
  # for binary cases
  if (length(unique(x.obs))==2){
    if (binary.inputs=="0/1"){
      x <- (x-min(x.obs))/(max(x.obs)-min(x.obs))
      return (x)
    }
    else if (binary.inputs=="-0.5,0.5"){
      return (x-0.5)
    }
    else if (binary.inputs=="center"){
      return (x-mean(x.obs))
    }
    else if (binary.inputs=="full"){
      return ((x-mean(x.obs))/(2*sd(x.obs)))
    }
  }      
  else {
    return ((x-mean(x.obs))/(2*sd(x.obs)))
  }
}

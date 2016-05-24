descriptive <- function(x)
{
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

  mode <- function(x){
    if(all(is.wholenumber(x))==TRUE){
      matrix_number = matrix(0,nrow=length(x),ncol=2)
      for(i in 1:length(x)){
        id = which(x==x[i])
        matrix_number[i,] = c(x[i],length(id))
      }
      mode = unique(x[which(as.numeric(matrix_number[,2])==max(as.numeric(matrix_number[,2])))])
      if(all(matrix_number[,2]==matrix_number[1,2])==TRUE){
        mode = NA
      }
    }
    if(all(is.wholenumber(x))==FALSE){
      histogram = hist(x)
      if(all(histogram$counts==histogram$counts[1])==TRUE) mode = NA
      else{
        id = which(histogram$counts==max(histogram$counts))
        mode = histogram$mids[id]
      }
    if(all(histogram$counts==histogram$counts[1])==TRUE) mode = NA
    }
    return(mode)
  }
  
  skewness <- function(x){
    (1/length(x) * sum((x-mean(x))^3))/(1/length(x) * sum((x-mean(x))^2))^(3/2)
  }
  
  kurtosis <- function(x){
    (1/length(x) * sum((x-mean(x))^4))/(1/length(x) * sum((x-mean(x))^2))^2 - 3 
  }
  
  statistics   <- list(mean = round(mean(x), 5),median = round(median(x), 5),
                       mode = round(mode(x), 5), variance = round(var(x), 5),
                       Skewness = round(skewness(x), 5), Kurtosis = round(kurtosis(x), 5), minimum = round(min(x), 5),
                       maximum = round(max(x), 5), n = length(x))
  return(statistics)
}
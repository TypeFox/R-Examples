uniPlot <- function(data, type = c("qqplot","histogram","box","scatter"), mfrow = NULL, ...){
  if (!is.data.frame(data) && !is.matrix(data) && !is.numeric(data)) stop(warning('Input must be one of classes \"vector\", \"data frame\" or \"matrix\"'))
  type = match.arg(type)
  
  if (is.data.frame(data) || is.matrix(data)){
    data = as.data.frame(data)
    if (nrow(data) < 2) stop(warning("Too few number of observations (n < 2)."))
    if (is.null(colnames(data))) varNames = paste("Column",1:ncol(data),sep="")
    if (!is.null(colnames(data))) varNames = colnames(data)
    
    if (is.null(mfrow)){
      nCol = ceiling(sqrt(ncol(data)))   
      nRow = ceiling(ncol(data) / nCol)
    }
    
    if (type != "box"){
      if (is.null(mfrow)) par(mfrow = c(nRow, nCol))
      else par(mfrow=mfrow)
    }
    
    if (type == "histogram"){
      for (i in 1:ncol(data)){
        hist(data[,i], xlab = varNames[i], freq=FALSE, main="", ...) 
        x <- NULL; rm(x);
        curve(dnorm(x, mean = mean(data[,i]), sd = sd(data[,i])), col="red", add=TRUE)
      }      
    }
    
    if (type == "qqplot"){
      for (i in 1:ncol(data)){
        qqnorm(data[,i], main = paste("Normal Q-Q Plot (",varNames[i],")", sep=""),...)
        qqline(data[,i])
      }      
    }
    
    if (type == "scatter"){
      if(nrow(data) == 1 || ncol(data) == 1) stop(warning("Not available for univariate input."))
      plot(data, ...)
    }
    
    if (type == "box"){
      warning("Box-Plots are based on standardized values (centered and scaled).")
      boxplot(scale(data), names = varNames, ...)
    }
  }
    
  if (is.null(ncol(data)) || is.null(nrow(data))){
    par(mfrow=c(1,1))
    data = as.numeric(data)
    if (type == "histogram"){
      hist(data, freq=FALSE, main="", ...) 
      curve(dnorm(x, mean = mean(data), sd = sd(data)), col="red", add=TRUE)
    }
    
    if (type == "qqplot"){
      qqnorm(data, ...)
      qqline(data)      
    }
    
    if (type == "box"){
      boxplot(data, ...)
    }
    
    if (type == "scatter"){
      stop(warning("Not available for univariate input."))
    }
  }
}

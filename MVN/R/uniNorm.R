uniNorm <- function(data, type = c("SW", "CVM", "Lillie", "SF", "AD"), desc = TRUE){
  if (!is.data.frame(data) && !is.matrix(data) && !is.numeric(data)) stop(warning('Input must be one of classes \"vector\", \"data frame\" or \"matrix\"'))
  type = match.arg(type)
  
  if (type == "AD") TestName = "Anderson-Darling"
  if (type == "CVM") TestName = "Cramer-von Mises"
  if (type == "Lillie") TestName = "Lilliefors (Kolmogorov-Smirnov)"
  if (type == "SW") TestName = "Shapiro-Wilk"
  if (type == "SF") TestName = "Shapiro-Francia"
  
  if (is.data.frame(data) || is.matrix(data)){ 
    varNames = colnames(data)
    dims = dim(data)
    
    if(desc){
    descriptives = describe(data)
    descriptives = descriptives[c(-1,-6,-7,-10,-13)]
    a=t(apply(data, 2, quantile))
    b=cbind(a[,2], a[,4])
    colnames(b)=c("25th", "75th")
    descriptives = cbind(descriptives, b)
    descriptives = round(cbind(descriptives[,1:6], descriptives[,9:10], descriptives[,7:8]),3)
    names(descriptives)=c("n", "Mean", "Std.Dev", "Median", "Min", "Max", "25th", "75th", "Skew", "Kurtosis")
    #name.width <- max(sapply(names(descriptives), nchar))
    #descriptives = format(descriptives, width = name.width, justify = "centre")
    }else descriptives = NULL
    
    if (is.matrix(data)){
      data = data.frame(data)
    }
    
    if (dims[2] > 1){
      if (nrow(data) < 2) stop(warning("Too few number of observations (n < 2)."))
      if (is.null(varNames)) varNames = paste("Column",1:ncol(data),sep="")
      
      res = data.frame(matrix(NA, nrow = ncol(data), ncol=4))
      colnames(res) = c("Variable", "Statistic", "p-value", "Normality")
      res[,"Variable"] = varNames
      
      if (type == "AD") res2 = apply(data, 2, ad.test)
      if (type == "CVM") res2 = apply(data, 2, cvm.test)
      if (type == "Lillie") res2 = apply(data, 2, lillie.test)
      if (type == "SW") res2 = apply(data, 2, shapiro.test)
      if (type == "SF") res2 = apply(data, 2, sf.test)
      
      res[, 2:3] = round(ldply(res2, .fun = function(x)cbind(x$stat, x$p.value))[,-1],4)
      res$Normality = ifelse(res[,3] > 0.05, "YES", "NO")
      name.width <- max(sapply(names(res), nchar))
      res = format(res, width = name.width, justify = "centre")
    }
    result = list(Descriptive = descriptives, Normality = res) 
  }
    
  if (!is.matrix(data) && !is.data.frame(data) && (is.null(nrow(data)) || is.null(ncol(data)))){ 
    if (type == "AD") res = ad.test(data)
    if (type == "CVM") res = cvm.test(data)
    if (type == "Lillie") res = lillie.test(data)
    if (type == "SW") res = shapiro.test(data)
    if (type == "SF") res = sf.test(data)
    
    
    if(desc){
      descriptives = describe(data)
      descriptives = descriptives[c(-1,-6,-7,-10,-13)]
      a = quantile(data)
      b = cbind(a[2], a[4])
      colnames(b) = c("25th", "75th")
      descriptives = cbind(descriptives, b)
      descriptives = round(cbind(descriptives[,1:6], descriptives[,9:10], descriptives[,7:8]),3)
      names(descriptives)=c("n", "Mean", "Std.Dev", "Median", "Min", "Max", "25th", "75th", "Skew", "Kurtosis")
      rownames(descriptives) = NULL
      #name.width <- max(sapply(names(descriptives), nchar))
      #descriptives = format(descriptives, width = name.width, justify = "centre")
    }else descriptives = NULL
    
    result = list(Descriptive = descriptives, Normality = res)
    
  }
  
  {
  if(is.numeric(data)) {return(result)
  names(result)[2] = paste("Normality(", TestName, ")", sep="")}
  else {
    
    #cat("\n","  ", TestName, "'s test of Normality", sep="", "\n\n")
    names(result)[1] = c("Descriptive Statistics")
    names(result)[2] = paste(TestName, "'s Normality Test", sep="")
    return(result)
  }
  }
}



labelled_output <- function(label, digits = 3){
    function(x, ...){
        cat(paste0(label, ":\n"))
        cat(round(x[[1]], getOption('lambda.digits', digits)), "\n")	
    }
}

#' @S3method print angoff
print.angoff<-function(x, ...){
	cat("Angoff's Coefficient \n")
	cat(round(x[[1]], getOption('lambda.digits', 3)))
	cat("\n\nSplit \n")
	cat(x[[2]])
}

#' @S3method print feldt
print.feldt<-labelled_output("Feldt's Coefficient")

#' @S3method print kristof 
print.kristof<-labelled_output("Kristof's Coefficient")

#' @S3method print lambda1
print.lambda1<-labelled_output("Guttman's Lambda 1 Coefficient")

#' @S3method print lambda2
print.lambda2<-labelled_output("Guttman's Lambda 2 Coefficient")

#' @S3method print lambda3
print.lambda3<-function(x, ...){
  cat("Coefficient Alpha (Guttman's Lambda 3 Coefficient) \n\n")
  
  cat("Unstandardized \n")
  cat(round(x$lambda3[[1]], getOption('lambda.digits', 3)))
  cat("\n\n")
  cat("Standardized \n")
  cat(round(x$lambda3[[2]], getOption('lambda.digits', 3)))
  cat("\n\n")
  
  if(x$items <= x$item.stats.max){
    cat("Item Statistics \n")
    print(round(x$item.stats, getOption('lambda.digits', 3)))
  }
  
}

#' @S3method print user.lambda4
print.user.lambda4<-function(x, ...){
  
  cat("User Specified Lambda 4 \n")
  cat(round(x$lambda4, getOption('lambda.digits', 3)))
  
  cat("\n\nSplit \n")
  cat(x$Split)
  cat("\n")
  
  if(!is.null(x$Item.Statistics)){
  cat("\nItem Statistics\n")
  print(round(x$Item.Statistics, getOption('lambda.digits', 3))) 
  }
  
}

#' @S3method print lambda5
print.lambda5<-labelled_output("Guttman's Lambda 5 Coefficient")

#' @S3method print lambda6
print.lambda6<-labelled_output("Guttman's Lambda 6 Coefficient")

#' @S3method print guttman
print.guttman<-function(x, ...){
  cat("Guttman's Lambda Coefficients\n\n")
  
  cat("Lambda1      ")
  cat(round(x$Lambda1, getOption('lambda.digits', 3)))
  cat("\nLambda2      ")
  cat(round(x$Lambda2, getOption('lambda.digits', 3)))
  cat("\nLambda3      ")
  cat(round(x$Lambda3, getOption('lambda.digits', 3)))
  cat("\nLambda4(max)\t")
  cat(round(x$Lambda4, getOption('lambda.digits', 3)))
  cat("\nLambda5      ")
  cat(round(x$Lambda5, getOption('lambda.digits', 3)))
  cat("\nLambda6      ")
  cat(round(x$Lambda6, getOption('lambda.digits', 3)))
}

#' @S3method print omega.tot
print.omega.tot<-labelled_output("McDonald's Omega")

#' @S3method print raju
print.raju<-labelled_output("Raju's Coefficient")

#' @S3method print cov.lambda4
print.cov.lambda4<-function(x, ...){
  if(x$method=="Hunt"){
    cat("Covariance Maximized Lambda 4 \n")
    cat("Mean ")
    cat(round(x$lambda4[[1]], getOption('lambda.digits', 3)))
    cat("\nMax  ")
    cat(round(x$lambda4[[2]], getOption('lambda.digits', 3)))
  
    if(x$show.splits==TRUE){
      cat("\n \n")
      cat("Median split \n")
      cat(x$Splits[,2])
      cat("\nMaximum split \n")
      cat(x$Splits[,3])
    }
  
    if(x$show.lambda4s==TRUE){
      cat("\n \n")
      cat("All Maximized Lambda 4 Estimates \n")
      cat(x$lambda4s)
    }
  }
    
  if(x$method=="Osburn"){
    cat("Osburn's Maximized Lambda4 \n")
    cat(round(x$l4[[1]], getOption('lambda.digits', 3)))
    cat("\n \nSplit \n")
    cat(x$Splits)
  }
}
  
#' @S3method print quant.lambda4
print.quant.lambda4<-function(x, ...){
    cat("Quantile Lambda 4 \n")
    for(i in 1:length(x$lambda4.quantile)){
      cat(names(x$lambda4.quantile)[i])
      cat(" \n")
      cat(round(x$lambda4.quantile[i], getOption('lambda.digits', 3)))
      cat("\n\n")
    }
    
    if(x$show.lambda4s==TRUE){
      cat("\n \n")
      cat("All Optimized Lambda 4 Estimates \n")
      cat(x$lambda4s)
    }    
}



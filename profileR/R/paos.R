#' Profile Analysis for One Sample with Hotelling's T-Square
#'
#' The \code{paos} function implements profile analysis for one sample using Hotelling's T-square.
#'
#' The \code{paos} function runs profile analysis for one sample based on the Hotelling's T-square test and 
#' tests the two htypothesis. First, the null hypothesis that all the ratios of the variables in the data are 
#' equal to 1. After rejecting the first hypothesis, a secondary null hypothesis that all of the ratios of the 
#' variables in the data equal to one another (not necessarily equal to 1) is tested.
#'
#' @importFrom stats pf sd
#' @param data A data matrix or data frame; rows represent individuals, columns represent variables.
#' @param scale If TRUE (default), variables are standardized by dividing their standard deviations.
#' @return A summary table is returned, listing the following two hypothesis:
#' \itemize{
#'  \item Hypothesis 1 - Ratios of the means of the variables over the hypothesized mean are equal to 1.
#' \item Hypothesis 2 - All of the ratios are equal to each other.
#'}
#' @examples
#' \dontrun{
#' data(nutrient) 
#' paos(nutrient, scale=TRUE)
#'}
#' @seealso \code{\link{cpa}}, \code{\link{pr}}
#'@export

paos <- function(data, scale=TRUE) {
  
  if(scale){
    x1 <- scale(data, center = FALSE, scale = apply(data, 2, sd, na.rm = TRUE))
    x1 <- as.matrix(x1, ncol=ncol(data))
  } else 
    {x1 <- as.matrix(data, ncol=ncol(data))}
  
  #The difference matrix for Part II
  x2 <- matrix(, nrow=nrow(data), ncol=(ncol(data)-1))
  
  for (i in 1:(ncol(x1)-1)) {
    x2[,i] <- x1[,(i+1)]-x1[,i]
  }

  
  #PART I
  mu0 <-as.matrix(rep(1,ncol(x1)),ncol=1)
  one <- as.matrix(rep(1,nrow(x1)),ncol=1)
  ident <- diag(nrow(x1))
  ybar <- t(x1)%*%one/nrow(x1)
  s <- t(x1)%*%(ident-one%*%t(one)/nrow(x1))%*%x1/(nrow(x1)-1)

#   print(mu0)
#   print(ybar)
#   print(s)

  t2 <- nrow(x1)%*%t(ybar-mu0)%*%solve(s)%*%(ybar-mu0)
  f <- (nrow(x1)-ncol(x1))*t2/ncol(x1)/(nrow(x1)-1)
  df1=ncol(x1)
  df2=nrow(x1)-ncol(x1)
  p <- 1-pf(f, df1, df2, lower.tail = TRUE, log.p = FALSE)
  result1 <- data.frame(t2=t2,F=f,df1=df1,df2=df2,p.value=p);


  #PART II
  mu0 <-as.matrix(rep(0,ncol(x2)),ncol=1)
  one <- as.matrix(rep(1,nrow(x2)),ncol=1)
  ident <- diag(nrow(x2))
  ybar <- t(x2)%*%one/nrow(x2)
  s <- t(x2)%*%(ident-one%*%t(one)/nrow(x2))%*%x2/(nrow(x2)-1)

#   print(mu0)
#   print(ybar)
#   print(s)

  t2 <- nrow(x2)%*%t(ybar-mu0)%*%solve(s)%*%(ybar-mu0)
  f <- (nrow(x2)-ncol(x2))*t2/ncol(x2)/(nrow(x2)-1)
  df1=ncol(x2)
  df2=nrow(x2)-ncol(x2)
  p <- 1-pf(f, df1, df2, lower.tail = TRUE, log.p = FALSE)
  result2 <- data.frame(t2=t2,F=f,df1=df1,df2=df2,p.value=p);

  #Final results
  result <- rbind(result1,result2)
  colnames(result) <- c("T-Squared","F","df1","df2","p-value")
  rownames(result) <- c("Ho: Ratios of the means over Mu0=1","Ho: All of the ratios are equal to each other")
  
  cat("\nProfile Analysis for One Sample with Hotelling's T-Square:\n")
  cat("\n")
  return(result)
}






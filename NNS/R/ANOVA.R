#' VN ANOVA
#'
#' Analysis of variance (ANOVA) based on lower partial moment CDFs for multiple variables.  Returns a degree of certainty the samples belong to the same population, not a p-value.
#' @param A Matrix of variables.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' A<-cbind(x,y)
#' \dontrun{VN.ANOVA(A)}
#' @export


VN.ANOVA<- function(A){

  mean.of.means <- mean(colMeans(A))
  n<- ncol(A)

  LPM_ratio = numeric(0L)
  MAD.CDF = numeric(0L)


  for (i in 1:n){
  #Continuous CDF for each variable from Mean of Means
  LPM_ratio[i] <- LPM(1,mean.of.means,A[,i])/
    (LPM(1,mean.of.means,A[,i])+UPM(1,mean.of.means,A[,i]))


  #Continuous CDF Deviation from 0.5
  MAD.CDF[i]<- abs(LPM_ratio[i] - 0.5)
  }


  Mean.MAD.CDF <- mean(MAD.CDF)


  #Certainty associated with samples
  VN.ANOVA.rho <- (0.5 - Mean.MAD.CDF)/0.5



  #Graphs
  boxplot(A, las=2, xlab= "Means", ylab="Variable",horizontal = TRUE,
           main="ANOVA", col=rainbow(n))


  #For ANOVA Visualization
  abline(v=mean.of.means,col="red",lwd=4)


  return(c("Certainty of Same Population"=VN.ANOVA.rho))

  }

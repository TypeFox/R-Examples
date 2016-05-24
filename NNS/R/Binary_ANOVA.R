#' LPM VaR
#'
#' Generates a VaR based on the Lower Partial Moment ratio
#' @param percentile The percentile for VaR
#' @param degree Degree = 0 for discrete distributions, degree = 1 for continuous distributions
#' @param x Variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' \dontrun{LPM.VaR(0.95,0,x)}
#' @export

LPM.VaR <- function(percentile,degree,x){

  x.sort <- sort(x, decreasing=FALSE)


  for (i in (1:length(x))){

    if(LPM(degree,x.sort[i],x)/(LPM(degree,x.sort[i],x)+UPM(degree,x.sort[i],x))>=(1-percentile)){

      return((x.sort[i-1]))

    }


  }}

#' UPM VaR
#'
#' Generates an upside VaR based on the Upper Partial Moment ratio
#' @param percentile The percentile for VaR
#' @param degree Degree = 0 for discrete distributions, degree = 1 for continuous distributions
#' @param x Variable
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' \dontrun{UPM.VaR(0.95,0,x)}
#' @export

UPM.VaR <- function(percentile,degree,x){

  x.sort <- sort(x, decreasing=TRUE)


  for (i in (1:length(x))){

    if(UPM(degree,x.sort[i],x)/(LPM(degree,x.sort[i],x)+UPM(degree,x.sort[i],x))>=(1-percentile)){

      return((x.sort[i-1]))

    }}}



#' VN ANOVA Binary
#'
#' Performs an analysis of variance (ANOVA) for two variables: control and treatment.  Returns the effect size of the treatment for a specified confidence interval.
#' @param control The control group sample
#' @param treatment The treatment group sample
#' @param confidence.interval The confidence interval surrounding the control mean.
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.ANOVA.bin(x,y,0.95)}
#' @export

VN.ANOVA.bin<- function(control,treatment,confidence.interval){

        mean.of.means <- mean(c(mean(control),mean(treatment)))

  #Continuous CDF for each variable from Mean of Means
        LPM_ratio.1 <- LPM(1,mean.of.means,control)/(LPM(1,mean.of.means,control)+UPM(1,mean.of.means,control))

        LPM_ratio.2 <- LPM(1,mean.of.means,treatment)/(LPM(1,mean.of.means,treatment)+UPM(1,mean.of.means,treatment))


  #Continuous CDF Deviation from 0.5
        MAD.CDF<- mean(c(abs(LPM_ratio.1 - 0.5),abs(LPM_ratio.2 - 0.5)))


  #Certainty associated with samples
        VN.ANOVA.rho <- (0.5 - MAD.CDF)/0.5

    print(c("Control Mean" = mean(control),"Treatment Mean" = mean(treatment),"Grand Mean" = mean.of.means,"Continuous CDF of Control" =LPM_ratio.1,"Continuous CDF of Treatment" = LPM_ratio.2))

  #Graphs
        boxplot(list(control,treatment), las=2, names=c("Control","Treatment"),
              xlab= "Means", horizontal = TRUE, main= "ANOVA and Effect Size",
              col=c("grey","white"),
              cex.axis= 0.75)

        #For ANOVA Visualization
        abline(v=mean.of.means,col="red",lwd=4)
            text(mean.of.means,pos=4, 2.5, "Mean of means", col = "red")

        #Upper end of CDF confidence interval for control mean
            a=UPM.VaR((confidence.interval+(1-confidence.interval)/2),1,control)
            b=UPM.VaR(.5,1,control)
        abline(v=max(a,b),
              col="green",lwd=4, lty=3)
            text(max(a,b),
              pos=4,0.75,"mu+",col="green")

        #Lower end of CDF confidence interval for control mean
            c=LPM.VaR((confidence.interval+(1-confidence.interval)/2),1,control)
            d=LPM.VaR(.5,1,control)
        abline(v=min(c,d),
              col="blue",lwd=4, lty=3)
            text(min(c,d),
              pos=2,0.75,"mu-",col="blue")

  #Effect Size Lower Bound
        Lower.Bound.Effect=min(mean(treatment)-max(a,b),0)


  #Effect Size Upper Bound
        Upper.Bound.Effect=max(mean(treatment)-min(c,d),0)

  #Certainty Statistic and Effect Size Given Confidence Interval
        return(c(Certainty.of.Same.Population=VN.ANOVA.rho,Lower.Bound.Effect=Lower.Bound.Effect,Upper.Bound.Effect=Upper.Bound.Effect))

}


#' Bootstrap Alternative to Vasicek Test.
#' @description Creates 100(1-alpha) percent bootstrap percentile confidence intervals as an alternative to the Vasicek test.
#' @import boot  
#' @param segment - The score groupings.
#' @param target - The binary target variable.
#' @param bootsamp - How many bootstrap samples to be computed. When bootsamp is too low, a warning will be produced.
#' @param grp - An integer value of how the segments are grouped.
#' @param alpha - Confidence level.
#' @param lower - Set equal to 1 if calculating a lower bound.
#' @examples
#' data("vasdata")
#' vas.level=vastol(segment=vasdata$segment,target=vasdata$response
#' ,bootsamp=500,grp=10,alpha=0.99,lower=1)
#' @references
#' [1] Glennon, D., ``Managing Model Risk in Retail Scoring," Credit Risk Analysis Division Office of the Comptroller of the Currency, September 2012.
#' @export

vastol <-function(segment,target,bootsamp,grp,alpha,lower)
{
  
  a="seg"
  b="trg"
  data1=setNames(data.frame(segment,target), c(a,b))
  
  dataset <- data.frame(Segment=integer(),
                        ConfidenceBound=double(),
                        ObservedRate=double(),
                        stringsAsFactors=FALSE) 
  
  for(j in seq(from=min(segment), to=max(segment), by=grp))
  {
    ###############################################################################
    ################ SEGMENT DATA SET       #######################################
    ###############################################################################
    
    data2=data1[which(data1$seg==j),]
    
    ###############################################################################
    ################ TOLERANCE FOR RATE     #######################################
    ###############################################################################

    Rate <- function(d, i)
    {
      d2 <- d[i,]
      return(sum(d2$trg)/length(d2$trg))
      }
      
      Rateboot = boot(data2,Rate,bootsamp)
      RateCI=boot.ci(Rateboot, conf = alpha, type ="perc")
    if(lower==1){
      bound_Rate=RateCI$percent[4]
    } else {
      bound_Rate=RateCI$percent[5]
    }
      obsRate=Rateboot$t0
      info=c(j,bound_Rate,obsRate)

    for(x in 1:ncol(dataset))
    {
      dataset[(j+1),x]=info[x]
    }

  }
  na.omit(dataset)
}

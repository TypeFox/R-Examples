
HL.ratio <- function(x, y, conf.level=0.95, alternative="two.sided", ...)

{

x <- na.omit(as.numeric(x))
y <- na.omit(as.numeric(y))

if(any(x<=0) | any(y<=0)) 
  {
  conf.int <- c(NA,NA); estimate <- NA; attr(conf.int, which="methodname") <- "Ratio of location"
  warning("Interval can not be computed \n because zero or negative values sample occured in one of the samples")
  }
   else{
logx <- log(x); logy <- log(y); active<-TRUE

lxy <- c(logx,logy)
fxy <- factor(rep(c("x","y"),c(length(x), length(y))))
dxy <- data.frame("lxy"=lxy, "fxy"=fxy)

addargs<-list(...)

addargs$formula <- as.formula("lxy ~ fxy")
addargs$data <- dxy
addargs$alternative <- alternative
addargs$conf.level <- conf.level
addargs$conf.int <- TRUE

if(is.null(addargs$distribution)){addargs$distribution <- "exact"}

 temp <- do.call(what="wilcox_test", args=addargs)
 tempCI <- confint(temp)

 conf.int <- exp(tempCI$conf.int) 
 estimate <- exp(tempCI$estimate)
 METHOD <- "Ratio of location (Hodges-Lehmann estimator)"
 attr(conf.int, which="methodname") <- METHOD
}

return(list(
conf.int=conf.int,
estimate=estimate
))


}


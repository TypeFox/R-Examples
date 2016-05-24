
HL.diff <- function(x, y, conf.level=0.95, alternative="two.sided", ...)

{

x <- na.omit(as.numeric(x))
y <- na.omit(as.numeric(y))

xy <- c(x,y)
fxy <- factor(rep(c("x","y"),c(length(x), length(y))))
dxy <- data.frame("xy"=xy, "fxy"=fxy)
addargs<-list(...)

addargs$formula <- as.formula("xy ~ fxy")
addargs$data <- dxy
addargs$alternative <- alternative
addargs$conf.level <- conf.level
addargs$conf.int <- TRUE

if(is.null(addargs$distribution)){addargs$distribution <- "exact"}

 temp <- do.call(what="wilcox_test", args=addargs)
 tempCI <- confint(temp)

 conf.int <- tempCI$conf.int 
 estimate <- tempCI$estimate
 METHOD <- "Difference in location (Hodges-Lehmann estimator)"

attr(conf.int, which="methodname")<-METHOD

return(list(
conf.int=conf.int,
estimate=estimate
))

}

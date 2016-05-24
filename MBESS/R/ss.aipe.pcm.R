ss.aipe.pcm <- function(true.variance.trend, error.variance,  variance.true.minus.estimated.trend=NULL, duration, frequency, width, conf.level=.95, trend="linear", assurance=NULL)
{

# true.variance.trend: the true variance of the individuals' slopes
# error.variance: the level 1 variance of the errors
#  variance.true.minus.estimated.trend: if  variance.true.minus.estimated.trend is used, then error.variance need not be used.  variance.true.minus.estimated.trend=error.variance/(sum of the squared C weights)


M <- frequency*duration+1



if(trend=="intercept"|trend=="Intercept"|trend==0|trend=="INTERCEPT")
{
Trend <- 0
K.p <- 1
}

if(trend=="linear"|trend=="Linear"|trend==1|trend=="LINEAR")
{
Trend <- 1
K.p <- 1/12
}
if(trend=="quadratic"|trend=="Quadratic"|trend==2|trend=="QUADRATIC")
{
Trend <- 2
K.p <- 1/720
}
if(trend=="cubic"|trend=="Cubic"|trend==3|trend=="CUBIC")
{
Trend <- 3
K.p <- 1/1000800
}


sum.c2.pm <- K.p * (factorial(M + Trend)/factorial(M - Trend - 1)) # K.p and this expression from Raudenbush and Liu (2001).
                             
if(!missing( variance.true.minus.estimated.trend))
{
if(missing(error.variance)) stop("You need to specify either 'error.variance' or ' variance.true.minus.estimated.trend'.")
if(round(error.variance/sum.c2.pm, 3) !=  variance.true.minus.estimated.trend) stop("You specified both 'error.variance' and ' variance.true.minus.estimated.trend' and they are not consistent ( variance.true.minus.estimated.trend=error.variance/(sum of appropriate squared weights)).")
}

if(missing( variance.true.minus.estimated.trend))
{
 variance.true.minus.estimated.trend <- error.variance/sum.c2.pm
}

n.i.plus.1 <- M+2 
dif <- 1
counter <- 0
past.ss <- 0
while(dif != 0)
{
counter <- counter + 1
past.ss <- c(past.ss, n.i.plus.1)

n.i <- n.i.plus.1 
nu.i <- 2*n.i-2
critival.t.i <- qt(df=nu.i, (1-(1-conf.level)/2))
n.i.plus.1 <- ceiling((8*(true.variance.trend+ variance.true.minus.estimated.trend)*critival.t.i^2)/(width^2))
dif <- n.i.plus.1 - n.i

if(counter==100)
{
n.i.plus.1 <- max(past.ss[95:100])
dif <- 0
}

}

if(is.null(assurance))
{
print("Results for expected width to be sufficiently narrow")
return(n.i.plus.1)
}

if(!is.null(assurance))
{

n.i.plus.1 <- n.i.plus.1
dif <- 1
counter <- 0
past.ss <- 0
while(dif != 0)
{
counter <- counter + 1
past.ss <- c(past.ss, n.i.plus.1)
n.i <- n.i.plus.1 
nu.i <- 2*n.i-2
critival.chi.square.i <- qchisq(df=nu.i, assurance)
variance.pi.hat.gamma <- (critival.chi.square.i*(true.variance.trend+ variance.true.minus.estimated.trend))/nu.i
n.i.plus.1 <- ceiling((8*(variance.pi.hat.gamma)*critival.t.i^2)/(width^2))
dif <- n.i.plus.1 - n.i


if(counter==100)
{
	
n.i.plus.1 <- max(past.ss[95:100])
dif <- 0
}

}	
	
print("Results for Assurance")	
return(n.i.plus.1)
}
}



"NHS" <-
function(nx, ny, X, Y, quantile=qnorm(0.95), alternative="two.sided" )
{

m<-nx
n<-ny

# Wilson Score CI for one proportion:

 WilsonScore <- function(n,Y,quant=qnorm(0.95),alternative="two.sided") {
  
  t=Y/n
   if(alternative =="two.sided") 
    {   
est.int=(Y+(quant^2)/2)/(n+(quant)^2)
w.se=((quant)*sqrt(n*t*(1-t)+(quant^2)/4))/(n+quant^2)
KI=c( est.int-w.se, est.int+w.se )
    }
  else{
   if(alternative=="less")
    {   
est.int=(Y+(quant^2)/2)/(n+(quant)^2)
w.se=((quant)*sqrt(n*t*(1-t)+(quant^2)/4))/(n+quant^2)
KI=c( 0, est.int+w.se )
    } 
  else{
   if(alternative=="greater")
    {  
est.int=(Y+(quant^2)/2)/(n+(quant)^2)
w.se=((quant)*sqrt(n*t*(1-t)+(quant^2)/4))/(n+quant^2)
KI=c( est.int-w.se , 1 )
    }
  else{stop("argument alternative misspecified")}}}

  conf.int = KI
  conf.int
 }

pX <- X/m
pY <- Y/n
estimate <- pX-pY


 if(alternative=="two.sided")
 {
  quant<-quantile
  CIX <- WilsonScore(n=m, Y=X, quant=quant, alternative="two.sided")
  lX <- CIX[1]
  uX <- CIX[2]

  CIY <- WilsonScore(n=n, Y=Y, quant=quant, alternative="two.sided") 
  lY <- CIY[1]
  uY <- CIY[2]

  conf.int <- c(estimate - quant * sqrt( (lX*(1-lX)/m) + (uY*(1-uY)/n) ) , 
                estimate + quant * sqrt( (uX*(1-uX)/m) + (lY*(1-lY)/n) ) )
  }

 if(alternative=="less")
  {
  quant<-quantile
  CIX <- WilsonScore(n=m, Y=X, quant=quant, alternative="less")
  uX <- CIX[2]

  CIY <- WilsonScore(n=n, Y=Y, quant=quant, alternative="greater") 
  lY <- CIY[1]


  conf.int <- c( -1 , 
                estimate + quant * sqrt( (uX*(1-uX)/m) + (lY*(1-lY)/n) ) )
  }


 if(alternative=="greater")
  {
  quant<- (-quantile)
  CIX <- WilsonScore(n=m, Y=X, quant=quant, alternative="greater")
  lX <- CIX[1]
  

  CIY <- WilsonScore(n=n, Y=Y, quant=quant, alternative="less") 
  uY <- CIY[2]

  conf.int <- c(estimate - quant * sqrt( (lX*(1-lX)/m) + (uY*(1-uY)/n) ) , 1 )
  }

# if(all("two.sided", "less", "greater") != alternative)
#  { stop("argument alternative mis-specified") }


list(conf.int = conf.int,
     estimate = estimate)

# end of HybScore
}


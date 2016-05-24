"Prop.or" <- 
function(x, y, conf.level=0.95, alternative="two.sided", CImethod=c("Exact","Woolf"), ...)
{

CImethod<-match.arg(CImethod)

switch(CImethod,
Woolf={

 if( is.data.frame(x) && is.data.frame(y) )
  {
   colsx<-colSums(x)
   colsy<-colSums(y)
   xa<-colsx[1]+0.5
   xb<-colsy[1]+0.5
   xc<-colsx[2]+0.5
   xd<-colsy[2]+0.5
  }
  else
   {
    if((is.numeric(x) && is.numeric(y)) && ( length(x)==2 && length(y)==2 ))
     {
     xa<-x[1]+0.5
     xb<-y[1]+0.5
     xc<-x[2]+0.5
     xd<-y[2]+0.5
     }
   else{stop("Prop.or needs two data.frames or two numeric vectors of length 2 as input")}
   }

 estimate <- (xa*xd)/(xc*xb)
 estI <- log( estimate )

 stderrlog <- sqrt( 1/xa + 1/xb + 1/xc + 1/xd )

 if(alternative=="two.sided")
  {
   zts <- qnorm(p = 1-(1-conf.level)/2 ) 
   lower <- estI - zts * stderrlog 
   upper <- estI + zts * stderrlog
  }

 if(alternative=="less")
  {
   zos <- qnorm(p = conf.level ) 
   lower <- (-Inf)
   upper <- estI + zos * stderrlog
  }

 if(alternative=="greater")
  { 
   zos <- qnorm(p = conf.level )
   lower <- estI - zos * stderrlog
   upper <- Inf
  }

 if(is.na(lower)){lower <- -Inf}
 if(is.na(upper)){upper <- Inf}

conf.int<-exp(c(lower, upper))

METHOD<-"Adjusted Woolf interval for the odds ratio"
attr(conf.int, which="methodname")<-METHOD
},

Exact={
 if( is.data.frame(x) && is.data.frame(y) )
  {
   colsx<-colSums(x)
   colsy<-colSums(y)
   tab<-rbind(colsx,colsy)
  }
  else
   {
    if((is.numeric(x) && is.numeric(y)) && ( length(x)==2 && length(y)==2 ))
     {
     tab<-rbind(x,y)
     }
   else{stop("Prop.or needs two data.frames or two numeric vectors of length 2 as input")}
   }

args<-list(...)
args$x<-tab
args$alternative<-alternative
args$conf.level<-conf.level

temp<-do.call("fisher.test", args)
conf.int<-temp$conf.int
estimate<-temp$estimate

METHOD<-"Exact confidence interval for the odds ratio"
attr(conf.int, which="methodname")<-METHOD

})



return(
list(conf.int=conf.int,
estimate=estimate)
)  
}

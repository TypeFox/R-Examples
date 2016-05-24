"Prop.diff" <-

function(x, y, conf.level=0.95, alternative="two.sided", CImethod=c("NHS","CC","AC"), ...)
{
 args<-list(...)

CImethod<-match.arg(CImethod)

 if( is.data.frame(x) && is.data.frame(y) )
  {
   colsx<-colSums(x)
   colsy<-colSums(y)
   n <- as.numeric(c( sum(colsx), sum(colsy) ))
   x <- as.numeric(c( colsx[1], colsy[1] ))
  }
  else
   {
    if((is.numeric(x) && is.numeric(y)) && ( length(x)==2 && length(y)==2 ))
     {
     n <- as.numeric(c( sum(x), sum(y) ))
     x <- as.numeric(c( x[1], y[1] ))
     }
   else{stop("Prop.test needs two data.frames or two numeric vectors of length 2 as input")}
   }

switch(CImethod,

CC={

 args$n <- n
 args$x <- x
 args$alternative <- alternative
 args$conf.level <- conf.level
 args$alternative <- alternative

 temp<- do.call("prop.test", args)
 conf.int<-temp$conf.int
 estimate<-temp$estimate[[1]]-temp$estimate[[2]]
 METHOD<-"Continuity corrected interval for the difference of proportions"
 },

AC={

   args$nx <- n[1]
   args$ny <- n[2]
   args$X <- x[1]
   args$Y <- x[2]

switch(alternative,
"two.sided"={args$quantile<-qnorm( 1-(1-conf.level)/2 ) },
"less"={args$quantile<-qnorm(conf.level)},
"greater"={args$quantile<-qnorm(1-conf.level)}
)

 args$alternative <- alternative

temp<- do.call("Add4", args)

estimate<-temp$estimate
conf.int<-temp$conf.int

METHOD<-"Agresti-Caffo interval for the difference of proportions"

},


NHS={


   args$nx <- n[1]
   args$ny <- n[2]
   args$X <- x[1]
   args$Y <- x[2]

switch(alternative,
"two.sided"={args$quantile<-qnorm( 1-(1-conf.level)/2 ) },
"less"={args$quantile<-qnorm(conf.level)},
"greater"={args$quantile<-qnorm(1-conf.level)}
)

 args$alternative <- alternative

temp<- do.call("NHS", args)

estimate<-temp$estimate
conf.int<-temp$conf.int

METHOD<-"Newcombes Hybrid Score interval for the difference of proportions"

}
)



attr(conf.int, which="methodname")<-METHOD

return(
list(conf.int=conf.int,
estimate=estimate)   
)
}


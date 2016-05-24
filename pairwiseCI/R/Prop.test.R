
"Prop.test" <-
function(x, y, alternative="two.sided", test=c("prop.test", "fisher.test"), ...)
{
args<-list(...)

test <- match.arg(test)

switch(test,
"prop.test"={
 if(( is.data.frame(x) && is.data.frame(y)) )
  {
   colsx<-colSums(x); colsy<-colSums(y)
   args$n <- as.numeric(c( sum(colsx), sum(colsy) ))
   args$x <- as.numeric(c( colsx[1], colsy[1] ))
  }
  else
   {
    if((is.numeric(x) && is.numeric(y)) && ( length(x)==2 && length(y)==2 ))
     {
     args$n <- as.numeric(c( sum(x), sum(y) ))
     args$x <- as.numeric(c( x[1], y[1] ))
     }
   else{stop("Prop.test needs two data.frames or two numeric vectors of length 2 as input")}
   }

 args$alternative <- alternative
 temp<- do.call("prop.test", args)
},
"fisher.test"={
 if( (is.data.frame(x) && is.data.frame(y)) )
  {
   colsx<-colSums(x); colsy<-colSums(y)
   args$x <- rbind(colsx,colsy)
  }
  else
   {
    if((is.numeric(x) && is.numeric(y)) && ( length(x)==2 && length(y)==2 ))
     {
     args$x <- as.numeric(rbind(x,y))
     }
   else{stop("Prop.test needs two data.frames with two columns each, or two numeric vectors of length 2 as input")}
   }
 args$alternative <- alternative
 temp<- do.call("fisher.test", args)
})
return(temp)
}





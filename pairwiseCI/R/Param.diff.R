
Param.diff <- function(x, y, conf.level=0.95, alternative="two.sided", ...)

{

addargs<-list(...)

   addargs$x <- x
   addargs$y <- y
   addargs$alternative <- alternative
   addargs$conf.level <- conf.level

 if(is.null(addargs$var.equal))
  {addargs$var.equal<-FALSE}

 if(addargs$var.equal)
  {METHOD<-"Difference of means assuming Normal distribution and equal variances"}
 else
  {METHOD<-"Difference of means assuming Normal distribution, allowing unequal variances"}

   temp <- do.call(what="t.test", args=addargs)
   conf.int <- temp$conf.int
   attr(conf.int, which="methodname")<-METHOD
   estimate <- temp$estimate[1]-temp$estimate[2]
   names(estimate) <- "difference of means"
  

return(list(
conf.int=conf.int,
estimate=estimate
))

}

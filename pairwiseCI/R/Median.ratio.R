
Median.ratio<-function(x, y, conf.level=0.95, alternative="two.sided", ...)

{

x<-as.numeric(x)
y<-as.numeric(y)

addargs<-list(...)

if(alternative=="two.sided"){conf.levelI=conf.level}
else{conf.levelI=1-(1-conf.level)*2}

ni<-c(length(x), length(y))

data<-data.frame(resp=c(x,y), trt=as.factor(rep(1:2, ni)) )

Medianratio <- function(d, i)
{
 ind1 <- i[1:ni[1]]
 ind2 <- i[-(1:ni[1])]

 x1<-d[ind1,1]
 x2<-d[ind2,1]

 median(x1)/median(x2)
}

if(is.null(addargs$R))
 {addargs$R<-999}
if(is.null(addargs$sim))
 {addargs$sim<-"ordinary"}

# stype in boot must be always == "i":
# statistics must be always HD50ratio

bootargs<-addargs
bootargs$stype<-"i"
bootargs$statistic<-Medianratio
bootargs$data<-data
bootargs$strata<-data[,2]

boot.out<-do.call("boot", bootargs)   

conf.int <- boot.ci(boot.out=boot.out, conf = conf.levelI,type =c("perc"))$perc[4:5]

if(alternative=="less")
 {conf.int[1]<-0}
else
 {if(alternative=="greater")
  {conf.int[2]<-Inf}
 }

estimate <- median(x,0.5)[[1]]/median(y,0.5)[[1]]

METHOD <- "Ratio of medians (percentile bootstrap)"
attr(conf.int, which="methodname")<-METHOD

return(list(
conf.int=conf.int,
estimate=estimate
))

}


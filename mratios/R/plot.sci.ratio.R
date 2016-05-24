"plot.sci.ratio" <-
function(x, rho0 = 1, rho0lty=2, rho0lwd=1, rho0col="black", CIvert=FALSE, CIlty = 1, CIlwd=1, CIcex=1, CIpch=16, main=NULL, ylab=NULL, xlab=NULL, sub=NULL, length=NULL, ...)
{

old.par <- par(no.readonly=TRUE)

method <- x$method
conf.int <- x$conf.int
esti <- x$estimate
compn <- x$compnames
num <- 1:length(esti)
args <- list(...)
alternative <- x$alternative
conf.level <- x$conf.level
mymai <- par("mai")

if(is.null(sub))
{
 if(x$type=="User defined")
  {sub <- paste("User defined contrasts")}
 else
  {
   if(any(c("Marcus", "McDermott", "Williams")==x$type))
    {sub <- paste(x$type, "-type contrasts for ratios")}
   else
    {

     if(method=="UmbrellaWilliams"){sub <- paste("Umbrella-protected Williams contrasts")}
     if(method=="Changepoint"){sub <- paste("Changepoint contrasts")}
     if(method=="GrandMean"){sub <- paste("Comparisons to grand mean")}
     if(method=="AVE"){sub <- paste("Comparisons to average of others")}
     if(method=="Tukey"){sub <- paste("All pairwise comparisons")}
     if(method=="Dunnett"){sub <- paste("Many-to-one comparisons")}
     if(method=="Sequen"){sub <- paste("Sequence contrasts")}
    }
  }
}


if(x$method=="Plug"){mI <- "\n(method: Plug-in)"; mcp<-"simultaneous" }

if(x$method=="Bonf"){mI <- "\n(method: Bonferroni)"; mcp<-"simultaneous"   }

if(x$method=="MtI")
 {
  if(alternative=="two.sided")
   {mI <- "\n(method: Sidak)"; mcp<-"simultaneous" }
  else
   {mI <- "\n(method: Slepian)"; mcp<-"simultaneous" }
 } 

if(x$method=="Unadj"){mI <- ""; mcp<-"unadjusted"   }


switch(alternative,
"two.sided"={
 lower <- conf.int[,1]; upper <- conf.int[,2] 
 if(any(lower=="NSD") || any(upper=="NSD")) {stop("Mean of control not significantly different from 0, no CI available")}
 lplot <- min(lower, rho0)
 uplot <- max(upper, rho0)
},

"less"={
 upper <- conf.int[,1] 
 if(any(upper=="NSD")){stop("Mean of control not significantly different from 0, no CI available")}  
 lplot <- min(rho0, esti) 
 uplot <- max(upper, rho0)
},

"greater"={
 lower <- conf.int[,1] 
 if(any(lower=="NSD") ){stop("Mean of control not significantly different from 0, no CI available")}
 lplot <- min(rho0, lower) 
 uplot <- max(esti, rho0)
})

 llplot <- lplot - 0.1*abs(lplot-uplot)
 uuplot <- uplot + 0.1*abs(lplot-uplot)
  
if(is.null(main)){
switch(alternative,
"two.sided"={ main <- paste(conf.level*100, "%",mcp,"CI (two-sided) for ratios",mI)},
"less"={ main <- paste(conf.level*100,  "%",mcp, "upper confidence limits for ratios",mI)},
"greater"={ main <- paste(conf.level*100, "%",mcp, "lower confidence limits for ratios",mI)})
 }
 
if (is.null(ylab)) {ylab=""}
if (is.null(xlab)) {xlab=""}
### produce the plot:

# vertical CI:

if(CIvert==TRUE)
{

 plot.new()
  args <- list(...)
 # the default margin size in inches
  mymai <- par("mai")

 # adjust margin under the x axis according to length of comparison names
  xwidth<- 1.5 * max(strwidth(compn, units = "inches", cex = par("cex.axis"))) 

 if (mymai[1] < xwidth) 
        mymai[1] <- xwidth
 par(mai=mymai, new=TRUE)


if(is.null(length)){
  mypin <- par("pin")
  xin <- mypin[1]
  len <- xin/(3*length(num))
}else{
if(is.numeric(length)){len<-length}else{stop("Argument length must be numeric")}
}


plot(x = num, y = esti, axes = FALSE, ylim = c(llplot, uuplot), xlim=c(1-1/3, length(num)+1/3),
 type="p", pch=CIpch, cex=CIcex,
 main=main,
 xlab="",
 ylab=ylab,
 sub=sub
 )


axis(side = 1, at = num, labels=compn, las=2, ... )
axis(side=2, ...)
box()

switch(alternative,
"two.sided"={arrows(x0=num, x1=num, y0=lower, y1=upper, lty = CIlty, lwd=CIlwd, code=3, angle=90, length=0.5/length(num))},
"less"={arrows(x0=num, x1=num, y0=llplot, y1=upper, lty = CIlty, lwd=CIlwd, code=2, angle=90, length=0.5/length(num))},
"greater"={arrows(x0=num, x1=num, y0=lower, y1=uuplot, lty = CIlty, lwd=CIlwd, code=1, angle=90, length=0.5/length(num))})

abline(h=rho0, lty=rho0lty, lwd=rho0lwd, col=rho0col)

}


# horizontal CI:


if(CIvert==FALSE)
{

 plot.new()
  args <- list(...)
 # the default margin size in inches
  mymai <- par("mai")


 # adjust margin under the x axis according to length of comparison names
  ywidth<- 1.5 * max(strwidth(compn, units = "inches", cex = par("cex.axis"))) 

 if (mymai[2] < ywidth) 
        mymai[2] <- ywidth
 par(mai=mymai, new=TRUE)

if(is.null(length)){
  mypin <- par("pin")
  yin <- mypin[2]
  len <- yin/(3*length(num))
}else{
if(is.numeric(length)){len<-length}else{stop("Argument length must be numeric")}
}

plot(y = rev(num), x = esti, axes = FALSE, xlim = c(llplot, uuplot), ylim=c(1-1/3, length(num)+1/3),
 type="p", pch=CIpch, cex=CIcex,
 main=main,
 xlab=xlab,
 ylab="",
 sub=sub
 )

axis(side = 2, at = rev(num), labels=compn, las=2, ...)
axis(side = 1, ...)
box()

switch(alternative,
"two.sided"={arrows(y0=rev(num), y1=rev(num), x0=lower,x1=upper, lty = CIlty, lwd=CIlwd, code=3, angle=90, length=len)},
"less"={arrows(y0=rev(num), y1=rev(num), x0=llplot,x1=upper, lty = CIlty, lwd=CIlwd, code=2, angle=90, length=len)},
"greater"={arrows(y0=rev(num), y1=rev(num), x0=lower,x1=uuplot, lty = CIlty, lwd=CIlwd, code=1, angle=90, length=len)})

abline(v=rho0, lty=rho0lty, lwd=rho0lwd, col=rho0col)

}


par(old.par)


}


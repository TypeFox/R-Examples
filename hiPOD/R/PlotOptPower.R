PlotOptPower <-
function(opt.design.results, save.contour=FALSE, contour.filename=NA, plot.3d=TRUE)
{
	
cost <- opt.design.results[[1]]
sample.size <- opt.design.results[[2]]
constraint.set <- opt.design.results[[3]]
scenario.set <- opt.design.results[[4]]
newdata <- opt.design.results[[5]]
	
newdata.P <- sort(unique(newdata$P))
newdata.N.p <- sort(unique(newdata$N.p))

P<-N.p<-NULL
rm(P,N.p)

pred.power <- matrix(newdata$pred.power, nrow=length(newdata.P), ncol=length(newdata.N.p))

sample.good <- subset(newdata, subset=(newdata$is.valid.design & newdata$upper.sample.good), select=c(P, N.p))
cost.good <- subset(newdata, subset=(newdata$Xmean.good), select=c(P, N.p))

sample.N.p.temp <- sort(unique(sample.good[,2]))
cost.N.p.temp <- sort(unique(cost.good[,2]))

if(nrow(sample.good) * nrow(cost.good) == 0) {
	write("There  is no valid design!", stdout())
	return()
}

for(cur.N.p in sample.N.p.temp)
{
	sample.P.min.temp <- min(subset(sample.good, subset=(N.p==cur.N.p))$P)
	sample.P.max.temp <- max(subset(sample.good, subset=(N.p==cur.N.p))$P)
	if(cur.N.p == sample.N.p.temp[1])
	{
		sample.P.min <- sample.P.min.temp
		sample.P.max <- sample.P.max.temp	
	}
	else
	{
		sample.P.min <- c(sample.P.min, sample.P.min.temp)
		sample.P.max <- c(sample.P.max, sample.P.max.temp)
	}
}


for(cur.N.p in cost.N.p.temp)
{
	cost.P.min.temp <- min(subset(cost.good, subset=(N.p==cur.N.p))$P)
	cost.P.max.temp <- max(subset(cost.good, subset=(N.p==cur.N.p))$P)
	if(cur.N.p == sample.N.p.temp[1])
	{
		cost.P.min <- cost.P.min.temp
		cost.P.max <- cost.P.max.temp	
	}
	else
	{
		cost.P.min <- c(cost.P.min, cost.P.min.temp)
		cost.P.max <- c(cost.P.max, cost.P.max.temp)	}
}

sample.N.p <- c(sample.N.p.temp, rev(sample.N.p.temp), sample.N.p.temp[1])
sample.P <- c(sample.P.min, rev(sample.P.max), sample.P.min[1])
cost.N.p <- c(cost.N.p.temp, rev(cost.N.p.temp), cost.N.p.temp[1])
cost.P <- c(cost.P.min, rev(cost.P.max), cost.P.min[1])




## plot all the constraints:
## P*N.p>=lower.sample ==> ln(N.p)>=ln(lower.sample)-ln(P)  
## P*N.p<=sample.size ==> ln(N.p)<=ln(sample.size)-ln(P)
## Xmean, P, N.p satisfies the cost function

if(save.contour==TRUE)
{
	bmp(filename=contour.filename, width=720, height=720, antialias="default")	
}


#############################
# # # # # # Color Prints!!!
#############################
# # {
# # filled.contour(x=log(newdata.N.p), y=log(newdata.P), pred.power, plot.title = title(main = "Grid Search of Optimal Power", xlab = expression("ln("*N[p]*")"), ylab = expression("ln("*P*")")), color = terrain.colors, key.title = title(main="\npower"),  plot.axes={ axis(1); axis(2);  polygon(x=log(sample.N.p), y=log(sample.P), border=2, lwd=5); polygon(x=log(cost.N.p), y=log(cost.P), border=4, lwd=5); legend(x=min(log(newdata.N.p))+0.2, y=min(log(newdata.P))+0.6, legend=c("sample size and validity constraints", "cost constraint"), col=c(2,4), lty=c(1,1), lwd=c(5,5), cex=0.8)}, ylim=c(min(log(newdata.P)), max(log(newdata.P))))
# # }

#############################
# # # # # # Grey Prints!!!
#############################
filled.contour(x=log(newdata.N.p), y=log(newdata.P), pred.power, plot.title = title(main = "Grid Search of Optimal Power", xlab = expression("ln("*N[p]*")"), ylab = expression("ln("*P*")")), col=grey.colors(n=20, start=0.9, end=0.3), levels=seq(floor(min(pred.power)*100)/100, ceiling(max(pred.power)*100)/100, length=21), key.title = title(main="\npower"),  plot.axes={ axis(1); axis(2);  polygon(x=log(sample.N.p), y=log(sample.P), border=1, lwd=5, lty=1); polygon(x=log(cost.N.p), y=log(cost.P), border=1, lwd=3, lty=2); legend(x=min(log(newdata.N.p))+0.1, y=max(log(newdata.P))-0.1, legend=c("sample size and validity constraints", "cost constraint"), lwd=c(6,3), lty=c(1,2), cex=0.8)}, ylim=c(min(log(newdata.P)), max(log(newdata.P))))


description <- bquote(paste("Scenario: VAF=", .(as.numeric(scenario.set["MAF"])), " OR=", .(round(as.numeric(scenario.set["OR"]), digits=2)), " ", epsilon, "=", .(round(as.numeric(scenario.set["error"]), digits=3)), "    Constraints: cost=$", .(round(cost)), " sample size=", .(round(sample.size))))

mtext(description, adj=0, padj=-0.5 )


if(save.contour==TRUE)
{
	dev.off()
}

library(rgl)

if(plot.3d==TRUE)
{
ln.Xmean <- matrix(log(newdata$Xmean), nrow=length(newdata.P), ncol=length(newdata.N.p))
pred.power.3d <- newdata$pred.power
## Set-up the rgl device
plot3d(log(newdata.N.p), log(newdata.P), ln.Xmean, type = "n", xlab="ln(Np)", ylab="ln(P)", zlab="ln(lambda)")

## Need a scale for pred.power to display as colours
## Here I choose 25 equally spaced colours from a palette
cols <- terrain.colors(25)

## Break pred.power into 25 equal regions
cuts <- cut(pred.power.3d, breaks = 25)

## Add in the surface, colouring by pred.power
surface3d(log(newdata.N.p), log(newdata.P), ln.Xmean, color = cols[cuts], shininess=80)

## refine the vector x, with length(x)>1
RefineVector <- function(x, num)
{
	x.origin <- x[-length(x)]
	x.diff <- diff(x)/num
	for(i in 1:(num-1))
	{
		x.temp <- x.origin+i*x.diff
		if(i==1) x.final <- rbind(x.origin, x.temp)	else x.final <- rbind(x.final, x.temp)
	}
	c(as.vector(x.final), x[length(x)])
}
for(cur.z in 1:length(sample.N.p))
{
	index.temp <- with(newdata, which(P==sample.P[cur.z] & N.p==sample.N.p[cur.z]))
	sample.points.z.temp <- newdata[index.temp, "Xmean"]
	if(cur.z==1) sample.points.z <- sample.points.z.temp
	else sample.points.z <- c(sample.points.z, sample.points.z.temp)
}
lines3d(RefineVector(log(sample.N.p),10), RefineVector(log(sample.P),10), RefineVector(log(sample.points.z),10)+0.2, color=2, alpha=0.6, lwd=8)
for(cur.z in 1:length(cost.N.p))
{
	index.temp <- with(newdata, which(P==cost.P[cur.z] & N.p==cost.N.p[cur.z]))
	cost.points.z.temp <- newdata[index.temp, "Xmean"]
	if(cur.z==1) cost.points.z <- cost.points.z.temp
	else cost.points.z <- c(cost.points.z, cost.points.z.temp)
}
lines3d(log(cost.N.p), log(cost.P), log(cost.points.z)+0.1, color=4, alpha=0.6, lwd=8)
bg3d(color="snow")
title3d(main="3D optimal power")
}

}


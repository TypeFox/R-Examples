hist.CAMAN.object <- function(x, nbreaks=NULL, mixdens=TRUE, mixdens.col="red", 
                              return.mixdens=FALSE, data.plot=NULL, singleDistr=TRUE, 
                              main="", xlab="", plotlegend=TRUE, ...){
  object <- x
	#plots an histogram and the determined distribution of the mixture model 
	if (is.null (nbreaks)) {nbreaks = min(60,max(30,object@num.obs))}
	manualData <- FALSE
	if (is.null(data.plot)){
		if (object@family == "poisson") data.plot <- object@dat[,1]/object@dat[,2]
		else if (object@family == "binomial") data.plot <- object@dat[,1]/object@dat[,3]
		else if (object@family == "gaussian") data.plot <- object@dat[,1]         
	}
	else {manualData <- TRUE}
	tmp = hist(data.plot, breaks=nbreaks, freq=FALSE, main=main, xlab=xlab,...) 
	idx = tmp$breaks
	
	if(object@family != "gaussian") idx = round(idx)
	
	x = rep(0,length(idx))
	for(i in 1:object@num.k){ 
		if (object@family == "gaussian" && class(object)[[1]] == "CAMAN.object") {
			x = x + object@p[i] * dnorm(idx, mean=object@t[i], sd=sqrt(object@component.var))
			if(singleDistr) lines(idx, object@p[i] * dnorm(idx, mean=object@t[i], sd=sqrt(object@component.var)), lwd=1, lty=1)
		}
		else if (object@family == "poisson"){ 
			x = x + object@p[i] * dpois(idx,object@t[i])
			if(singleDistr) lines(idx, object@p[i] * dpois(idx,object@t[i]), lwd=1, lty=2)
		}
		else if (object@family == "gaussian" && class(object)[[1]] == "CAMAN.glm.object"){ 
			x = x + object@p[i] * dnorm(idx, mean=object@t[i], sd=sqrt(object@residVar))
			if(singleDistr) lines(idx, object@p[i] * dnorm(idx, mean=object@t[i], sd=sqrt(object@residVar)), lwd=1, lty=2)
		}		
	}
	lines(idx, x, col=mixdens.col, lwd=2)
	if (return.mixdens || manualData) return(x)
	if (plotlegend){
	   if(singleDistr)	legend("topright",c("mixture density", "single components"),lty=c(1,1),col=c(mixdens.col,"black"), lwd=c(2,1))
	   else legend("topright",c("mixture density"),lty=c(1),col=c(mixdens.col), lwd=c(2))
    }
}


plot.CAMAN.BIEM.object <- function(x, ellipse = TRUE, ...){
object<-x
if (!ellipse){
plot(object@Mat[,1] ,object@Mat[,2] , xlab = "x1", ylab = "x2",pch=19,col="blue",cex=0.4,main="scatter plot");}
if (ellipse){
colors <- c("red","green", "blue","brown","yellow","black")
plot(object@Mat, xlab = "x1", ylab = "x2", type = "n",main="scatter plot with ellipse")
for (a in 1:length(object@RESULT[,1])){
points(object@Mat[object@Mat[,3] == a, ], col = colors[a], pch = 19, cex = 0.4)
points(object@Z[,,a],col = colors[a],cex=0.4, pch=19)}
}
return(invisible(NULL))
}

plot.CAMAN.BIMIXALG.object <- function(x, ellipse = TRUE, ...){
object<-x
if (!ellipse){
  plot(object@Mat[,1] ,object@Mat[,2] , xlab = "x1", ylab = "x2",pch=19,col="blue",cex=0.4,main="scatter plot");
  }
if (ellipse){
  colors <- c("red","green", "blue","brown","yellow","black")
  plot(object@Mat, xlab = "x1", ylab = "x2", type = "n",main="scatter plot with ellipse")
  for (a in 1:length(object@RESULT[,1])){
    points(object@Mat[object@Mat[,3] == a, ], col = colors[a], pch = 19, cex = 0.4)
    points(object@Z[,,a],col = colors[a],cex=0.4, pch=19)}
  }
}

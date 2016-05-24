local(envir=.PBSmodEnv,expr={
cVor <- function() {
	getWinVal(scope="L"); resetGraph();
	x <- switch(xdis,
			runif(n=n,min=a1x,max=a2x),
			rnorm(n=n,mean=a1x,sd=a2x),
			rgamma(n=n,shape=a1x,rate=a2x),
			rlnorm(n=n,meanlog=a1x,sdlog=a2x),
			rlogis(n=n,location=a1x,scale=a2x),
			rpois(n=n,lambda=a1x) );
	y <- switch(ydis,
			runif(n=n,min=a1y,max=a2y),
			rnorm(n=n,mean=a1y,sd=a2y),
			rgamma(n=n,shape=a1y,rate=a2y),
			rlnorm(n=n,meanlog=a1y,sdlog=a2y),
			rlogis(n=n,location=a1y,scale=a2y),
			rpois(n=n,lambda=a2y) );
	events   <- as.EventData(data.frame(EID=1:n,X=x,Y=y),projection=1);
	polys    <- calcVoronoi(events);
	polyData <- calcArea(polys)
	names(polyData)[is.element(names(polyData), "area")] <- "Z"
	colSeq   <- c("navy","blue","skyblue","lightblue1");
	brks     <- quantile(polyData$Z,c(0,.25,.5,.75,1));
	if (length(brks)!=length(unique(brks))) {
		frame(); N<<-N+1; addLabel(.5,.5,"TRY AGAIN",col=4-N%%3,font=8,cex=3); }
	polyData <- makeProps(polyData, breaks=brks,propName="col", propVals=colSeq)

	resetGraph(); expandGraph(mar=c(2.5,2.5,1,1));
	plotMap(polys, polyProps=polyData, plt=NULL)   #--- plot the tesselation
	addPoints(events, pch=20, col="orangered")     #--- plot the points
}

require(PBSmodelling); require(PBSmapping)
if(!require(deldir,quietly=TRUE)) {
	if (getYes("Load package `deldir` from CRAN?","Package Needed")) install.packages("deldir")
	else stop("Package `deldir` needed for this example",call.=FALSE) }
createWin("CalcVorWin.txt"); N<-0;

}) # end local scope


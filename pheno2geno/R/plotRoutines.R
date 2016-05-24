#
# plottingRoutines.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
# first written Nov, 2011
# Contains: plotParentalExpression, plotChildrenExpression, plotMarkerDistribution
#           projectOldMarkers
#

# plotParentalExpression
#
# DESCRIPTION:
#  plot red points for expression values for parent of type 0, blue for parent 1 and green lines for means of rows
# PARAMETERS:
#   - population - Ril type object, must contain founders phenotypic data.
#   - markers - markers to be printed numbers or names 
#   - groupLabels - Specify which column of founders data belongs to group 0 and which to group 1
# OUTPUT:
#  Plot
#
plotParentalExpression <- function(population, markers=1:100, groupLabels=c(0,0,1,1)){
  check.population(population) # Check
	#*******remove too short chromosomes*******
	if(nrow(population$founders$phenotypes)<length(markers)){
		cat("WARNING: There are only",nrow(population$founders$phenotypes),"markers in population$founders$phenotypes, function will plot them all.\n")
		markers <- 1:nrow(population$founders$phenotypes)
	}
	for(i in markers) if(!(i %in% 1:nrow(population$founders$phenotypes))){
		stop("ERROR: There is no marker number: ",i,", stopping \n")
	}
	founders <- population$founders$phenotypes[markers,]
	plot(x=markers[1], y=founders[1,1], xlim=c(min(markers),max(markers)), ylim=c(min(founders),max(founders)), col="red",
	xlab="Marker", ylab="Expression value", main="Parental gene expression data")
	for(i in 1:nrow(founders)){
		for(j in which(groupLabels==0)){
			if(i%%50==0) cat(i,j,founders[i,j],"\n")
			points(x=markers[i],y=founders[i,j],col="red")
		}
		for(k in which(groupLabels==1)){
			if(i%%50==0) cat(i,k,founders[i,k],"\n")
			points(x=markers[i],y=founders[i,k],col="blue")
		}
	}
	points(x=markers,y=apply(founders,1,mean),col="green", pch=95, cex=3)
}

############################################################################################################
#									*** plotChildrenExpression ***
#
# DESCRIPTION:
#	boxplot of data for selected markers + points of founders mean for each marker
# 
# PARAMETERS:
# 	population - Ril type object, must contain founders phenotypic data.
# 	markers - markers to be printed numbers or names 
# 
# OUTPUT:
#	plot
#
############################################################################################################
plotChildrenExpression <- function(population, markers=1:100){
	### checks
	check.population(population)
	
	### function itself
	if(nrow(population$offspring$phenotypes)<length(markers)){
		cat("WARNING: There are only",nrow(population$offspring$phenotypes),"markers in population$offspring$phenotypes, function will plot them all.\n")
		markers <- 1:nrow(population$offspring$phenotypes)
	}
	if(nrow(population$founders$phenotypes)<length(markers)){
		cat("WARNING: There are only",nrow(population$founders$phenotypes),"markers in population$founders$phenotypes, function will plot them all.\n")
		markers <- 1:nrow(population$founders$phenotypes)
	}
	for(i in markers) if(!(i %in% 1:nrow(population$offspring$phenotypes))){
		stop("ERROR: There is no marker number: ",i,", stopping \n")
	}
	for(i in markers) if(!(i %in% 1:nrow(population$founders$phenotypes))){
		stop("ERROR: There is no marker number: ",i,", stopping \n")
	}
	children <- population$offspring$phenotypes[markers,]
	founders <- population$founders$phenotypes[which(rownames(population$founders$phenotypes) %in% rownames(children)),]
	rownames(children) <- markers
	boxplot(t(children), ylim=c(min(children), max(children)),	xlab="Marker", ylab="Expression value", main="Children gene expression data")
	points(apply(founders,1,mean),col="green", pch=95, cex=3)
	points(apply(founders,1,max),col="red", pch=24, cex=1)
	points(apply(founders,1,min),col="blue", pch=25, cex=1)
}

############################################################################################################
#									*** plotMarkerDistribution ***
#
# DESCRIPTION:
#	plotting histogram of distribution of values for single marker and specified number of normal distribution 
#	curves, fitted to data using EM algorithm
# 
# PARAMETERS:
# 	population - population type object, must contain founders phenotypic data.
# 	marker - number or name of marker to be printed 
# 	nrDistributions - numbers of normal distributions to be fitted
# 	logarithmic - TRUE - log(data) will be used instead of raw data
#
# OUTPUT:
#	plot
#
############################################################################################################
plotMarkerDistribution <- function(population,marker,nrDistributions,logarithmic=FALSE){
	if(missing(population)) stop("population object not provided.\n")
	if(missing(marker)) stop("marker not specified.\n")
	if(missing(nrDistributions)) stop("nrDistributions not specified.\n")
	if(length(marker)!=1) stop("plotMarkerDistribution can plot only one marker art once.\n")
	check.population(population)
	phenotypeRow <- population$offspring$phenotypes[marker,]
	if(logarithmic) phenotypeRow <- log(phenotypeRow)
	EM<-normalmixEM(phenotypeRow, k=nrDistributions)
	if(logarithmic){
		xlab <- "log(expression values)"
	}else{
		xlab <- "expression values"
	}
	h <- NULL
	len <- vector(mode="numeric",length=nrDistributions)
	for(i in 1:nrDistributions){
		len[i]<-length(phenotypeRow)*EM$lambda[which(EM$mu==sort(EM$mu)[i])]
		startVal <- sum(len[1:i-1])
		x <- sort(phenotypeRow)[startVal:(startVal+len[i])]
		h[i] <- hist(x,breaks=50)
	}
 	h0<-hist(phenotypeRow,breaks=50,col="grey62",border="grey70",xlab=xlab,ylab="Number of counts",main="Distribution of expression values for selected marker")
	colorP <- vector(mode="character",length=3)
	colorP[1] <- rgb(1,0,0)
	colorP[2] <- rgb(0,1,0)
	colorP[3] <- rgb(0,0,1)
	for(i in 1:nrDistributions){
		abline(v=EM$mu[i],col=colorP[i%%3+1])
		abline(v=c(EM$mu[i]-EM$sigma[i],EM$mu[i]+EM$sigma[i]),col=colorP[i%%3+1],lty=2)
		startVal <- sum(len[1:i-1])
		x <- sort(phenotypeRow)[startVal:(startVal+len[i])]
		xfit<-seq(min(x),max(x),length=40) 
		yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
		yfit <- yfit*diff(h[[i]][1:2])*length(x) 
		lines(xfit, yfit, col=colorP[i%%3+1], lwd=2)
	}
}

############################################################################################################
#									*** projectOldMarkers ***
#
# DESCRIPTION:
#	function calculating matrix of corelations between chromosomes
# 
# PARAMETERS:
# 	cross - object of R/qtl cross type
# 	cur_map - map to be used in comparison
# 
# OUTPUT:
#	matrix of corelations
#
############################################################################################################
projectOldMarkers <- function(cross,population,map=c("genetic","physical"),label=c("positions","names","no"),...){
	if(missing(cross)) stop("No cross object provided!\n")
	if(missing(population)) stop("No population object provided!\n")
  map <- match.arg(map)
	if(map=="genetic"){
    cur_map <- population$maps$genetic
  }else{
    cur_map <- population$maps$physical
  }
  if(is.null(cur_map)) stop("no ",map," map provided!")  
  label <- match.arg(label)
	qp_ <- NULL
  cross <- jittermap(cross)
	for(i in 1:nchr(cross)){
		qp_ <- c(qp_,cross$geno[[i]]$map)
	}
  presentMarkers <- as.character(rownames(cur_map)[which(rownames(cur_map) %in% names(qp_))])
  qc <- cur_map[presentMarkers,1]
	qp <- qp_[presentMarkers]
	if(label=="positions"){
		qn <- round(cur_map[presentMarkers,2])
	}else if(label=="names"){
    qn <- presentMarkers
	}else if(label=="no"){
		qn <- rep('',length(presentMarkers))
	}
	cross <- sim.geno(cross)
	qtl <- makeqtl(cross,qc,qp,qn)
	plot(qtl,...)
}

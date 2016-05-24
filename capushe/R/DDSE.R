setClass(
	Class="DDSE",
        representation=representation(
	model="character",
	kappa="numeric",
	ModelHat="list",
	interval="list",
	graph="list")
)

setMethod("print","DDSE",
	function(x,...){
		print(x@model)
	}
)

setMethod("show","DDSE",
	function(object){
		print(object@model)
	}
)

setMethod("summary","DDSE",
	function(object,checking=FALSE,checkgraph=FALSE){
		cat("\nSelected model: ",as.character(object@model),"\n",sep="")
		cat("\nCorresponding Slope interval: [",as.character(signif(object@interval$interval[1],3)),";",as.character(signif(object@interval$interval[2],3)),"](",as.character(signif(object@interval$percent_of_points,3)*100),"% of points)\n",sep="")
		cat("\nSelected model complexity: ",as.character(object@graph$Complexity),"\n",sep="")
	if (checking){
		result=list(object@model,object@kappa,object@ModelHat,object@interval)
		names(result)=c("model","kappa","ModelHat","interval")
		if (checkgraph==TRUE){
			Tempo=list(object@graph)
			names(Tempo)=c("graph")
			result=c(result,Tempo)
			}
		result
		}
	}
)


setMethod(
	f="plot",
	signature="DDSE",
	definition=function(x,y,newwindow=TRUE,...){
	  
	  if(!is.logical(newwindow)){
	    stop("newwindow must be logical")
	  }
	  refgraph=list(ask=par()$ask,
	                mgp=par()$mgp,
	                oma=par()$oma,
	                xaxs=par()$xaxs,
	                mfrow=par()$mfrow,
	                cex.axis=par()$cex.axis)


	plength=length(x@graph$model)-1
	p=plength-x@interval$point_using+2

	Model=x@ModelHat$model_hat[x@ModelHat$point_breaking[x@ModelHat$imax]]

	if (newwindow){	dev.new(width=14)}
	par(mgp = c(3, 0.5, 0))
	layout(matrix(c(1,1,1,1,1,1,1,1,2,3,2,3,2,3),ncol=7))
	Couleur=c(rep("blue",(p-1)),rep("red",(plength-p+2)))
	plot(x=x@graph$pen,y=-x@graph$contrast,main="Contrast representation",xlab=expression(paste("\n\n",pen[shape](m)," (labels : Model names)")),ylab=expression(-gamma[n] (hat(s)[m])),col=Couleur,xaxt = "n")
	legend("bottomright",legend=c(paste("The regression line is computed with",as.character(plength-p+2),"points")),lty=1,col=c("red"))

	labelpen=x@graph$pen
	labelmodel=x@graph$model
	pas=(max(labelpen)-min(labelpen))/60*11.5/par("din")[1]
	leng=length(labelpen)
	Coord=c()
	i=1
	j=2
	while (j<leng){
		while ((j<leng)*((labelpen[j]-labelpen[i])<pas)){
			Coord=c(Coord,-j)
			j=j+1
			}
		i=j
		j=j+1
		}

	if (length(Coord)>0){
		labelpen=labelpen[Coord]
		labelmodel=labelmodel[Coord]
		}
	axis(1,labelpen,labelmodel,las=2)
	
	Cof=x@graph$reg$coefficients
	abline(a=Cof[1],b=Cof[2],col="red")
	abscisse=seq(plength+1,2)
	plot(x@kappa,main="Successive slope values",xlab=expression(paste("Number of points (",pen[shape](m),",",-gamma[n] (hat(s)[m]),") for the regression")),ylab=expression(kappa),type="l",xaxt = "n")

	abscisse2=seq(1,plength+1,length=floor(13*11.5/par("din")[1]))
	axis(1,abscisse2,abscisse[abscisse2])

	points(x@kappa,col="red",pch=18)
	plot(x@ModelHat$model_hat,main="Selected models with respect to the successive slope values",xlab=expression(paste("Number of points (",pen[shape](m),",",-gamma[n] (hat(s)[n]),") for the regression")),ylab="Model",pch=20,yaxt="n",xaxt = "n")
	axis(1,abscisse2,abscisse[abscisse2])

	labelhat=x@ModelHat$model_hat
	pas=(max(labelhat)-min(labelhat))/61*5.75/par("din")[2]
	leng=length(labelhat)
	Coord=c()
	i=1
	j=2
	while (j<leng){
		while ((j<leng)*((labelhat[j]-labelhat[i])<pas)){
			Coord=c(Coord,-j)
			j=j+1		
			}
		i=j
		j=j+1
		}
	Coord=c(Coord,-which(duplicated(labelhat)==TRUE))
	if (length(Coord)>0){
		labelhat=labelhat[Coord]
		}

	axis(2,labelhat,x@graph$model[labelhat],las=1)
	lines(x=x@ModelHat$point_breaking[x@ModelHat$imax]:(x@ModelHat$point_breaking[x@ModelHat$imax]+x@ModelHat$number_plateau[x@ModelHat$imax]-1),y=rep(Model,x@ModelHat$number_plateau[x@ModelHat$imax]),col="blue")
	par(refgraph)
	}
)


DDSE = function(data,pct=0.15,point=0,psi.rlm=psi.bisquare,scoef=2){

if(!(is.numeric(pct))){
	stop("pct must be numeric")
	}
if(!(is.numeric(point))){
	stop("point must be numeric")
	}
if (!(floor(point)==point)|(point<0)){
	stop("point must be an positive integer")
	}
if(!(is.numeric(scoef))){
	stop("scoef must be numeric")
	}

if (is.character(data)){data=read.table(data)}



if (length(data[1,])!=4){
	stop("data must have 4 columns with name, penshape, complexity and contrast")
	}
data=data.frame(data)
names(data)=c("model","pen","complexity","contrast")

if(any(is.na(data))){
	bad=which(is.na(data),arr.ind=TRUE)[,1]
	repbad=which(duplicated(bad),arr.ind=TRUE)
	if (length(repbad)!=0){bad=bad[-repbad]}
	data=data[-bad,]
	if (length(bad)==1){
		warning(paste("1 line has been removed by DDSE"))
		}
		else{warning(paste(as.character(length(bad)),"lines have been removed by DDSE"))}
	}
mlength=length(data$model)
if (mlength<10){
	stop("At least 10 observations are needed")
	}

if (!is.numeric(data[,2])){
	stop("pen must be numeric")
	}
if (!is.numeric(data[,3])){
	stop("complexity must be numeric")
	}
if (!is.numeric(data[,4])){
	stop("contrast must be numeric")
	}

if ((length(data$pen)!=mlength)|(length(data$complexity)!=mlength)|(length(data$contrast)!=mlength)){
	stop("The lengths of the columns are not equal")
	}

if ((pct<0)|(pct>1)){
	stop("pct must be between 0 and 1")
	}
if (point>mlength){
	stop("point must be smaller than the number of models")
	}
if (!prod(data$complexity>=0)){
	stop("Complexity must be positive")
	}

data=data[order(data$pen,data$contrast),]
plength=length(data$pen)
for (i in plength:2){
	if (data$pen[i]==data$pen[i-1]){
		data=data[-i,]
		mlength=mlength-1
		}
	}
Tempoligne=data$pen[order(data$complexity)]
if (!prod((Tempoligne[2:mlength]-Tempoligne[1:(mlength-1)])>0)){
	stop("Penshape should be an increasing function of the complexity")
	}

P=data$pen
plength=length(P)-1
kappa=numeric(plength)
couplepen=P
couplegamma=-data$contrast
if (!is.function(psi.rlm)){
	warning("The function lm is used instead of rlm")
	for (p in 1:(plength)){
		kappa[p]=lm(couplegamma~couplepen)$coefficients[2]    
		couplepen=couplepen[-1]
		couplegamma=couplegamma[-1]
		}
	}
else{
	options("warn"=-1)
	for (p in 1:(plength)){
		kappa[p]=rlm(couplegamma~couplepen,psi=psi.rlm)$coefficients[2]    
		couplepen=couplepen[-1]
		couplegamma=couplegamma[-1]
		}
	options("warn"=0)
	}
if (!prod(kappa>0)){
	warning("Some elements in Kappa are negative")
	}

mhat=numeric(plength)
for (p in 1:plength){
	mhat[p]=which.min(data$contrast+scoef*kappa[p]*data$pen)
	}

Pi=c(1)
N=c()
number=1
for (p in 2:plength){
	if (mhat[p]!=mhat[p-1]){
		Pi=c(Pi,p)
		N=c(N,number)
		number=0
		}
	number=number+1
	}
N=c(N,number)
Pi=c(Pi,p)

Imax=length(N)
Ntot=sum(N)
if (point==0){
	const=pct*Ntot
	while (N[Imax]<const){
		Imax=Imax-1
		if (Imax==0){
			stop("pct is too high")
			}
		}
	}
else {
	while (N[Imax]<point){
		Imax=Imax-1
			if (Imax==0){
			stop("point is too high")
			}
		}
	}

No=N[Imax]
p1=Pi[Imax]
p=p1+floor(No/2)
Model=mhat[p]


result=list(data$model[Model])
names(result)=c("Model")

ModelHat=list(mhat,Pi,N,Imax)
names(ModelHat)=c("model_hat","point_breaking","number_plateau","imax")

Interval=kappa[Pi[Imax]:(Pi[Imax+1]-1)]
inter=c(min(Interval),max(Interval))
names(inter)=c("min","max")
pourcent=(No+1)/(Ntot+1)
interval = list(plength-p+2,inter,pourcent)
names(interval)=c("point_using","interval","percent_of_points")

if (!is.function(psi.rlm)){
	reg=lm(-data$contrast[p:(plength+1)]~data$pen[p:(plength+1)])
	}
else{
	options(warn=-1)
	reg=rlm(-data$contrast[p:(plength+1)]~data$pen[p:(plength+1)],psi=psi.rlm)
	options(warn=0)
}
graph=list(reg,data$pen,data$contrast,data$model,data$complexity[mhat[Pi[Imax]]])
names(graph)=c("reg","pen","contrast","model","Complexity")


new(Class="DDSE",model=as.character(data$model[Model]),kappa=kappa,ModelHat=ModelHat,interval=interval,graph=graph)

}


setMethod(
	f="validation",
	signature="DDSE",
	definition=function(x,data2,newwindow=TRUE,...){
	  
	  if(!is.logical(newwindow)){
	    stop("newwindow must be logical")
	  }
	  refgraph=list(ask=par()$ask,
	                mgp=par()$mgp,
	                oma=par()$oma,
	                xaxs=par()$xaxs,
	                mfrow=par()$mfrow,
	                cex.axis=par()$cex.axis)

	if (is.character(data2)){data2=read.table(data2)}
	if (length(data2[1,])!=4){
		stop("data2 must have 4 columns with name, penshape, complexity and contrast")
		}
	if(any(is.na(data2))){
		bad=which(is.na(data2),arr.ind=TRUE)[,1]
		repbad=which(duplicated(bad),arr.ind=TRUE)
		if (length(repbad)!=0){bad=bad[-repbad]}
		data2=data2[-bad,]
		if (length(bad)==1){
			warning(paste("1 line in data2 has been removed by validation"))
			}
			else{warning(paste(as.character(length(bad)),"lines in data2 have been removed by validation"))}
		}
	data2=data.frame(data2)
	names(data2)=c("model","pen","complexity","contrast")

	plength2=length(data2[,1])
	plength=length(x@graph$model)-1
	p=plength-x@interval$point_using+2
	pen=c(x@graph$pen,data2$pen)
	contrast=c(x@graph$contrast,data2$contrast)
	model=c(x@graph$model,data2$model)

	if (newwindow){	dev.new(width=14) }
	par(mgp = c(3, 0.5, 0))
	Couleur=c(rep("blue",(p-1)),rep("red",(plength-p+2)))
	Cof=x@graph$reg$coefficients
	plot(xlim=c(min(pen),max(pen)),ylim=c(min(-contrast),max(-contrast,Cof[1]+Cof[2]*data2$pen)),x=x@graph$pen,y=-x@graph$contrast,main="Contrast representation",xlab=expression(paste(pen[shape](m)," (labels : Model names)")),ylab=expression(-gamma[n] (hat(s)[m])),col=Couleur,xaxt = "n")
	points(x=data2$pen,y=-data2$contrast,col="black",pch="*")

	legend("bottomright",legend=c(paste("The regression line is computed with",as.character(plength-p+2),"points"),"Validation points"),lty=c(1,0),col=c("red","black"),pch=c(".","*"))

	labelpen=pen
	labelmodel=model
	pas=(max(labelpen)-min(labelpen))/100*11.5/par("din")[1]
	leng=length(labelpen)
	Coord=c()
	i=1
	j=2
	while (j<leng){
		while ((j<leng)*((labelpen[j]-labelpen[i])<pas)){
			Coord=c(Coord,-j)
			j=j+1
			}
		i=j
		j=j+1
		}

	if (length(Coord)>0){
		labelpen=labelpen[Coord]
		labelmodel=labelmodel[Coord]
		}
	axis(1,labelpen,labelmodel,las=2)

	abline(a=Cof[1],b=Cof[2],col="red")
	abscisse=seq(plength+1,2)
	par(refgraph)
	}
)
#######################dependences###########################

# requrired by the plot-functions
library(ggplot2)
library(MASS)
# required by independence test bootstrap
library(foreach)
library(matrixStats)
######################localgauss###################################
   
  	# x,y:               Datavectors
    # b1,b2;             Bandwidth in x and y direction    	
  	# gsize              Gridsize (used only if xy.mat is NULL)
    # hthresh            Leave out points where the kernel density estimate is lower than this.
    #                    (used only if xy.mat is NULL)
    # xy.mat             Points where the local Gaussian correlation is computed.
    #                    If xy.mat is not NULL, this option overrides the gsize and hthresh options




localgauss <- function(x,y,b1=1,b2=1,gsize=15,hthresh=0.001,xy.mat=NULL){

# work out if default or custom points (i.e. xy.mat is given) are to be used in
# esitmation
if(is.null(xy.mat)){
  #Create xy.mat by selecting points with density>hthresh	
	xy.est=kde2d(x, y, n=gsize, h=c(bandwidth.nrd(x), bandwidth.nrd(y)), lims=c(min(x),max(x),min(y),max(y)))  #density estimation
	xy.est=con2tr(xy.est)                                                                                      #convert list to dataframe
	xy.est=subset(xy.est,xy.est[,3]>hthresh)                                                                   #remove points with to low density
	xy.mat.internal=as.matrix(xy.est[,1:2])     
} else {
  # xy.mat is provided
  xy.mat.internal=xy.mat
}
# handle data in format suitable for the fortran function
n = length(x)
if(is.matrix(xy.mat.internal)){
  nrep = dim(xy.mat.internal)[1];
  x0 = xy.mat.internal[,1];
  y0 = xy.mat.internal[,2];
}
else{
  nrep = 1;
  x0 = as.vector(xy.mat.internal[1]);
  y0 = as.vector(xy.mat.internal[2]);
}
# do not alter function arguments
hx = abs(b1);
hy = abs(b2);
# the fortran function takes separate bandwidths for each row in
# xy.mat
if(length(hx)==1){
  hxx = rep(hx,nrep);
}
else{
  hxx = as.vector(hx);
}

if(length(hy)==1){
  hyy = rep(hy,nrep);
}
else{
  hyy = as.vector(hy);
}
# the fortran function takes an additional argument - rhof -
# that is not used by this package. Essentially, the fortran function
# has the ability to estimate means and std.devs. while keeping correlation
# fixed at some value rhof. Passing rhof that is not a valid correlation causes
# rho to be estimated. However we hide this functionality as it has not
# been tested thoroughly and have not proven useful in practice, and  
# thus a vector of -2.0s is passed to cause rho always to be estimated.

rhofix = rep(-2.0,nrep);

# allocate space for output
res.in = matrix(0.0,nrow=nrep,ncol=7)
hess.in = array(0.0,dim=c(nrep,5,5))
# actual call to the fortran function
f.out <- .Fortran("localgauss",
                  as.integer(n),
                  as.double(as.vector(x)),
                  as.double(as.vector(y)),
                  as.integer(nrep),
                  as.double(as.vector(x0)),
                  as.double(as.vector(y0)),
                  as.double(as.vector(hxx)),
                  as.double(as.vector(hyy)),
		  as.double(as.vector(rhofix)),
                  as.double(res.in),
		  as.double(hess.in),PACKAGE="localgauss")


# format output for easy return
res.out = as.matrix(matrix((f.out[10])[[1]],nrow=nrep,ncol=7));
hess.out = array(f.out[11][[1]],dim=c(nrep,5,5));
par.est= matrix(res.out[,1:5],nrow=nrep,ncol=5);
colnames(par.est) = c("mu_1","mu_2","sig_1","sig_2","rho");
eflag = res.out[,7]
ret = list(x=x,y=y,xy.mat=xy.mat.internal,b1=hx,b2=hy,par.est=par.est,eflag=eflag,hessian=hess.out)
attr(ret,"class") <- "localgauss"
return(ret);
}

######################plot.localgauss###################################

# lgobj:                "localgauss" object from function localgauss
# plot.text:            Should the numerical local correlation value be printed on each tile?
# plot.points:          Should the original observation points be overlain?
# tsize:                Font size used if plot.text is true.
# lowcol:               Colour used to indicate small values of correlation.                                 <---updated
# highcol:              Colour used to indicate large values of corelation.                                  <---updated
# point.col:            Color used for observations points (if plot.points is true).
# point.size:           Size of observations points (if plot.points is true).
# xlab:                 label of x-axis.
# ylab:                 label of y-axis.
# divergent.col.grad: 


plot.localgauss=function(x,..., plot.text=TRUE,plot.points=FALSE,tsize=3,
  lowcol="cyan",highcol="magenta",point.col="black",
  point.size=NULL,xlab="",ylab="",divergent.col.grad=T){
	# to provide consistency with S3 generic functions
	lgobj = x
	# Dummies to trick R CMD check 
	x.obs <- NULL; rm(x.obs);
	y.obs <- NULL; rm(y.obs);
	#Objects to plot.
	dat.obs = data.frame(x.obs=lgobj$x,y.obs=lgobj$y)
	rho = lgobj$par.est[,5]
	x = lgobj$xy.mat[,1]
	y = lgobj$xy.mat[,2]
	
	#Format rho with leading sign.
	rho.label = formatC(rho,digits=2, format="f",flag="+")                      
	#Dataframe for plotting rho.
	d = data.frame(x=x,y=y,rho=rho,rho.label=rho.label,stringsAsFactors=FALSE)  
	
	#Make a ggplot using the geom "tile".
	g = ggplot() + layer(data=d, mapping=aes(x=x, y=y, fill=rho, label=rho.label),geom="tile",stat="identity",position="identity",params=list(na.rm=TRUE)) # Plot estimated rho.
	
	# Add a new colour gradient to graph with possible different colours. Either divergent from 0, or continously from lowcol to highcol.
	lims=c(floor(min(rho)*10),ceiling(max(rho)*10))/10      
	breaks = seq(lims[1],lims[2],by=(lims[2]-lims[1])/10)  #adjust breaks after the range of rho
	if(divergent.col.grad)
	gr = scale_fill_gradient2(midpoint=0,low=lowcol, high=highcol, space="Lab", limits=lims,breaks=breaks)	
	else
	gr=scale_fill_gradient(low = lowcol, high = highcol, space = "Lab",limits=lims, breaks=breaks)	
	
    g = g + gr 
	
	#Should we overlay the original observations?
	if(plot.points) 
		if(is.null(point.size))
			g = g + layer(data=dat.obs, mapping=aes(x=x.obs,y=y.obs), geom="point",
                stat="identity",position="identity",params=list(na.rm = TRUE,colour=point.col))
		else
			g = g + layer(data=dat.obs, mapping=aes(x=x.obs,y=y.obs), geom="point", stat="identity",position="identity",params=list(colour=point.col,size=point.size,na.rm = TRUE))
	
	#Add numerical values of the local correlation to the graph if wanted
	if(plot.text)
		g = g + layer(stat = "identity",position = "identity",data=d, mapping=aes(x=x, y=y, label=rho.label),geom="text",
            params=list(size=tsize,na.rm = FALSE))
	
	#Add labels
	g = g + scale_x_continuous(name=xlab)
	g = g + scale_y_continuous(name=ylab)
		
	# Return plot
	return(g)
	
}

#localgauss.indtest     Pointwise Independence test based on LGC.

#locobj           "localgauss" object from function localgauss
# R:              bootstrap size
# alpha:          significance level (note: two sided test)
# seed:           random seed in bootstrap


localgauss.indtest=function(locobj,R=10,alpha=0.10,seed=1){
		
	#compute the sample correlation
	sample.rho=locobj$par.est[,5]	
	# Dummies to trick R CMD check 
	i <- NULL; rm(i);
	#Then the bootstrap begins:
	######################################################################################################	
	
	set.seed(seed)
	
	bootstrap.rho=foreach(i=seq(R),.inorder=FALSE,.combine="rbind",
	.packages='localgauss') %dopar%{                   #foreach supports a parallel backend
		
		# set seed for each iteration to ensure reproducibility
		set.seed((i-1)*R+seed)
		#Sample from the estimated null distribution, i.e. the product of the marginal empirical dfs:
		
		x.bootstrap=locobj$x[sample(seq(locobj$x),replace=TRUE)]                               #bootstrap sample of x nr.i
		y.bootstrap=locobj$y[sample(seq(locobj$x),replace=TRUE)]                               #bootstrap sample of y nr.i
		
		#local gaussian correlation estimate based on bootstrap sample nr.i
		localgauss(x=x.bootstrap,y=y.bootstrap,xy.mat=locobj$xy.mat,b1=locobj$b1,b2=locobj$b2)$par.est[,5]	
		
	} #End of bootstrapping
	
	#qcheck. Used to evaluate if sample LGC exceeds the quantiles of the bootstrap distribution   
	qcheck=function(a){
		if(a[1]<a[3]){q=-1}
		if(a[1]>a[2]){q=1}
		if((a[1]>=a[3])&(a[1]<=a[2])){q=0}
		return(q)
	} #end of qcheck
	
	#find lower and upper quantiles for each grid
	upper=colQuantiles(bootstrap.rho,probs=1-alpha/2)
	lower=colQuantiles(bootstrap.rho,probs=alpha/2)
	
	#apply qcheck to the sample LGC. Assigns +1 if signifcantly positive, -1 if significantly negative and 0 if neither.
	test.results=apply(cbind(sample.rho,upper,lower),1,qcheck)
	
	#return objects suitable for plotting along with auxiliary data
	ret=list("localgauss"=locobj,"upper"=upper,"lower"=lower,"test.results"=test.results)
	attr(ret,"class")<-"localgauss.indtest"
	return(ret)
	
}   #end of localgauss.indtest



#plot.localgauss.indtest             Plots the result of a local independence test

# indtestobj:   "localgauss.indtest" object from function localgauss.indtest
# plot.points: Should the original observation points be overlain?
# poscol:      Colour used to indicate significantly positive correlation
# negcol:      Colour used to indicate significantly negative correlation
# point.col:   Color used for observations points (if plot.points is true).
# point.size:  Size of observations points (if plot.points is true).
# xlab:        label of x-axis.
# ylab:        label of y-axis.


plot.localgauss.indtest=function(x,...,plot.points=F,
poscol="magenta",negcol="cyan",point.col="black",
point.size=NULL,xlab="",ylab=""){
	# for consistency with S3 generic functions
	indtestobj = x
	# Dummies to trick R CMD check 
	rho <- NULL; rm(rho); 
	x.obs <- NULL; rm(x.obs);
	y.obs <- NULL; rm(y.obs); 
	#First get the objects to plot.
	test = indtestobj$test.results
	x = indtestobj$localgauss$xy.mat[,1]
	y = indtestobj$localgauss$xy.mat[,2]

	dat.obs = data.frame(x.obs=indtestobj$localgauss$x,y.obs=indtestobj$localgauss$y)
	
	#Attach colour strings to the results of the test
	test[test==1]=poscol
	test[test==-1]=negcol
	test[test==0]="white"
	
	col=c(poscol,"white",negcol)
	
	#Put the test results with corresponding coordinates into a dataframe
	dat = data.frame(x=x,y=y,rho=test)
	
	#Plot test results.
	g = ggplot() + layer(data=dat, mapping=aes(x=x, y=y, fill=rho,labels=rho),geom="tile",
        stat="identity",position="identity",params=list(na.rm=TRUE))
	
	#Adjust label
	g=g+scale_fill_identity(labels = c("+","0","-"), name = "Test results",breaks=col,guide="legend")     
	
	#Should we overlay the original observations?
	if(plot.points) 
		if(is.null(point.size))
			g = g + layer(data=dat.obs, 
			mapping=aes(x=x.obs,y=y.obs),
			 geom="point", params=list(colour=point.col,na.rm=TRUE),
             stat="identity",position="identity")
		else
			g = g + layer(data=dat.obs, 
			mapping=aes(x=x.obs,y=y.obs),
			 geom="point", params=list(colour=point.col,size=point.size,na.rm=TRUE),
             stat="identity",position="identity")
	
    #Add x,y labels
	g = g + scale_x_continuous(name=xlab)
	g = g + scale_y_continuous(name=ylab)
	
	#Return plot
	return(g)
	
	} #end of plot.localgauss.indtest



dcircmix<-function(x,model=NULL,dist=NULL,param=NULL){
	x <- conversion.circular(x, units = "radians", zero = 0, rotation = "counter")
	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	if (is.null(model) & is.null(dist)) stop("No model specified")
	x.na <- is.na(x)
	x <- x[!x.na]
	n <- length(x)
	if (sum(x.na)>0) warning("Missing values were removed")
	if (n==0) stop("No observations (at least after removing missing values)")
	if (!is.null(dist)){
		nobj <- length(param)
		if (nobj<5){
			for (i in 1:length(param)){
				if (length(dist)!=length(param[[i]])) stop ("Length of the objects of the list 'param' must be equal to the length of argument 'dist'")
			}
		}else{
			for (i in 1:4){
				if (length(dist)!=length(param[[i]])) stop ("Length of the first four objects of the list 'param' must be equal to the length of argument 'dist'. Length of fifth object must be 1.")
			}
		}
		if (sum(param[[1]])!=1){
			warning ("Proportions must sum 1. Proportions were normalized by the sum")
			ndist<-length(dist)
			param[[1]]<-param[[1]]/sum(param[[1]])
		}
	}
	if (!is.null(model)){
		if (!is.numeric(model)) stop("argument 'model' must be numeric")
		if (!any(model==1:20)) stop("Value specified for argument 'model' is not valid")
		if (model==1){ 
			dist<-"unif"
			param <- list(p=1,mu=0,con=0)
		}else if (model==2){
			dist<-"vm"
			param <- list(p=1,mu=pi,con=1)
		}else if (model==3){
			dist<-"wn"
			param <- list(p=1,mu=pi,con=0.9)
		}else if (model==4){
			dist<-"car"
			param <- list(p=1,mu=pi,con=0.5)
		}else if (model==5){
			dist<-"wc"
			param <- list(p=1,mu=pi,con=0.8)
		}else if (model==6){
			dist<-"wsn"
			param <- list(p=1,mu=pi,con=1,sk=20)
		}else if (model==7){
			dist<-c("vm","vm")
			param <- list(p=c(1/2,1/2),mu=c(0,pi),con=c(4,4))
		}else if (model==8){
			dist<-c("vm","vm")
			param <- list(p=c(1/2,1/2),mu=c(2,4),con=c(5,5))
		}else if (model==9){
			dist<-c("vm","vm")
			param <- list(p=c(1/4,3/4),mu=c(0,pi/sqrt(3)),con=c(2,2))
		}else if (model==10){
			dist<-c("vm","wc")
			param <- list(p=c(4/5,1/5),mu=c(pi,4*pi/3),con=c(5,0.9))
		}else if (model==11){
			dist<-c("vm","vm","vm")
			param <- list(p=c(1/3,1/3,1/3),mu=c(pi/3,pi,5*pi/3),con=c(6,6,6))
		}else if (model==12){
			dist<-c("vm","vm","vm")
			param <- list(p=c(2/5,1/5,2/5),mu=c(pi/2,pi,3*pi/2),con=c(4,4,4))
		}else if (model==13){
			dist<-c("vm","vm","vm")
			param <- list(p=c(2/5,2/5,1/5),mu=c(0.5,3,5),con=c(6,6,24))
		}else if (model==14){
			dist<-c("vm","vm","vm","vm")
			param <- list(p=c(1/4,1/4,1/4,1/4),mu=c(0,pi/2,pi,3*pi/2),con=c(12,12,12,12))
		}else if (model==15){
			dist<-c("vm","wc","wn","wsn")
			param <- list(p=c(1/4,3/10,1/4,1/5),mu=c(pi+2,pi-1,pi+0.5,6),con=c(3,0.6,0.9,1),sk=c(0,0,0,1))
		}else if (model==16){
			dist<-c("vm","vm","vm","vm","vm")
			param <- list(p=c(1/5,1/5,1/5,1/5,1/5),mu=c(pi/5,3*pi/5,pi,7*pi/5,9*pi/5),con=c(18,18,18,18,18))
		}else if (model==17){
			dist<-c("car","wc")
			param <- list(p=c(2/3,1/3),mu=c(pi,pi),con=c(0.5,0.9))
		}else if (model==18){
			dist<-c("vm","vm","vm","vm")
			param <- list(p=c(1/2,1/6,1/6,1/6),mu=c(pi,pi-0.8,pi,pi+0.8),con=c(1,30,30,30))
		}else if (model==19){
			dist<-c("vm","vm","vm","vm","vm")
			param <- list(p=c(4/9,5/36,5/36,5/36,5/36),mu=c(2,4,3.5,4,4.5),con=c(3,3,50,50,50))
		}else if (model==20){
			dist<-c("wc","wc","wsn","wsn")
			param <- list(p=c(1/6,1/6,1/3,1/3),mu=c(3*pi/4,7*pi/4,0,pi),con=c(0.9,0.9,0.7,0.7),sk=c(0,0,20,20))
		}
	}
	mix<-0
	for (i in 1:length(dist)){
		distribution <- dist[i] 
		if (distribution=="unif"){
			mix <- mix + param[[1]][i]*dcircularuniform(x)
		}else if (distribution=="vm"){
			mix <- mix + param[[1]][i]*dvonmises(x,mu=circular(param[[2]][i]),kappa=param[[3]][i])
		}else if (distribution=="car"){
			mix <- mix + param[[1]][i]*dcardioid(x,mu=circular(param[[2]][i]),rho=param[[3]][i])
		}else if (distribution=="wc"){
			mix <- mix + param[[1]][i]*dwrappedcauchy(x,mu=circular(param[[2]][i]),rho=param[[3]][i])
		}else if (distribution=="wn"){
			mix <- mix + param[[1]][i]*dwrappednormal(x,mu=circular(param[[2]][i]),rho=param[[3]][i])
		}else if (distribution=="wsn"){
			mix <- mix + param[[1]][i]*dwsn(x,xi=circular(param[[2]][i]),eta=param[[3]][i],lambda=param[[4]][i])
		}
	}
	return(mix)
}



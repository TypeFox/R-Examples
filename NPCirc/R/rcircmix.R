rcircmix<-function(n,model=NULL,dist=NULL,param=NULL){
	if (!is.numeric(n)) stop("argument 'n' must be numeric")
	if (is.null(model) & is.null(dist)) stop("No model specified")
	if (!is.null(dist)){
		if (length(param)<3 | length(param)>4) stop ("Length of argument 'param' must be 3 or 4")
		for (i in 1:length(param)){
			if (length(dist)!=length(param[[i]])) stop ("Length of the objects of the list 'param' must be equal to the length of argument 'dist'")
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
			param <- list(p=c(2/3,1/3),mu=c(pi,pi),con=c(0.5,0.9),sk=c(0,0))
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

	p <- param[[1]]
    	ind <- as.numeric(cut(runif(n), c(0, cumsum(p)), include.lowest=TRUE))
    	pos <- split(seq_len(n), ind)
	nms <- names(pos)
	result <- rep(NA, n)
	for (i in seq_along(pos)){
		j <- as.numeric(nms[i])
		distribution <- dist[j]
		position <- pos[[i]]
		npos<-length(position)
		if (distribution=="unif"){
			 result[position] <- rcircularuniform(npos)
		}else if (distribution=="vm"){
			result[position] <- rvonmises(npos,mu=circular(param[[2]][j]),kappa=param[[3]][j])
		}else if (distribution=="car"){
			result[position] <- rcardioid(npos, mu=circular(param[[2]][j]), rho=param[[3]][j])
		}else if (distribution=="wc"){
			result[position] <- rwrappedcauchy(npos,mu=circular(param[[2]][j]),rho=param[[3]][j])
		}else if (distribution=="wn"){
			result[position] <- rwrappednormal(npos,mu=circular(param[[2]][j]),rho=param[[3]][j])
		}else if (distribution=="wsn"){
			result[position] <- rwsn(npos,xi=circular(param[[2]][j]),eta=param[[3]][j],lambda=param[[4]][j])
		}
	}
	return(circular(result))
}

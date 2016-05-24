##' @export
smcure <-
function(formula,cureform,offset=NULL,data,na.action=na.omit,model= c("aft", "ph"),link="logit", Var=TRUE,emmax=50,eps=1e-7,nboot=100)
{
	call <- match.call()
	model <- match.arg(model)
	cat("Program is running..be patient...")
	## prepare data
	data <- na.action(data)
	n <- dim(data)[1]
	mf <- model.frame(formula,data)
	cvars <- all.vars(cureform)
	Z <- as.matrix(cbind(rep(1,n),data[,cvars]))
	colnames(Z) <- c("(Intercept)",cvars)
	if(!is.null(offset)) {
	offsetvar <- all.vars(offset)
	offsetvar<-data[,offsetvar]}
	else offsetvar <- NULL
	Y <- model.extract(mf,"response")
	X <- model.matrix(attr(mf,"terms"), mf)
	if (!inherits(Y, "Surv")) stop("Response must be a survival object")
   	Time <- Y[,1]
   	Status <- Y[,2]
	bnm <- colnames(Z)
	nb <- ncol(Z)
	if(model == "ph") { 
			betanm <- colnames(X)[-1]
			nbeta <- ncol(X)-1}
	if(model == "aft"){ 
			betanm <- colnames(X)
			nbeta <- ncol(X)}
	## initial value 
	w <- Status
	b <- eval(parse(text = paste("glm", "(", "w~Z[,-1]",",family = quasibinomial(link='", link, "'",")",")",sep = "")))$coef
	if(model=="ph") beta <- coxph(Surv(Time, Status)~X[,-1]+offset(log(w)), subset=w!=0, method="breslow")$coef
	if(model=="aft") beta <- survreg(Surv(Time,Status)~X[,-1])$coef
     	## do EM algo
	emfit <- em(Time,Status,X,Z,offsetvar,b,beta,model,link,emmax,eps) 
		b <- emfit$b
		beta <- emfit$latencyfit
		s <- emfit$Survival
		logistfit <- emfit$logistfit
		  if(Var){
		   if(model=="ph") {b_boot<-matrix(rep(0,nboot*nb), nrow=nboot)
					beta_boot<-matrix(rep(0,nboot*(nbeta)), nrow=nboot)
					iter <- matrix(rep(0,nboot),ncol=1)}
  
		   if(model=="aft") {b_boot<-matrix(rep(0,nboot*nb), nrow=nboot)
					beta_boot<-matrix(rep(0,nboot*(nbeta)), nrow=nboot)}
		  tempdata <- cbind(Time,Status,X,Z)
		  data1<-subset(tempdata,Status==1);data0<-subset(tempdata,Status==0)
		  n1<-nrow(data1);n0<-nrow(data0)  
		  i<-1
		while (i<=nboot){
             id1<-sample(1:n1,n1,replace=TRUE);id0<-sample(1:n0,n0,replace=TRUE)
			bootdata<-rbind(data1[id1,],data0[id0,])
			bootZ <- bootdata[,bnm]
			if(model=="ph") bootX <- as.matrix(cbind(rep(1,n),bootdata[,betanm]))
			if(model=="aft") bootX <- bootdata[,betanm]
			bootfit <- em(bootdata[,1],bootdata[,2],bootX,bootZ,offsetvar,b,beta,model,link,emmax,eps)
			b_boot[i,] <- bootfit$b
		   	beta_boot[i,] <- bootfit$latencyfit
   			if (bootfit$tau<eps) i<-i+1}
		b_var <- apply(b_boot, 2, var)
		beta_var <- apply(beta_boot, 2, var)
		b_sd <- sqrt(b_var)
		beta_sd <- sqrt(beta_var)
		}
	fit<-list()
	class(fit) <- c("smcure")
	fit$logistfit <- logistfit	
	fit$b <- b
	fit$beta <- beta
	if(Var){
	fit$b_var <- b_var
	fit$b_sd <- b_sd
	fit$b_zvalue <- fit$b/b_sd
	fit$b_pvalue <- (1-pnorm(abs(fit$b_zvalue)))*2
	fit$beta_var <- beta_var
	fit$beta_sd <- beta_sd
	fit$beta_zvalue <- fit$beta/beta_sd
	fit$beta_pvalue <- (1-pnorm(abs(fit$beta_zvalue)))*2	}
	cat(" done.\n")
	fit$call <- call
	fit$bnm <- bnm
	fit$betanm <- betanm
	fit$s <- s
	fit$Time <- Time
	if(model=="aft"){
	error <- drop(log(Time)-beta%*%t(X))
	fit$error <- error}
	fit
	printsmcure(fit,Var)
	}


benchmark.pls<-function(X,y,m=ncol(X),R=20,ratio=0.8,verbose=TRUE,k=10,ratio.samples=1,use.kernel=FALSE,criterion="bic",true.coefficients=NULL){
    n<-floor(nrow(X)*ratio.samples)
    m.crash.krylov<-m.crash.lanczos<-vector(length=R)
    m.cv<-m.krylov<-m.naive<-m.lanczos<-vector(length=R)
    ntrain<-floor(n*ratio)
    time.cv<-time.krylov<-time.naive<-time.lanczos<-vector(length=R)
    DoF.cv<-DoF.krylov<-DoF.naive<-DoF.lanczos<-vector(length=R)
    mse.cv<-mse.naive<-mse.lanczos<-mse.krylov<-mse.null<-vector(length=R)
    DoF.complete<-matrix(,R,m+1)
    sigmahat.krylov<-sigmahat.naive<-sigmahat.lanczos<-sigmahat.cv<-sigmahat.null<-vector(length=R)
    mse.null<-vector(length=R)
	if (is.null(true.coefficients)==FALSE){
	model.cv<-model.naive<-model.lanczos<-model.krylov<-model.null<-vector(length=R)
	}
    for (i in 1:R){
        if (verbose==TRUE){
            cat(paste("iteration no ",i," \n"))
        }
        ###################################
        # split of training and test data #
        ###################################
        samples<-sample(nrow(X),n,replace=FALSE)
        XX<-X[samples,]
        yy<-y[samples]
        train<-sample(n,ntrain,replace=FALSE)
        Xtrain<-XX[train,,drop=FALSE]
        Xtest<-XX[-train,,drop=FALSE]
        ytrain<-yy[train]
        ytest<-yy[-train]
        ######################
        # Fitting the models #
        ######################
        #
        # null model
        mse.null[i]=mean((mean(ytrain)-ytest)^2)
	res<-ytrain -mean(ytrain)
	sigmahat.null[i]<-sqrt(sum(res^2)/(length(ytrain)-1))
        compute.jacobian=FALSE
        time.krylov[i]<-system.time(krylov.object<-pls.ic(Xtrain,ytrain,m=m,naive=FALSE,criterion=criterion,use.kernel=use.kernel,compute.jacobian=compute.jacobian))[3]
        sigmahat.krylov[i]<-krylov.object$sigmahat[krylov.object$m.opt+1]
        m.crash.krylov[i]<-krylov.object$m.crash
        m.krylov[i]<-krylov.object$m.opt
        # lanczos
        compute.jacobian=TRUE
        time.lanczos[i]<-system.time(lanczos.object<-pls.ic(Xtrain,ytrain,m=m,naive=FALSE,criterion=criterion,use.kernel=use.kernel,compute.jacobian=compute.jacobian))[3]
	sigmahat.lanczos[i]<-as.vector(lanczos.object$sigmahat[lanczos.object$m.opt+1])
        #sigmahat.lanczos[i,]<-as.vector(lanczos.object$sigmahat)
        DoF.complete[i,]=as.vector(lanczos.object$DoF)
        m.crash.lanczos[i]<-lanczos.object$m.crash
        m.lanczos[i]<-lanczos.object$m.opt
        # naive
        time.naive[i]<-system.time(naive.object<-pls.ic(Xtrain,ytrain,m=m,naive=TRUE,criterion=criterion,use.kernel=use.kernel,compute.jacobian=FALSE))[3]
	 sigmahat.naive[i]<-as.vector(naive.object$sigmahat[naive.object$m.opt+1])
      
        m.naive[i]<-naive.object$m.opt
        # cross-validation
        time.cv[i]<-system.time(cv.object<-pls.cv(Xtrain,ytrain,use.kernel=use.kernel,m=m))[3]
        m.cv[i]<-cv.object$m.opt
	res<-ytrain - rep(cv.object$intercept,length(ytrain)) - Xtrain%*%cv.object$coefficients
	sigmahat.cv[i]<-sqrt(sum(res^2)/(length(ytrain)-m.cv[i]-1))
        ######################
        # compute test error #  
        ######################
        m.max<-max(m.naive[i],m.krylov[i],m.lanczos[i],m.cv[i],1)
        pls.object=pls.model(Xtrain,ytrain,Xtest=Xtest,ytest=ytest,m=m.max,compute.DoF=FALSE,compute.jacobian=FALSE,use.kernel=use.kernel)
        mse.naive[i]<-pls.object$mse[m.naive[i]+1]
        mse.cv[i]<-pls.object$mse[m.cv[i]+1]
        mse.krylov[i]<-pls.object$mse[m.krylov[i]+1]
        mse.lanczos[i]<-pls.object$mse[m.lanczos[i]+1]
        ##############################
        # compute Degrees of Freedom #
        ##############################
        DoF.naive[i]=DoF.complete[i,m.naive[i]+1]
        DoF.krylov[i]=DoF.complete[i,m.krylov[i]+1]
        DoF.lanczos[i]=DoF.complete[i,m.lanczos[i]+1]
        DoF.cv[i]<-DoF.complete[i,m.cv[i]+1]
	#########################################################################
	# compute model error, if the true regression coefficients are provided #
	#########################################################################
	if (is.null(true.coefficients)==FALSE){
	model.cv[i]<-sum((true.coefficients-pls.object$coefficients[,m.cv[i]+1])^ 2)
	model.naive[i]<-sum((true.coefficients-pls.object$coefficients[,m.naive[i]+1])^ 2)
	model.null[i]<-sum((true.coefficients-pls.object$coefficients[,1])^ 2)
	model.krylov[i]<-sum((true.coefficients-pls.object$coefficients[,m.krylov[i]+1])^ 2)
	model.lanczos[i]<-sum((true.coefficients-pls.object$coefficients[,m.lanczos[i]+1])^ 2)
	
	}
        }
        #################
        # store results #
        #################
        namen<-c("CV","KRYLOV","LANCZOS","NAIVE","NULL")
        # mse
        MSE<-matrix(,R,5)
        MSE[,1]<-mse.cv
        MSE[,2]<-mse.krylov
        MSE[,3]<-mse.lanczos
        MSE[,4]<-mse.naive
        MSE[,5]<-mse.null
        MSE<-data.frame(MSE)
        colnames(MSE)=namen
        # number of components
        M<-matrix(,R,5)
        M[,1]<-m.cv
        M[,2]<-m.krylov
        M[,3]<-m.lanczos
        M[,4]<-m.naive
        M[,5]<-rep(0,R)
        M<-data.frame(M)
        colnames(M)=namen
        # Degrees of Freedom
        DoF<-matrix(,R,5)
        DoF[,1]<-DoF.cv
        DoF[,2]<-DoF.krylov
        DoF[,3]<-DoF.lanczos
        DoF[,4]<-DoF.naive
        DoF[,5]<-rep(1,R)
        DoF<-data.frame(DoF)
        colnames(DoF)=namen
        # time
        TIME<-matrix(,R,4)
        TIME[,1]<-time.cv
        TIME[,2]<-time.krylov
        TIME[,3]<-time.lanczos
        TIME[,4]<-time.naive
        TIME<-data.frame(TIME)
        colnames(TIME)<-namen[-5]
        # Sigmahat
        SIGMAHAT<-matrix(,R,5)
        SIGMAHAT[,1]<-sigmahat.cv
        SIGMAHAT[,2]<-sigmahat.krylov
        SIGMAHAT[,3]<-sigmahat.lanczos
        SIGMAHAT[,4]<-sigmahat.naive
        SIGMAHAT[,5]<-sigmahat.null
	SIGMAHAT<data.frame(SIGMAHAT)
	colnames(SIGMAHAT)<-namen
        #m.crash
        M.CRASH<-matrix(,R,2)
        M.CRASH[,1]<-m.crash.krylov
        M.CRASH[,2]<-m.crash.lanczos
        M.CRASH<-data.frame(M.CRASH)
        colnames(M.CRASH)<-namen[c(2,3)]
        # misc
        #DoF.complete=data.frame(DoF.complete)
        #colnames(DoF.complete)=0:m
	ME<-NULL
	if (is.null(true.coefficients)==FALSE){
		 ME<-matrix(,R,5)
		ME[,1]<-model.cv
		ME[,2]<-model.krylov
		ME[,3]<-model.lanczos
		ME[,4]<-model.naive
		ME[,5]<-model.null
		ME<-data.frame(ME)
        	colnames(ME)=namen	
	}
        
        

return(list(criterion=criterion,MSE=MSE,M=M,DoF=DoF,TIME=TIME,M.CRASH=M.CRASH,ME=ME,SIGMAHAT=SIGMAHAT))
}

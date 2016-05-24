bestfit.search.exhaustive.true<-function(statevar,lagstate,covariate,P,R,Q,indexBGlobal,indexCGlobal,...){

bestGlobalB<-matrix(0,nrow=P,ncol=P)
bestGlobalC<-matrix(0,nrow=P,ncol=R)

ttinfo<-"-\n-"
ttinfo2<-""
tt<-tktoplevel()
tkwm.geometry(tt,"300x200-30+30")
tkwm.title(tt,"Search Progress")
ttstatus<-tklabel(tt,text=ttinfo)
ttstatus2<-tklabel(tt,text=ttinfo2)
pb<-ttkprogressbar(parent=tt,orient="horizontal",
	length=200,mode="determinate",maximum=P,value=0)
pb2<-ttkprogressbar(parent=tt,orient="horizontal",
	length=100,mode="determinate",maximum=100,value=0)
tkgrid(tklabel(tt,text=""))
tkgrid(ttstatus)
tkgrid(tklabel(tt,text=""))
tkgrid(pb)
tkgrid(pb2)
tkgrid(ttstatus2)
tkgrid.columnconfigure(tt,pb,minsize=300)

all.models<-vector("list",P)

for(dv in 1:P){	#  *************** START Main loop ***************
# (for each of the variates...)

indexB<-indexBGlobal[dv,]
indexC<-indexCGlobal[dv,]
indexBC<-c(indexB,indexC)

E<-matrix(0,nrow=Q,ncol=1)	# residuals
A<-matrix(0,nrow=1,ncol=1)	# intercepts for the variates
B<-matrix(0,nrow=1,ncol=P)	# parameters for the variates
C<-matrix(0,nrow=1,ncol=R)	# parameters for the covariates
Yhat<-matrix(0,nrow=Q,ncol=1)	# estimates

# varind = matrix for the variates : which interactions to include
varind<-indexB
varind[varind==.5]<-1

# covarind = matrix for the covariates : which interactions to include
covarind<-indexC
covarind[covarind==.5]<-1

Y<-statevar[,dv]						# dependent variate
lv<-which(varind==1)
cv<-which(covarind==1)

nvar<-sum(varind)

X<-cbind(1,lagstate[,lv],covariate[,cv])		# predictors
X<-X[,-1,drop=F]

ttinfo<-paste("Searching for best-fit model for\nvariate",dv,"of",P,"...")
tkconfigure(pb,value=(dv-1))
tkconfigure(pb2,maximum=2^ncol(X)-1)
tkconfigure(ttstatus,text=ttinfo)

###################################################
    y <- Y
    n <- length(y)
    X <- X
    p <- ncol(X)
    M <- 2^p - 1
    NumDF <- rep(1, p)
    TotalSize <- sum(NumDF)


ones<-which(indexBC[indexBC!=0]==1)
AllModels <- data.frame(AIC=NA,t(matrix(rep(NA,p))))[-1,]
if(length(ones)==0)  AllModels[1,] <- c(AIC(glm(y~1)),rep(0,p))
		
		###### LOOP
		for (i in 1:M) {
		# <<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>
			ttinfo2<-paste(i, "\nof", M, "possible models")
			tkconfigure(pb2,value=(i))
			tkconfigure(ttstatus2,text=ttinfo2)
		# <<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>
            vars <- as.logical(rev(to.binary(i, p)))

		if(!all(vars[ones])) next

            k <- sum(vars)
            Xi <- as.matrix(X[, vars, drop = FALSE])
		Xi <- cbind(1,Xi)

            #ans <- glm(y ~ ., data = Xi)

# least squares estimates
beta<-solve( (t(Xi)%*%Xi),(t(Xi)%*%y) )

# calculate residuals for the dependent variate
yhat<-Xi%*%beta
e<-y-yhat

# parameter estimates
a<-beta[1,1]

# calculate log-likelihood
sigma<-t(e)%*%e/nrow(Xi)
lnlike<- -nrow(Xi)*(1/2)*log(2*pi)-(nrow(Xi)/2)*log(det(sigma))-nrow(Xi)/2

par<-ncol(Xi)+1

            vars[vars==T]<-beta[-1,]
			Criterion <-  -2*lnlike+2*par
			tvars<-data.frame(t(vars))
			names(tvars)<-names(AllModels)[-1]
            AllModels <- rbind(AllModels,data.frame(AIC=Criterion,tvars))

        }
		###### END LOOP

## get coefficients of lowest AIC model
AllModels<-AllModels[order(AllModels[,1]),]
indexBC[indexBC!=0]<-as.numeric(AllModels[1,-1])
indexB<-indexBC[1:P]
indexC<-indexBC[(P+1):length(indexBC)]

bestGlobalB[dv,]<-as.numeric(indexB)
bestGlobalC[dv,]<-as.numeric(indexC)

			all.models[[dv]]<-matrix(c(999,indexBGlobal[dv,],indexCGlobal[dv,]),
				nrow=nrow(AllModels),ncol=P+R+1,byrow=T)
			all.models[[dv]][all.models[[dv]]!=0]<-as.matrix(AllModels)

		
}			#  **************** END Main loop ****************


tkconfigure(pb,value=dv)
ttinfo<-"COMPLETE\n - "
ttinfo2<-"-"
tkconfigure(ttstatus,text=ttinfo)
tkconfigure(ttstatus2,text=ttinfo2)
Sys.sleep(2)
tkdestroy(tt)

list(bestGlobalB=bestGlobalB,bestGlobalC=bestGlobalC,all.models=all.models)

}


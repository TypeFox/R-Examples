fitting<-function(mydata,n,property)
	{
	dim<-dim(mydata)
	print("###### linear regression  ######")
	print("MLR - original data")
	switch(n,
	## 1 desc
fit <- lm(mydata[,1] ~ mydata[,2]),
	## 2 desc	
fit <- lm(mydata[,1] ~ mydata[,2] + mydata[,3]),
	## 3 desc	
fit <- lm(mydata[,1] ~ mydata[,2] + mydata[,3] + mydata[,4]),
	## 4 desc	
fit <- lm(mydata[,1] ~ mydata[,2] + mydata[,3] + mydata[,4] + mydata[,5]),
	## 5 desc	
fit <- lm(mydata[,1] ~ mydata[,2] + mydata[,3] + mydata[,4] + mydata[,5] + mydata[,6]),
	## 6 desc	
fit <- lm(mydata[,1] ~ mydata[,2] + mydata[,3] + mydata[,4] + mydata[,5] + mydata[,6]+ mydata[,7]),
	## 7 desc	
fit <- lm(mydata[,1] ~ mydata[,2] + mydata[,3] + mydata[,4] + mydata[,5] + mydata[,6]+ mydata[,7] + mydata[,8]),
	)

	print(summary(fit))
	print(names(mydata))
	R<-cor(mydata[,1], fit$fitted.values)**2
	
	## writting file
	mat <- matrix(c(mydata[,1],fit$fitted.values), nrow=dim[1], ncol=2)
	colnames(mat) <- c(paste('Experimental',property,sep=" "), paste('Predicted',property,sep=" "))
	write.table(mat, paste('prediction_TrainSet_',property,'.csv',sep=""), sep = ";")

	## Creation Graph
	tiff(paste(property,'_TrainingSet','.tiff',sep=""), res=400, width = 8, height = 8, units = "cm", pointsize=8)
	plot(mydata[,1],fit$fitted.values, xlab=paste('Experimental',property,sep=" ") , ylab=paste('Predicted',property,sep=" "), col="blueviolet", pch=16, xlim=c(min(fit$fitted.values,mydata[,1]),max(fit$fitted.values,mydata[,1])), ylim=c(min(fit$fitted.values,mydata[,1]),max(fit$fitted.values,mydata[,1])))
	abline(0,1)
	dev.off()

	print("# R^2 model #")
	print(R)

	# Calculation of RMSE
	rss<-sum(fit$residuals**2)
	rmse<-sqrt(rss/(dim[1]-(dim(as.table(fit$coefficients))
-1)-1))
	print("# RMSE model #")
	print(rmse)

	return(fit)
	}

graphe_3Sets <- function(fit,mydata,mynewdata,mynewdata2,n)
	{
	print("###### Prediction ######")
	pred<-NULL
	pred2<-NULL
	var<-fit$coeff
	rmse<-0

	# Prediction of test/validation set
	for(i in 1:dim(mynewdata)[1])
		{
			switch(n,
	## 1 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]),
	## 2 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]+ var[3]*mynewdata[i,3]),
	## 3 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]+ var[3]*mynewdata[i,3] + var[4]*mynewdata[i,4]),
	## 4 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]+ var[3]*mynewdata[i,3] + var[4]*mynewdata[i,4] + var[5]*mynewdata[i,5]),
	## 5 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]+ var[3]*mynewdata[i,3] + var[4]*mynewdata[i,4] + var[5]*mynewdata[i,5]+ var[6]*mynewdata[i,6]),
	## 6 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]+ var[3]*mynewdata[i,3] + var[4]*mynewdata[i,4] + var[5]*mynewdata[i,5]+ var[6]*mynewdata[i,6]+ var[7]*mynewdata[i,7]),
	## 7 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]+ var[3]*mynewdata[i,3] + var[4]*mynewdata[i,4] + var[5]*mynewdata[i,5]+ var[6]*mynewdata[i,6]+ var[7]*mynewdata[i,7]+ var[8]*mynewdata[i,8])
				)
		pred<-c(pred,y_pred)
		}

	# Prediction of external test set 
	for(i in 1:dim(mynewdata2)[1])
		{
			switch(n,
	## 1 desc	
y_pred2<-(var[1] + var[2]*mynewdata2[i,2]),
	## 2 desc	
y_pred2<-(var[1] + var[2]*mynewdata2[i,2]+ var[3]*mynewdata2[i,3]),
	## 3 desc	
y_pred2<-(var[1] + var[2]*mynewdata2[i,2]+ var[3]*mynewdata2[i,3] + var[4]*mynewdata2[i,4]),
	## 4 desc	
y_pred2<-(var[1] + var[2]*mynewdata2[i,2]+ var[3]*mynewdata2[i,3] + var[4]*mynewdata2[i,4] + var[5]*mynewdata2[i,5]),
	## 5 desc	
y_pred2<-(var[1] + var[2]*mynewdata2[i,2]+ var[3]*mynewdata2[i,3] + var[4]*mynewdata2[i,4] + var[5]*mynewdata2[i,5]+ var[6]*mynewdata2[i,6]),
	## 6 desc	
y_pred2<-(var[1] + var[2]*mynewdata2[i,2]+ var[3]*mynewdata2[i,3] + var[4]*mynewdata2[i,4] + var[5]*mynewdata2[i,5]+ var[6]*mynewdata2[i,6]+ var[7]*mynewdata2[i,7]),
	## 7 desc	
y_pred2<-(var[1] + var[2]*mynewdata2[i,2]+ var[3]*mynewdata2[i,3] + var[4]*mynewdata2[i,4] + var[5]*mynewdata2[i,5]+ var[6]*mynewdata2[i,6]+ var[7]*mynewdata2[i,7]+ var[8]*mynewdata2[i,8])
				)
		pred2<-c(pred2,y_pred2)
		}

	Rext<-cor(mynewdata[,1], pred)**2
	print("Validation set")
	print("# R^2 ext #")
	print(Rext)

	Rext2<-cor(mynewdata2[,1], pred2)**2
	print("External validation set")
	print("# R^2 ext - jeu externe #")
	print(Rext2)

	## Creation Graph
	tiff("Graphe_3sets.tiff", res=400, width = 8, height = 8, units = "cm", pointsize=8)

	plot(mynewdata[,1], pred, xlab="Experimental Y" , ylab="Predicted Y", col="salmon", pch=2 , xlim=c(min(fit$fitted.values,pred,pred2),max(fit$fitted.values,pred,pred2)), ylim=c(c(min(fit$fitted.values,pred2),max(fit$fitted.values,pred2))))
	points(mydata[,1], fit$fitted.values, col="blueviolet", pch=16)
	points(mynewdata2[,1], pred2, col="chartreuse4", pch=22)
	abline(0,1)
	legend("topleft", c("Training set", "Test set","External set"), pch=c(16,2,22), col=c("blueviolet","salmon","chartreuse4"), bty="n")
	dev.off()
	
	return(c(Rext,Rext2))
	}


LMO <- function(mydata,cv,n)
	{
	pred<-NULL
	print("### Leave many out validation ###")
	dim_y<-dim(mydata)[1]

	## Calculation by group of predicted y  
	for(k in 1:cv)
		{
		list<-NULL
		for(i in seq(k,dim_y,cv)) {list<-c(list,i)}	# decoupage
		print(list)

			switch(n,
		## x=desc[-list,], y=mydata[list,1]
		## 1 desc	
fit <- lm(mydata[-list,1] ~ mydata[-list,2]),
		## 2 desc	
fit <- lm(mydata[-list,1] ~ mydata[-list,2]+mydata[-list,3]),
		## 3 desc	
fit <- lm(mydata[-list,1] ~ mydata[-list,2]+mydata[-list,3]+mydata[-list,4]),
		## 4 desc	
fit <- lm(mydata[-list,1] ~ mydata[-list,2]+mydata[-list,3]+mydata[-list,4]+mydata[-list,5]),
		## 5 desc	
fit <- lm(mydata[-list,1] ~ mydata[-list,2]+mydata[-list,3]+mydata[-list,4]+mydata[-list,5]+mydata[-list,6]),
		## 6 desc	
fit <- lm(mydata[-list,1] ~ mydata[-list,2]+mydata[-list,3]+mydata[-list,4]+mydata[-list,5]+mydata[-list,6] +mydata[-list,7]),
		## 7 desc	
fit <- lm(mydata[-list,1] ~ mydata[-list,2]+mydata[-list,3]+mydata[-list,4]+mydata[-list,5]+mydata[-list,6] +mydata[-list,7]+ mydata[-list,8])
				)

		var<-fit$coeff

			for(i in list)
			{
				switch(n,
			## 1 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]),
			## 2 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]+ var[3]*mydata[i,3]),
			## 3 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]+ var[3]*mydata[i,3] + var[4]*mydata[i,4]),
			## 4 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]+ var[3]*mydata[i,3] + var[4]*mydata[i,4]+ var[5]*mydata[i,5]),
			## 5 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]+ var[3]*mydata[i,3] + var[4]*mydata[i,4] + var[5]*mydata[i,5]+ var[6]*mydata[i,6]),
			## 6 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]+ var[3]*mydata[i,3] + var[4]*mydata[i,4] + var[5]*mydata[i,5]+ var[6]*mydata[i,6]+ var[7]*mydata[i,7]),
			## 7 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]+ var[3]*mydata[i,3] + var[4]*mydata[i,4] + var[5]*mydata[i,5]+ var[6]*mydata[i,6]+ var[7]*mydata[i,7]+ var[8]*mydata[i,8])
				)

			pred[i]<-y_pred
			}
		}

	# Calculation of mean quadratic error
	res<-0
	A<-0
	for(i in 1:dim_y)
		{
		res<-(res + (mydata[i,1]-pred[i])**2)
		A<-(A + (mydata[i,1]-mean(mydata[,1]))**2)
		}
	PRESS<-mean(res)

	q<-(1-(PRESS/A))
	print("Number of groups")
	print(cv)
	print("Q^2 LMO")
	print(q)

	return(q)
	}
LOO <- function(mydata,n)
	{
	pred<-NULL
	print("### Leave one out validation ###")
	dim_y<-dim(mydata)[1]

	for(i in 1:dim_y)
		{
			switch(n,
		# Re-fitting of the group without the k molecule
	## 1 desc	
fit <- lm(mydata[-i,1] ~ mydata[-i,2]),
	## 2 desc	
fit <- lm(mydata[-i,1] ~ mydata[-i,2]+mydata[-i,3]),
	## 3 desc	
fit <- lm(mydata[-i,1] ~ mydata[-i,2]+mydata[-i,3]+mydata[-i,4]),
	## 4 desc	
fit <- lm(mydata[-i,1] ~ mydata[-i,2]+mydata[-i,3]+mydata[-i,4]+mydata[-i,5]),
	## 5 desc	
fit <- lm(mydata[-i,1] ~ mydata[-i,2]+mydata[-i,3]+mydata[-i,4]+mydata[-i,5]+mydata[-i,6]),
	## 6 desc	
fit <- lm(mydata[-i,1] ~ mydata[-i,2]+mydata[-i,3]+mydata[-i,4]+mydata[-i,5]+mydata[-i,6]+mydata[-i,7]),
	## 7 desc	
fit <- lm(mydata[-i,1] ~ mydata[-i,2]+mydata[-i,3]+mydata[-i,4]+mydata[-i,5]+mydata[-i,6]+mydata[-i,7]+mydata[-i,8])
				)
		# Prediction for the molecule k using the new equation
		var<-fit$coeff
			switch(n,
	## 1 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]),
	## 2 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]+ var[3]*mydata[i,3]),
	## 3 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]+ var[3]*mydata[i,3] + var[4]*mydata[i,4]),
	## 4 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]+ var[3]*mydata[i,3] + var[4]*mydata[i,4] + var[5]*mydata[i,5]),
	## 5 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]+ var[3]*mydata[i,3] + var[4]*mydata[i,4] + var[5]*mydata[i,5]+ var[6]*mydata[i,6]),
	## 6 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]+ var[3]*mydata[i,3] + var[4]*mydata[i,4] + var[5]*mydata[i,5]+ var[6]*mydata[i,6]+ var[7]*mydata[i,7]),
	## 7 desc	
y_pred<-(var[1] + var[2]*mydata[i,2]+ var[3]*mydata[i,3] + var[4]*mydata[i,4] + var[5]*mydata[i,5]+ var[6]*mydata[i,6]+ var[7]*mydata[i,7]+ var[8]*mydata[i,8])
				)
		pred<-c(pred,y_pred)
		}

	# Calculation of mean quadratic error
	res<-0
	A<-0
	for(i in 1:dim_y)
		{
		res<-(res + (mydata[i,1]-pred[i])**2)
		A<-(A + (mydata[i,1]-mean(mydata[,1]))**2)
		}
	PRESS<-mean(res)

	q<-(1-(PRESS/A))
	print("Q2 LOO")
	print(q)

	return(q)
	}

prediction <- function(fit,mydata,mynewdata,n)
	{
	print("###### Prediction ######")
	pred<-NULL
	var<-fit$coeff

	# Prediction for the molecule k using the new equation
	for(i in 1:dim(mynewdata)[1])
		{
			switch(n,
	## 1 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]),
	## 2 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]+ var[3]*mynewdata[i,3]),
	## 3 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]+ var[3]*mynewdata[i,3] + var[4]*mynewdata[i,4]),
	## 4 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]+ var[3]*mynewdata[i,3] + var[4]*mynewdata[i,4] + var[5]*mynewdata[i,4]),
	## 5 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]+ var[3]*mynewdata[i,3] + var[4]*mynewdata[i,4] + var[5]*mynewdata[i,5]+ var[6]*mynewdata[i,6]),
	## 6 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]+ var[3]*mynewdata[i,3] + var[4]*mynewdata[i,4] + var[5]*mynewdata[i,5]+ var[6]*mynewdata[i,6]+ var[7]*mynewdata[i,7]),
	## 7 desc	
y_pred<-(var[1] + var[2]*mynewdata[i,2]+ var[3]*mynewdata[i,3] + var[4]*mynewdata[i,4] + var[5]*mynewdata[i,5]+ var[6]*mynewdata[i,6]+ var[7]*mynewdata[i,7]+ var[8]*mynewdata[i,8])
				)
		pred<-c(pred,y_pred)
		}

	## writting file
	mat <- matrix(c(mynewdata[,1],pred), nrow=dim(mynewdata)[1], ncol=2)
	colnames(mat) <- c("Y_exp", "Y_pred")
	write.table(mat, "predictions_TestSet.csv", sep = ";")

	Rext<-cor(mynewdata[,1], pred)**2
	#print(pred)
	print("# R^2 ext #")
	print(Rext)

	# Calculation of RMSE
        residuals<-(mynewdata[,1]-pred)
	rss<-sum(residuals**2)
	rmse<-sqrt(rss/(dim[1]-(dim(as.table(fit$coefficients))
-1)-1))
	print("# RMSE ext #")
	print(rmse)

	## Creation Graph
	tiff(paste("Exp.vs.Pred.tiff"), res=400, width = 8, height = 8, units = "cm", pointsize=8)

	plot(mynewdata[,1], pred, xlab="Experimental Y" , ylab="Predicted Y", col="salmon", pch=2 , xlim=c(min(fit$fitted.values,pred,mydata[,1],mynewdata[,1]),max(fit$fitted.values,pred,mydata[,1],mynewdata[,1])), ylim=c(min(fit$fitted.values,pred,mydata[,1],mynewdata[,1]),max(fit$fitted.values,pred,mydata[,1],mynewdata[,1])))
	points(mydata[,1], fit$fitted.values, col="blueviolet", pch=16)
	abline(0,1)
	legend("topleft", c("Training set", "Test set"), pch=c(16,2), col=c("blueviolet","salmon"), bty="n")
	dev.off()
	
	return(Rext)
	}

preselection<-function(desc){

        # variance 
	print("Analysis variance")
	ColVar <- apply(desc, 2, var)

	# Remove descriptors with missing values
	d<-desc[,!is.na(ColVar)]

	# Select descriptor with value more than 0.0001
	ColVar <- apply(d, 2, var)
	desc<-d[,(ColVar>0.001) ]

        #print(desc)
	return(desc)
}
scramb <- function(mydata,k,n,cercle=FALSE)
	{
	print("######### Y-Scrambling #########")
	Ry<-c()
	Rmodel<-c()
        ys<-c()

	## Fitting
	#print("MLR - original data")
				switch(n,
		## 1 desc	
fit <- lm(mydata[,1]  ~ mydata[,2]  ),
		## 2 desc	
fit <- lm(mydata[,1]  ~ mydata[,2]  + mydata[,3] ),
		## 3 desc	
fit <- lm(mydata[,1]  ~ mydata[,2]  + mydata[,3]  + mydata[,4] ),
		## 4 desc	
fit <- lm(mydata[,1]  ~ mydata[,2]  + mydata[,3]  + mydata[,4]  + mydata[,5] ),
		## 5 desc	
fit <- lm(mydata[,1]  ~ mydata[,2]  + mydata[,3]  + mydata[,4]  + mydata[,5]  + mydata[,6]  ),
		## 6 desc	
fit <- lm(mydata[,1]  ~ mydata[,2]  + mydata[,3]  + mydata[,4]  + mydata[,5]  + mydata[,6] + mydata[,7] ),
		## 7 desc	
fit <- lm(mydata[,1]  ~ mydata[,2]  + mydata[,3]  + mydata[,4]  + mydata[,5]  + mydata[,6] + mydata[,7]+ mydata[,8] )
					)
	R<-cor(mydata[,1] , fit$fitted.values)**2
	
	for (i in 1:k) 
		{
		#### Resampling
		s<-sample(mydata[,1] )
		r2<-cor(s,mydata[,1] )**2
		Ry[i]<-r2

		#### Fit
				switch(n,
		## 1 desc	
fitBis <- lm(s ~ mydata[,2]  ),
		## 2 desc	
fitBis <- lm(s ~ mydata[,2]  + mydata[,3] ),
		## 3 desc	
fitBis <- lm(s ~ mydata[,2]  + mydata[,3]  + mydata[,4] ),
		## 4 desc	
fitBis <- lm(s ~ mydata[,2]  + mydata[,3]  + mydata[,4]  + mydata[,5] ),
		## 5 desc	
fitBis <- lm(s ~ mydata[,2]  + mydata[,3]  + mydata[,4]  + mydata[,5]  + mydata[,6]  ),
		## 6 desc	
fitBis <- lm(s ~ mydata[,2]  + mydata[,3]  + mydata[,4]  + mydata[,5]  + mydata[,6] + mydata[,7]  ),
		## 7 desc	
fitBis <- lm(s ~ mydata[,2]  + mydata[,3]  + mydata[,4]  + mydata[,5]  + mydata[,6] + mydata[,7] + mydata[,8] )
					)
		R2<-cor(s, fitBis$fitted.values)**2
		Rmodel[i]<-R2
		}

	Rmodel[k+1]<-R
	Ry[k+1]<-1

	## writting file
	mat <- matrix(c(Ry,Rmodel), nrow=k+1, ncol=2)
	colnames(mat) <- c("R^2 Yrandom/Yexp", "R^2  random")
	write.table(mat, "Scramb.csv", sep = ";")

	## Creation of Graph
	tiff("Scramb.tiff", res=400, width = 7.9, height = 7.9, units = "cm", pointsize=8)
	plot(Ry,Rmodel, pch=4 ,xlab=(expression(paste("R"^"2","(Y"["random"],"/Y"["exp"],")",sep=""))) , ylab=(expression("R"^"2"["random"])), ylim=c(0,1))
        if(cercle==TRUE){
            points(x=Ry[k+1], y=Rmodel[k+1], pch = 1, col = "forestgreen", cex = 2.5,lwd = 1.5)
        }
	dev.off()

	## Mean and Standard deviation
	print("Mean R^2 new model")
	moyenne<-mean(Rmodel[1:k])
	print(moyenne)
	print("Standard deviation R^2 new model")
	ecart<-sd(Rmodel[1:k])
	print(ecart)

	ys<-c(moyenne,ecart)

	return(ys)
	}


select_variables<-function(id,y,d,ThresholdInterCor,auto=FALSE)
	{
        ### correlation/covariance matrix
	print("Intercorrelation")
	mat_cor<-as.matrix(cor(d, use="complete.obs", method="kendall"))
	cor_y<-as.vector(cor(y,d))

        ### deletion/removal of intercorrelated variables 
	col<-dim(d)[2]
	ligne<-dim(d)[1]
	l<-NULL # list of desc to delete
	desc<-NULL

	for( i in 1:(col-1)){	
                for( j in (i+1):col){
		    if( mat_cor[i,j] >= ThresholdInterCor){
                        if( (is.element(i,l)==FALSE) ){  
                            print("Intercorrelation between:")
                            print(paste("1)",names(d[i]), sep=" ")  )
                            print(paste("    correlation with property =",signif(cor_y[i],5), sep=" "))
                            print(paste("2)",names(d[j]), sep=" ") )
			    print(paste("    correlation with property =",signif(cor_y[j],5), sep=" "))
                            print("Select 1 or 2 to choose the variable to keep for next step")

                            if(auto == TRUE){
                                if(abs(cor_y[i])  > abs(cor_y[j])){ 
                                print(paste("Select 1",names(d[i]), sep=" "))
                                l<-c(l,j)
                                }
                                else{
                                print(paste("Select 2",names(d[j]), sep=" "))
                                l<-c(l,i)
                                }
                            }

                            else{
                                n<-scan(file="", what=integer(),nmax=1)
                                if(n==1) { l<-c(l,j) } 
                                     # add j to the l list of desc to remove     
                                else if(n==2) { l<-c(l,i) } 	
                                else if(n!=1 && n!=2){
                                    stop("\n ERROR: select 1 or 2 ! \n Run Select_variables a second time\n")
                                }
                            }	
	
                        }
                    }
		}
        }

	if(is.null(l)){desc<-d}	else {desc<-d[,-l]}
	listdesc<-as.data.frame(desc)

        ###  writting file
	mydata<-cbind(id,y,listdesc)
	write.table(mydata, 'variables_selected.csv', sep = ";")

	return(listdesc)
}


### Selection of descriptors and model development ###
select_MLR<-function(y,desc,n,method="forward") 
           #function(y,desc,n) 
	{
	variables<-NULL
	liste_var<-NULL

        # Control of the value option method
        if( !(method %in% c("backward", "forward", "seqrep"))){stop("unknown method")}
# "exhaustive", ???

	print("Model's development: Selection of descriptors")
        print(method)

	# Model selection by exhaustive search, forward or backward stepwise, or sequential replacement
       # if(method == "exhaustive"){really.big=TRUE}
       # else{really.big=FALSE}
	really.big<-FALSE

	a<-regsubsets(x=desc, y=y, nbest=1,nvmax=n,intercept=TRUE, method, weights=rep(1, length(y)),force.in=NULL, force.out=NULL,really.big ) 

	# summary file with of all MLR
	con <- file("_MLR", open = "w")
	for (i in 1:n){
		delta<-(summary(a)$rsq[i]-summary(a)$rsq[i-1] )
		cat("n =",i,"\nR^2 = ",summary(a)$rsq[i],"\nA(R^2) = ", delta ,"\n",names(coef(a,i)),"\n\n   ------------ \n\n", file = con)
		}
	close(con)


	for (i in 1:n) 
		{
		# If Delta R^2 < 0.02 then we choose the previous variables
		### also if Q^2 decrease
		delta<-(summary(a)$rsq[i+1]-summary(a)$rsq[i])
		if (delta < 0.02)	
			{
			variables<-names(coef(a,i))
			n<-length(variables)-1
			print(paste("Number of variables :", n, sep=" "))
			print(paste("Delta R^2 (n/n+1) =",signif(delta,4), sep=" "))
			print(coef(a,i))
			print("R^2 and ajusted R^2 ") 
			print(summary(a)$rsq[i])	# coeff of determination
			print(summary(a)$adjr2[i])	# ajusted coef of determination 
			k<-i
			break
			}
		} 


	for( i in 2:(dim(desc)[2]+1) ){
		if(summary(a)$which[n,][i] == TRUE){
			#print( names(summary(a)$which[n,][i]) )
			#print(i)
			liste_var<-c( liste_var,(i-1))
			}
		}

	return(desc[liste_var])
	}

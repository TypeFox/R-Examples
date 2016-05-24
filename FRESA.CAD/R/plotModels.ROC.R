plotModels.ROC <-
function(modelPredictions,number.of.models=0,specificities=c(0.975,0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05),theCVfolds=1,predictor="Prediction",...) 
{

	number.of.runs=number.of.models;
	par(mfrow=c(1,1),pty='s');
	rocadded = 0;
	if (number.of.runs == 0)
	{
		number.of.runs=max(modelPredictions[,"Model"]) %/% theCVfolds;
	}
	auclist <- vector()
	sumSen <- NULL;
	blindSen <- NULL;
	

	for (n in 1:number.of.runs)
	{
		mm = n;
		blindmodel <- modelPredictions[which(((modelPredictions[,3]-1) %/% theCVfolds) + 1  == mm),];
		if ( (sum(blindmodel[,"Outcome"]==1) > 3) && (sum(blindmodel[,"Outcome"]==0) > 3))
		{
			auclist <- append(auclist,pROC::roc(blindmodel[,"Outcome"],blindmodel[,predictor],auc=TRUE,plot=TRUE,col="lightgray",lty=4,lwd=1)$auc)
			par(new=TRUE)
			sen <- pROC::roc(blindmodel[,"Outcome"],blindmodel[,predictor],auc=TRUE,plot=FALSE,ci=TRUE,progress='none',of='se',specificities=specificities,boot.n=100,smooth=FALSE,lty=3,lwd=1)$ci[,2]
			if (n == 1) 
			{
				blindSen <- sen;
			}
			else
			{
				blindSen <- rbind(sen,blindSen);
			}
			rocadded = rocadded +1;
		}
	}
	if (rocadded>1)
	{
		par(new=TRUE)
		boxplot(blindSen,add=TRUE, axes = FALSE,boxwex=0.04,at=specificities);
		par(new=TRUE)
	}
	psta <- boxplot(modelPredictions[,predictor]~rownames(modelPredictions),plot=FALSE)
	outcomesta <- boxplot(modelPredictions$Outcome~rownames(modelPredictions),plot=FALSE)
	ensemble.auc <- pROC::roc(outcomesta$stats[3,],psta$stats[3,],col="black",auc=TRUE,plot=TRUE,smooth=FALSE,lty=3,lwd=3,...)$auc
	par(new=TRUE)
	auc1 <- pROC::roc(modelPredictions[,"Outcome"],modelPredictions[,predictor],col="darkblue",auc=TRUE,plot=TRUE,smooth=FALSE,...)$auc

	ensemblePrediction <- cbind(outcomesta$stats[3,],psta$stats[3,]);
	rownames(ensemblePrediction) <- psta$names
	colnames(ensemblePrediction) <- c("Outcome",predictor);
	thres = 0.5
	if (min(psta$stats[3,])<0) 
	{
		thres = 0;
	}				
	dtable <- table(psta$stats[3,]<thres,1-outcomesta$stats[3,])
	colnames(dtable) <- c("O(+)","O(-)")
	rownames(dtable) <- c("P(+)","P(-)")
	Sen=dtable[1,1]/(dtable[1,1]+dtable[2,1])
	Spe=dtable[2,2]/(dtable[1,2]+dtable[2,2])
	enauc = 0.5*(Sen+Spe)	

	ley.names <- c(paste("Coherence (",sprintf("%.3f",auc1),")"));
	ley.colors <- c("darkblue");
	ley.lty <- c(1);
	ley.names <- append(ley.names,paste("Ensemble (",sprintf("%.3f",ensemble.auc),")"));
	ley.colors <- append(ley.colors,"black");
	ley.lty <- append(ley.lty,3);

    lines(c(1,Spe,0),c(0,Sen,1),col="green",lwd=1.0,lty=1);
    
    ley.names <- append(ley.names,paste("Ensemble Prediction (",sprintf("%.3f",enauc),")"));
    ley.colors <- append(ley.colors,"green");
    ley.lty <- append(ley.lty,1);

	auc = 0;
	if (rocadded>1)
	{
#		boxplot(blindSen,add=TRUE, axes = FALSE,boxwex=0.04,at=specificities);
		sumSen <- colMeans(blindSen,na.rm = TRUE);
		sennames <- names(sumSen);
		sumSen <- append(0,sumSen);
		sumSen <- append(sumSen,1);
		sennames <- append("1",sennames);
		sennames <- append(sennames,"0");
		names(sumSen) <- sennames;
		spevalues <- as.numeric(names(sumSen));
		lines(spevalues,sumSen,col="red",lwd=1.0,lty=2);
		for (i in 2:length(spevalues))
		{
			auc = auc + (spevalues[i-1]-spevalues[i])*(sumSen[i-1]+(sumSen[i]-sumSen[i-1])/2)
		}
		ley.names <- append(ley.names,c("Fold ROCs",paste("Mean Sensitivities(",sprintf("%.3f",auc),")")));
		ley.colors <- append(ley.colors,c("lightgray","red"));
		ley.lty <- append(ley.lty,c(4,2));
	}

	legend(0.7,0.20, legend=ley.names,col = ley.colors, lty = ley.lty,bty="n")

    x=1.025
    y=1.025
    Acc=(dtable[1,1]+dtable[2,2])/(dtable[1,1]+dtable[2,2]+dtable[1,2]+dtable[2,1])
    F1=(2*dtable[1,1])/(2*dtable[1,1]+dtable[1,2]+dtable[2,1])
    text(x,y,paste("(TPR=",sprintf("%.3f",Sen),",TNR=",sprintf("%.3f",Spe),",ACC=",sprintf("%.3f",Acc),",F1=",sprintf("%.3f",F1),",AUC=",sprintf("%.3f",enauc),")"),adj = c(0,1),cex=0.7,col="dark green")
	par(new=TRUE,plt=c(0.6,0.8,0.37,0.57),pty='s')
	plot(t(dtable),main="Ensemble",ylab=predictor,xlab="Outcome")
	par(plt=c(0.0,1.0,0.09,0.91),pty='m')
	result <- list(ROC.AUCs=auclist,
	mean.sensitivities=sumSen,
	model.sensitivities=blindSen,
	specificities=specificities,
	senAUC=auc,
	ensemblePrediction=ensemblePrediction,
	predictionTable=dtable)
	return (result)
}

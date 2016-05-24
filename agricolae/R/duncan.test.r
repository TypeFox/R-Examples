`duncan.test` <-
		function (y, trt, DFerror, MSerror, alpha=0.05, group=TRUE,main = NULL,console=FALSE)
{
	name.y <- paste(deparse(substitute(y)))
	name.t <- paste(deparse(substitute(trt)))
	if(is.null(main))main<-paste(name.y,"~", name.t)
	clase<-c("aov","lm")
	if("aov"%in%class(y) | "lm"%in%class(y)){
		if(is.null(main))main<-y$call
		A<-y$model
		DFerror<-df.residual(y)
		MSerror<-deviance(y)/DFerror
		y<-A[,1]
		ipch<-pmatch(trt,names(A))
		nipch<- length(ipch)
		for(i in 1:nipch){
			if (is.na(ipch[i]))
				return(if(console)cat("Name: ", trt, "\n", names(A)[-1], "\n"))
		}
		name.t<- names(A)[ipch][1]
		trt <- A[, ipch]
		if (nipch > 1){
			trt <- A[, ipch[1]]
			for(i in 2:nipch){
				name.t <- paste(name.t,names(A)[ipch][i],sep=":")
				trt <- paste(trt,A[,ipch[i]],sep=":")
			}}
		name.y <- names(A)[1]
	}
	junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
	Mean<-mean(junto[,1])
	CV<-sqrt(MSerror)*100/Mean	
	means <- tapply.stat(junto[,1],junto[,2],stat="mean") # change
	sds <-   tapply.stat(junto[,1],junto[,2],stat="sd") #change
	nn <-   tapply.stat(junto[,1],junto[,2],stat="length") # change
	mi<-tapply.stat(junto[,1],junto[,2],stat="min") # change
	ma<-tapply.stat(junto[,1],junto[,2],stat="max") # change
	means<-data.frame(means,std=sds[,2],r=nn[,2],Min=mi[,2],Max=ma[,2])
	rownames(means)<- means[,1]
	means <- means[,-1]
	names(means)[1]<-name.y
#   row.names(means)<-means[,1]
	ntr<-nrow(means)
	Tprob <- qtukey((1-alpha)^(1:(ntr-1)),2:ntr, DFerror)
	nr <- unique(nn[,2])
	
#"Critical Value of Studentized Range")
	if(console){
	cat("\nStudy:", main)
	cat("\n\nDuncan's new multiple range test\nfor",name.y,"\n")
	cat("\nMean Square Error: ",MSerror,"\n\n")
	cat(paste(name.t,",",sep="")," means\n\n")
	print(means)}
	if(length(nr) == 1 ) sdtdif <- sqrt(MSerror/nr)
	else {
		nr1 <-  1/mean(1/nn[,2])
		sdtdif <- sqrt(MSerror/nr1)
	}
	DUNCAN <- Tprob * sdtdif
	names(DUNCAN)<-2:ntr
	if(console){cat("\nalpha:",alpha,"; Df Error:",DFerror,"\n")
	cat("\nCritical Range\n")
	print(DUNCAN)}
	if (length(nr) > 1) {
		if(console){cat("\nHarmonic Mean of Cell Sizes ", nr1 )
		cat("\n\nDifferent value for each comparison")}
	}
	if (group) {
		if(console)cat("\nMeans with the same letter are not significantly different.")
		if(console)cat("\n\nGroups, Treatments and means\n")
		groups <- order.group(rownames(means), means[,1], means[,3], MSerror, 1 ,means[,2], parameter=0.5,snk=2,DFerror,alpha,sdtdif,console=console)
	comparison=NULL
	groups <- data.frame(groups[,1:3])
	}
	if (!group) {
		Omeans<-order(means[,1],decreasing = TRUE)
		Ordindex<-order(Omeans)
		comb <-utils::combn(ntr,2)
		nn<-ncol(comb)
		dif<-rep(0,nn)
		DIF<-dif
		LCL<-dif
		UCL<-dif
		pvalue<-dif
		odif<-dif
		sig<-NULL
		for (k in 1:nn) {
			i<-comb[1,k]
			j<-comb[2,k]
#if (means[i, 2] < means[j, 2]){
#comb[1, k]<-j
#comb[2, k]<-i
#}
			dif[k]<-means[i,1]-means[j,1]
			DIF[k]<-abs(dif[k])
			nx<-abs(i-j)+1
			odif[k] <- abs(Ordindex[i]- Ordindex[j])+1
#sdtdif<-sqrt(MSerror * (1/means[i,4] + 1/means[j,4]))
			pvalue[k]<- round(1-ptukey(DIF[k]/sdtdif,odif[k],DFerror)^(1/(odif[k]-1)),4)
			LCL[k] <- dif[k] - DUNCAN[odif[k]-1]
			UCL[k] <- dif[k] + DUNCAN[odif[k]-1]
			sig[k]<-" "
			if (pvalue[k] <= 0.001) sig[k]<-"***"
			else  if (pvalue[k] <= 0.01) sig[k]<-"**"
			else  if (pvalue[k] <= 0.05) sig[k]<-"*"
			else  if (pvalue[k] <= 0.1) sig[k]<-"."
		}
		tr.i <- rownames(means)[comb[1,] ]
		tr.j <- rownames(means)[comb[2,] ]
		comparison<-data.frame("Difference" = dif, pvalue=pvalue,"sig."=sig,LCL,UCL)
		rownames(comparison)<-paste(tr.i,tr.j,sep=" - ")
		if(console){cat("\nComparison between treatments means\n\n")
		print(comparison)}
		groups=NULL
#		output<-data.frame(trt= means[,1],means= means[,2],M="",N=means[,4],std.err=means[,3])
	}
	parameters<-data.frame(Df=DFerror,ntr = ntr,alpha=alpha,test="Duncan",name.t=name.t)
	statistics<-data.frame(Mean=Mean,CV=CV,MSerror=MSerror)
	Duncan<-data.frame(Table=Tprob,CriticalRange=DUNCAN)
	rownames(parameters)<-" "
	rownames(statistics)<-" "
	output<-list(statistics=statistics,parameters=parameters, Duncan=Duncan,
	means=means,comparison=comparison,groups=groups)
	invisible(output)
}

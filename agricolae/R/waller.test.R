`waller.test` <-                                                                                     
		function (y, trt, DFerror, MSerror, Fc, K = 100, group=TRUE,main = NULL,console=FALSE)                             
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
		tabla<-anova(y)                                                                                  
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
		Fc<-tabla[name.t,4]                                                                              
		name.y <- names(A)[1]                                                                            
	}                                                                                                
	junto <- subset(data.frame(y, trt), is.na(y) == FALSE)                                           
	Mean<-mean(junto[,1])
	CV<-sqrt(MSerror)*100/Mean	
	means <- tapply.stat(junto[,1],junto[,2],stat="mean") # change                                   
	sds <-   tapply.stat(junto[,1],junto[,2],stat="sd")   # change                                   
	nn <-   tapply.stat(junto[,1],junto[,2],stat="length") # change                                  
	mi<-tapply.stat(junto[,1],junto[,2],stat="min") # change
	ma<-tapply.stat(junto[,1],junto[,2],stat="max") # change
	means<-data.frame(means,std=sds[,2],r=nn[,2],Min=mi[,2],Max=ma[,2])
	names(means)[1:2]<-c(name.t,name.y)                                                              
#    row.names(means)<-means[,1]                                                                     
	ntr<-nrow(means)                                                                                 
	Tprob <- waller(K,ntr-1,DFerror,Fc)                                                              
	nr <- unique(nn[,2]) 
	nr1 <-nr                                                                            
	nfila<-c("K ratio", "Error Degrees of Freedom", "Error Mean Square","F value",                       
			"Critical Value of Waller")                                                                          
	nvalor<-c( K,  DFerror, MSerror, Fc, Tprob)                                                          
	if(console){
	cat("\nStudy:", main)                                                                            
	cat("\n\nWaller-Duncan K-ratio t Test for",name.y,"\n")                                          
	cat("\nThis test minimizes the Bayes risk under additive")                                       
	cat("\nloss and certain other assumptions.\n") }                                                  
	xtabla<-data.frame("......"=nvalor)                                                                  
	row.names(xtabla)<-nfila                                                                             
	if(console){print(xtabla)                                                                                        
	cat("\n")                                                                                            
	cat(paste(name.t,",",sep="")," means\n\n")                                                           
	print(data.frame(row.names = means[,1], means[,-1]))}                                                 
	if (length(nr) == 1) {                                                                           
		MSD <- Tprob * sqrt(2 * MSerror/nr)                                                          
		if(console)cat("\nMinimum Significant Difference", MSD)                                                 

	}                                                                                                
	else {                                                                                           
		nr<- 1/mean(1/nn[,2])
		MSD <- Tprob * sqrt(2 * MSerror/nr)
        if(console)cat("\nHarmonic Mean of Cell Sizes ", nr )
	}                                                                                                

	if (group) {                                                                                         
		if(console){cat("\nMeans with the same letter are not significantly different.")                                 
		cat("\n\nComparison of treatments\n\nGroups, Treatments and means\n")}                               
		groups <- order.group(means[,1], means[,2], nr, MSerror, Tprob,means[,3],console=console)                     
		comparison=NULL
		groups <- data.frame(groups[,1:3])
statistics<-data.frame(Mean=Mean,CV=CV,MSerror=MSerror,F.Value=Fc,CriticalDifference=MSD)	
	}                                                                                                    
	if (!group) {                                                                                        
		comb <-utils::combn(ntr,2)                                                                                  
		nn<-ncol(comb)                                                                                       
		dif<-rep(0,nn)                                                                                       
		MSD1<-rep(0,nn)                                                                                      
		for (k in 1:nn) {                                                                                    
			i<-comb[1,k]                                                                                         
			j<-comb[2,k]                                                                                         
#if (means[i, 2] < means[j, 2]){                                                                     
#comb[1, k]<-j                                                                                       
#comb[2, k]<-i                                                                                       
#}                                                                                                   
			dif[k]<-means[i,2]-means[j,2]                                                                        
			MSD1[k]<-Tprob*sqrt(MSerror * (1/means[i,4] + 1/means[j,4]))                                         
		}                                                                                                    
		tr.i <- means[comb[1, ],1]                                                                           
		tr.j <- means[comb[2, ],1]                                                                           
		if (length(nr1) == 1)  {                                                                              
			significant = abs(dif) > MSD                                                                         
			comparison<-data.frame("Difference" = dif, significant)
statistics<-data.frame(Mean=Mean,CV=CV,MSerror=MSerror,F.Value=Fc,CriticalDifference=MSD)
		}                                                                                                    
		else  {                                                                                              
			significant = abs(dif) > MSD1                                                                             
			comparison<-data.frame("Difference" = dif, MSD=MSD1,significant)
statistics<-data.frame(Mean=Mean,CV=CV,MSerror=MSerror,F.Value=Fc, r.harmonic=nr)
		}                                                                                                    
		rownames(comparison)<-paste(tr.i,tr.j,sep=" - ")                                                         
		if(console){cat("\nComparison between treatments means\n\n")                                                     
		print(comparison)}                                                                                       
		groups=NULL
#		output<-data.frame(trt= means[,1],means= means[,2],M="",N=means[,4],std.err=means[,3])               
	}                                                                                                    
	parameters<-data.frame(Df=DFerror,ntr = ntr, K=K,Waller=Tprob,name.t=name.t)
	rownames(parameters)<-" "
	rownames(statistics)<-" "
	rownames(means)<-means[,1]
	means<-means[,-1]
	output<-list(statistics=statistics,parameters=parameters, 
			means=means,comparison=comparison,groups=groups)
	invisible(output)                                                                                    
}                                                                                                    

computeCellFrequencyDistributions<-function(dm, max_PM=6, precision, min_CellFreq=0.1){

print("Computing cell-frequency probability distributions...")
##add cell-frequency that best explains CN and AF to dm
dm=.addColumn(dm,"f",NA);

##compute cell-freq-probability for each mutation
freq=t(seq(min_CellFreq,1,by=precision/10));
densities=matrix(matrix(NA,nrow(dm),length(freq)),nrow=nrow(dm),ncol=length(freq),dimnames=list(1:nrow(dm),freq));
success=0; 
errors=c();
warnings=c();
for (k in 1:nrow(dm)){
	output=try(cellfrequency_pdf(dm[k,"AF_Tumor"],dm[k,"CN_Estimate"],
		dm[k,"PN_B"],freq, max_PM),silent=TRUE);
	if(class(output)=="try-error"){
         	errors=rbind(errors,output);
	}else{
	        output$p=output$p/sum(output$p,na.rm=T); ##under the assumption that relative rather than absolute probabilities matter for clustering
        	densities[k,]=output$p;
        	dm[k,"f"]=output$bestF;
        	success=success+1;
	} 
        
   if (mod(k,20)==0){
       print(paste("Processed", k, "out of ",nrow(dm),"SNVs --> success: ",
           success,"/",k))
   }
}
failure=nrow(dm)-success;
if (length(errors)>0){
	print(paste("Failed to find cell-frequency distribution for ",
		failure,"SNVs. Causes:"));
	errors=unique(errors);
	for (i in 1:length(errors)){
	    print(errors[i]);
	}
}
if (length(warnings)>0){
	warnings=unique(warnings);
	print(paste("Additional warnings:"))
	for (i in 1:length(warnings)){
	    print(warnings[i]);
	}
}
print("...Done.")
output=list("densities"=densities,"freq"=freq,"dm"=dm);
return(output);
}

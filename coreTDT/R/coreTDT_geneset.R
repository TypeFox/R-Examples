#library(rvTDT)

#setwd("/Users/yujiang/network_driver/chgv3/TDT_rare_variants/ATAV_interact")

#samplePed = read.table("ATAV_interact_output_recodeA.raw",header=TRUE,stringsAsFactors=F,strip.white=T,check.names=F);
#controlInf = read.csv("temp_evs_process.rvTDT_tmp",header=F,sep=",",row.names=1, stringsAsFactors=F, strip.white=T)

coreTDT_geneset = function(samplePed,controlInf,useControlMAF=TRUE,maf.threshold = 1, qc.proportion =0.8,geneList=c(),outputFile="coreTDT_analysis.out",chrX=FALSE,writeFile=FALSE){
	
	##formalize the sample to be ordered as the first n row are samples and second are mothers and the last are fathers	
	nsample_trio = dim(samplePed)[1]
	if(nsample_trio %%3){
		stop("sample is not exact trios");
	}
	nsample = nsample_trio/3
	
	child_sample = which(samplePed$PAT !=0 & samplePed$MAT !=0)
	childID = samplePed$IID[child_sample]
	motherID =samplePed$MAT[child_sample]
	fatherID =samplePed$PAT[child_sample]  
	
	
	##find the shared variants between 	
	##clear the control file, i.e., remove loci that not covered in evs, ----------------------should we?
	controlInf.clear = controlInf[which(controlInf[,3]),]
	
	sampleVar = colnames(samplePed)[-c(1:6)]
	controlVar = rownames(controlInf.clear)
	mergedVar = intersect(sampleVar, controlVar)
	
	rownames(samplePed) = samplePed$IID
	samplePed.new = samplePed[c(childID,motherID,fatherID),mergedVar]
	controlInf.new = controlInf[mergedVar,]
	
	if(useControlMAF ==T){
		controlMAF = apply(as.matrix(controlInf.new[,c(4,5,6)]),1,function(x) (x[1]+x[2]/2)/(x[1]+x[2]+x[3]))		
	}else{
		controlMAF =NULL	
	}
	
	
	results=list()
		
		if(length(geneList)==0){
			geneList = unique(controlInf.new[,1])
		}
			
		igene=0	
		for(gene in geneList){
			geneVar = rownames(controlInf.new[controlInf.new[,1] == gene,])
			if(length(geneVar) >1){
				ped = as.matrix(samplePed.new[, geneVar])
				results[[gene]] = coreTDT_single(ped, control.maf = controlMAF[geneVar] ,maf.threshold = maf.threshold, qc.proportion =qc.proportion)
				
			}else if(length(geneVar) ==1){
				ped = matrix(samplePed.new[, geneVar],ncol=1)
				results[[gene]] = coreTDT_single(ped,control.maf = controlMAF[geneVar] ,maf.threshold = maf.threshold, qc.proportion =qc.proportion)
			}
			
			igene = igene+1
			if(igene%%100 ==0 && writeFile==TRUE){
				result.df = coreTDT.results.2.df(results)
				write.table(result.df,file=outputFile,quote=F)		
			}
					
		}
		
	
	##output the result
	result.df = coreTDT.results.2.df(results)
	
	if(writeFile==TRUE){
		write.table(result.df,file=outputFile,quote=F)
	}
	
	return(result.df)
}

coreTDT.results.2.df <- function(coreTDTresults) {
    # convert list of rvTDT results to a more
    # output-friendly dataframe of results
    
    # assumed list keys for rvTDT output
    coreTDTdf.colnames = c("nfamily","nsnp","pvalue_pr",
        "pvalue_lr","pvalue_lr2","N11","N12","N112","N122","nmissing","nMedErr");
    # get sorted character array of gene names
    genes = sort(names(coreTDTresults));
    # build up output data.frame to store rvTDT results in
    coreTDTdf.rownames = genes;
    coreTDTdf = data.frame(matrix(nrow=length(coreTDTdf.rownames),
        ncol=length(coreTDTdf.colnames)))
    # assign row and column names to rvTDT output dataframe 
    rownames(coreTDTdf) = coreTDTdf.rownames
    colnames(coreTDTdf) = coreTDTdf.colnames
    # iterate through each gene/key combo and store to
    # rvTDT output dataframe
    for (gene in coreTDTdf.rownames) {
    	if(!is.null(coreTDTresults[[gene]])){
	        for (element in coreTDTdf.colnames) {
    	        coreTDTdf[gene,element] = coreTDTresults[[gene]][[element]];
    	    }
    	 }      
    }
    return(coreTDTdf);
}

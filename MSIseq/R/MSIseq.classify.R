## classify function
MSIseq.classify <- function(mutationNum, classifier = NGSclassifier, cancerType = NULL) {
  
  ## check mutationNum file
  if(!all(c("T.sns", "S.sns", "T.ind", "S.ind", "T", "S", "ratio.sns", "ratio.ind", "ratio")%in%colnames(mutationNum))){
    	stop("Wrong column names in mutationNum.")
  }
  
  ## check the cancerType file
  if (is.null(cancerType)) {
  	classifyset = mutationNum
  }
  else {
  	
  	if(!all(c("Tumor_Sample_Barcode", "cancer_type")%in%colnames(cancerType))){
    	stop("Wrong column names in cancerType.")
  	}
  	if (nrow(mutationNum)!=nrow(cancerType) | 
  	all(rownames(mutationNum)%in%cancerType$Tumor_Sample_Barcode) == FALSE){
    	stop("Samples do not match between classification and data.")
  	}
  	
  	## match the data with cancerType
  	classifyset = cbind(mutationNum, cancerType$cancer_type
  	[match(rownames(mutationNum), cancerType$Tumor_Sample_Barcode)])
  	colnames(classifyset)[ncol(classifyset)] = "cancer_type"
  }

  ## flag likely POLE deficent samples
  POLE = (mutationNum$T.sns>60)&(mutationNum$S.ind<0.18)
  POLE = as.character(POLE)
  POLE[POLE=="FALSE"] = "No"
  POLE[POLE=="TRUE"] = "Yes"
  
  classify.result = predict(classifier, newdata=classifyset, type="class")
  result = as.data.frame(cbind(rownames(classifyset), 
  as.character(classify.result), POLE))
  colnames(result) = c("Tumor_Sample_Barcode", "MSI_status", "Likely_POLE_deficiency")
  result
} 

## re-train function
MSIseq.train <- function(mutationNum, classification, cancerType = NULL) {
  
  ## check the classification file
  if(!all(c("Tumor_Sample_Barcode", "MSI_status")%in%colnames(classification))){
    stop("Wrong column names in classification.")
  }
  ## check if MSI_status has the right format: MSI-H and Non-MSI-H
  if(!all(classification$MSI_status%in%c("Non-MSI-H", "MSI-H"))){
    stop("Only MSI-H and Non-MSI-H are accepted in MSI_status.")  
  }
  if(length(unique(classification$MSI_status))!=2){
    stop("Both MSI-H and Non-MSI-H samples should be in training set.")  
  }
  if (nrow(mutationNum)!=nrow(classification) | 
  all(rownames(mutationNum)%in%classification$Tumor_Sample_Barcode)==FALSE){
    stop("Samples do not match between classification and data.")
  }
  
  ## match the data with classification
  trainset = mutationNum
  trainset = cbind(trainset, classification$MSI_status
  [match(rownames(trainset), classification$Tumor_Sample_Barcode)])
  colnames(trainset)[ncol(trainset)] = "MSI_status"

  ## check the cancerType file
  if (!is.null(cancerType)) {
  	
  	if(!all(c("Tumor_Sample_Barcode", "cancer_type")%in%colnames(cancerType))){
    	stop("Wrong column names in cancerType.")
  	}
  	if (nrow(mutationNum)!=nrow(cancerType) | 
  	all(rownames(mutationNum)%in%cancerType$Tumor_Sample_Barcode) == FALSE){
    	stop("Samples do not match between classification and data.")
  	}
  	
  	## match the data with cancerType
  	trainset = cbind(trainset, cancerType$cancer_type
  	[match(rownames(trainset), cancerType$Tumor_Sample_Barcode)])
  	colnames(trainset)[ncol(trainset)] = "cancer_type"
  }
  
  if(nrow(trainset)<50){
    warning('training set size < 50')  
  }
  
  ## build model with trainset
  tree.model = J48(MSI_status ~ . , data=trainset)
  tree.eval = evaluate_Weka_classifier(tree.model, numFolds = 5, seed = 1)
  cat('5 fold cross validation result: ',
      unname(tree.eval$detail[1]), '\n', sep='')
  tree.model
}

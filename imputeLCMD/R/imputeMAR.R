# .....................................................................................
# - this function performs missing values imputation under MAR/MCAR hypothesis
# - the imputation of MVs is performed for each protein containing MAR/MCAR missing values 
#   

# arguments ___________________________________________________________________________
#   : dataSet.mvs       - expression matrix containing abundances with MVs 
#                         (either peptides or proteins)
#   : model.selector    - binary vector; "1" indicates MAR/MCAR proteins
#   : method            - the method to be used for MAR/MCAR missing values
#                       - possible values: MLE (default), SVD, KNN

# output ______________________________________________________________________________
#   : dataSet.imputed   - dataset containing only MNAR (assumed to be left-censored)
#                         missing data

impute.MAR = function(dataSet.mvs, model.selector, method = "MLE"){
  
  if (length(which(model.selector[[1]]==1)) == 0){
    dataSet.imputed = dataSet.mvs
  }
  else{
    # ___________________________________________________________________________________
    # select MCAR proteins 
    # -----------------------------------------------------------------------------------
    dataSet.MCAR = dataSet.mvs[which(model.selector[[1]]==1),]
    
    # ___________________________________________________________________________________
    # perform imputation using the specified method
    # -----------------------------------------------------------------------------------
    
    switch(method,
           MLE = {
             dataSet.MCAR.imputed = impute.wrapper.MLE(dataSet.MCAR) 
           },
           
           SVD = {
             dataSet.MCAR.imputed = impute.wrapper.SVD(dataSet.MCAR, K = 2) 
           },
           
           KNN = {
             dataSet.MCAR.imputed = impute.wrapper.KNN(dataSet.MCAR, K = 15) 
           }     
    )
    
    # ___________________________________________________________________________________
    # replace imputed MCAR in the data matrix
    # -----------------------------------------------------------------------------------
    dataSet.imputed = dataSet.mvs
    dataSet.imputed[which(model.selector[[1]]==1),] = dataSet.MCAR.imputed
    
  }
  
  return(dataSet.imputed)
  
}
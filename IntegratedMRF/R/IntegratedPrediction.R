#' Integrated Prediction of Testing samples from integrated RF or MRF model 
#' 
#' Generates Random Forest or Multivariate Random Forest model for each subtype of dataset and predicts testing samples using the generated models. 
#' Subsequently, the prediction for different subtypes of dataset are combined using the Combination weights generated from 
#' Integrated Model which is based on Bootstrap error estimate
#'  
#' @param finalX List of Matrices where each matrix represent a specific data subtype (such as genomic characterizations for 
#' drug sensitivity prediction). Each subtype can have different types of features. For example, if there are three subtypes containing
#'  100, 200 and 250 features respectively,  finalX will be a list containing 3 matrices of sizes M x 100, M x 200 and M x 250 
#'  where M is the number of Samples.
#' @param finalY_train A M x T matrix of output features for training samples, where M is number of samples and
#' T is the number of output features. The dataset is assumed to contain no missing values. If there are missing values, an imputation method 
#' should be applied before using the function. A function 'Imputation' is included within the package.
#' @param Cell It contains a list of samples (the samples can be represented either numerically by indices or by names) for each data subtype. 
#' For the example of 3 data subtypes, it will be a list containing 3 arrays where each array contains the sample information for each data subtype.
#' @param finalY_train_cell Cell lines of output features for training samples
#' @param finalY_test_cell Cell lines of output features for testing samples
#' @param n_tree number of trees in the forest
#' @param m_feature Number of randomly selected features considered for a split in each regression tree node.
#' @param min_leaf minimum number of samples in the leaf node 
#' @param Serial Consists of a  list of all combinations of different subtypes of a dataset (except for the case with no dataset being selected). 
#' For example, if a 
#' dataset has 3 subtypes, then Serial is a list of size 2^3-1=7.  The ordering of the seven sets will be [1 2 3], [1 2], [1 3], [2 3], [1], [2], [3]
#' 
#' @return Final Prediction of testing samples based on provided testing sample names.
#' @details
#' Input matrix and output response of training samples have been used to build Random Forest or Multivariate Random Forest model for each subtype of 
#' a dataset. These models are used to calculate prediction of 
#' testing samples for each subtype separately. Subsequently Combination Weights are used to integrate the predictions from data subtypes. 
#' 
#' Combination Weight Generation: For M x N dataset, N number of bootstrap sampling sets are considered. For each bootstrap sampling set and each subtype, a Random Forest (RF) 
#' or, Multivariate Random Forest (MRF) model is generated, which is used for calculating the prediction performance for out-of-bag samples.  
#' The prediction performance for each subtype of the dataset is based on the averaging over different bootstrap training sets. 
#' The combination weights (regression coefficients) for each combination of subtypes are generated using least Square Regression from the 
#' individual subtype predictions and used integrate the predictions from data subtypes.
#'  
#' The specific set of combination weights to be used for testing samples will depend on the number of data subtypes available 
#' for the testing samples. Note that not all subtype information maybe available for all samples. 
#' As an example with three data subtypes, a testing sample with all subtype data available will using 
#' the combination weights corresponding to Serial [1 2 3] where if subtype 3 is not available, the function will 
#' using the combination weights corresponding to Serial [1 2].
#' @examples
#' #library(IntegratedPredictionUsingRandomForest)
#' #n_tree=10
#' #m_feature=5
#' #min_leaf=3
#' #Cell=NULL
#' #Expression=NULL
#' #finalX=NULL
#' #library(openxlsx)
#' #for (i in 1:5){#5=number_of_subtypes_in_dataset
#' #     Genome=read.xlsx("Subtype_filename.xlsx")
#' #     Expression[[i]]=Genome[complete.cases(Genome),]#remove all feature vector with NaN
#' #     Cell[[i]]=colnames(Expression[[i]], do.NULL = TRUE, prefix = "col")[-1]
#' #     # Taking the cell line names for that subtype of dataset
#' #     finalX[[i]]=matrix(as.numeric(t(Expression[[i]])[-1,]),nrow=length(Cell[[i]]))
#' #     #Input Matrix(MxN), with M number of samples and N number of features
#' #}
#' #Drug_Sen_train <- read.xlsx("Output_Response_FileName.xlsx", colNames = TRUE)
#' #for (j in 1:3){#3=Number_of_output_Response
#' #     XX=matrix(Drug_Sen_train[,Column_of_the_Response_for_Prediction],ncol=1)
#' #     finalY_train[,j]=matrix(Imputation(XX),ncol=1)
#' #}
#' #finalY_train_cell=Drug_Sen_train[,1]
#' #Serial=NULL
#' #library(caTools)
#' #for (p in length(Cell):1){
#' #       nk=combs(1:5,p)
#' #       sk=length(Serial)
#' #       for (q in 1:dim(nk)[1]){
#' #               Serial[[sk+q]]=nk[q, ]
#' #       }
#' #}
#' #Drug_Sen_test <- read.xlsx("Output_Response_Test.xlsx", sheet = 1,  startRow = 1, colNames = TRUE)
#' #finalY_test_cell=Drug_Sen_test[,1]
#' #Prediction=IntegratedPrediction(finalX,finalY_train,Cell,finalY_train_cell,finalY_test_cell,
#' #+n_tree,m_feature,min_leaf,Serial)
#'
#' @export

IntegratedPrediction <- function(finalX,finalY_train,Cell,finalY_train_cell,finalY_test_cell,n_tree,m_feature,min_leaf,Serial){
  ##
  Common_cell_train=finalY_train_cell
  for (q in 1:length(Cell)){
    Common_cell_train=intersect(Common_cell_train,Cell[[q]])
  }
  
  final=NULL
  for (q in 1:length(Cell)){
    Cell_ind=match(Common_cell_train,Cell[[q]])
    final[[q]]=finalX[[q]][Cell_ind, ]
  }
  finalY=NULL
  ia6=match(Common_cell_train,finalY_train_cell)
  finalY=matrix(finalY_train[ia6,],ncol=ncol(finalY_train))
  Variable_number=ncol(finalY)
  if (Variable_number>1){
    Command=2
  }else if(Variable_number==1){
    Command=1
  } 
  ################################## BSP ###############################
  if (nrow(finalY)<50){
    N=floor(0.75*nrow(finalY))
  }else if (nrow(finalY)>=50 && nrow(finalY)<101){
    N=floor(nrow(finalY)/2)
  }else if (nrow(finalY)>=101){
    N=floor(nrow(finalY)/5)
  }
  
  Y_hat_BSP=NULL
  for (q in 1:length(Cell)){
    Y_hat_BSP[[q]]=matrix(rep(0,length(finalY)),ncol=ncol(finalY))
  }
  bootsam_FF=NULL
  Index=NULL
  Index=1:nrow(finalY)
  #   library(bootstrap) 
  theta <- function(x){x}
  results <- bootstrap::bootstrap(Index,N,theta) #no indics, gives number
  bootsam=results$thetastar
  
  Store=rep( list(NULL), length(Cell) )
  for (q in 1:length(Cell)){
    Store[[q]]=rep( list(NULL), nrow(finalY) )
  }
  
  for (FF in 1:N){
    bootsam_FF=bootsam[,FF]
    Index_FF=unique(bootsam_FF)
    Index_pred=setdiff(Index, Index_FF)
    
    finalY_bsp=matrix(finalY[bootsam_FF,],ncol=Variable_number)
    finalY_bsp_pred=matrix(finalY[Index_pred,],ncol=Variable_number)
    
    final_genome=NULL
    finalX_bsp=NULL
    finalX_bsp_pred=NULL
    finalY_pred=NULL
    
    for (q in 1:length(Cell)){
      finalX_bsp[[q]]=final[[q]][bootsam_FF,]
      finalX_bsp_pred[[q]]=final[[q]][Index_pred,]
      finalY_pred[[q]]=build_forest_predict(finalX_bsp[[q]], finalY_bsp, n_tree, m_feature, min_leaf, finalX_bsp_pred[[q]])
      #final_genome=cbind(final_genome,finalY_pred[[q]])
      for (R in 1:length(Index_pred)){
        Store[[q]][[Index_pred[R]]]=rbind(Store[[q]][[Index_pred[R]]],finalY_pred[[q]][R,])
      }
    }
  }
  
  final_genome_BSP=NULL
  for (RR in 1:Variable_number){
    final_genome_BSP[[RR]]=matrix(rep(0,length(Cell)*nrow(finalY)),nrow=nrow(finalY))
  }
  for (RR in 1:Variable_number){
    for (q in 1:length(Cell)){
      for (FF in 1:nrow(finalY)){
        final_genome_BSP[[RR]][FF,q]=mean(Store[[q]][[FF]][,RR])
      }
    }
  }
  ## Find Combination Weight
  BSP_coeff=rep(list(NULL), length(Serial))
  for (S in 1:length(Serial)){
    BSP_coeff[[S]]=matrix(rep(0,length(Serial[[S]])*Variable_number),ncol=Variable_number)
  }
  
  for (RR in 1:Variable_number){
    for (q in 1:length(Cell)){
      Y_hat_BSP[[q]]=final_genome_BSP[[RR]][,q]
    }
    for (S in 1:length(Serial)){
      final_genome_BSP1=NULL
      W=Serial[[S]]
      for (q in 1:length(Serial[[S]])){
        final_genome_BSP1=cbind(final_genome_BSP1,matrix(Y_hat_BSP[[W[q]]],ncol=1))
      }
      BSP_coeff[[S]][,RR]=matrix(limSolve::lsei(A=final_genome_BSP1, B=finalY[,RR], E=rep(1,dim(final_genome_BSP1)[2]), F=1)$X, ncol=1)
    }
  }
  ############################
  Common_cell_dataset=NULL
  for (q in 1:length(Cell)){
    Common_cell_dataset[[q]]=intersect(finalY_train_cell,Cell[[q]])
  }
  
  final_dataset=NULL
  finalY_dataset=NULL
  for (q in 1:length(Cell)){
    Cell_ind=match(Common_cell_dataset[[q]],Cell[[q]])
    final_dataset[[q]]=finalX[[q]][Cell_ind, ]
    final_dataset[[q]]=matrix(as.numeric(final_dataset[[q]]),nrow = dim(final_dataset[[q]])[1], ncol = dim(final_dataset[[q]])[2])
    
    Cell_ind_Y=match(Common_cell_dataset[[q]],finalY_train_cell)
    finalY_dataset[[q]]=matrix(finalY_train[Cell_ind_Y,],ncol=Variable_number)
  }
  
  Y_hat_test=NULL
  final_test=NULL
  Common_cell_test=NULL
  for (q in 1:length(Cell)){
    Common_cell_test[[q]]=intersect(finalY_test_cell,Cell[[q]])
    Cell_ind=match(Common_cell_test[[q]],Cell[[q]])
    final_test[[q]]=finalX[[q]][Cell_ind, ]
    final_test[[q]]=matrix(as.numeric(final_test[[q]]),nrow = dim(final_test[[q]])[1], ncol = dim(final_test[[q]])[2])
    
    Y_hat_test[[q]] = build_forest_predict(final_dataset[[q]],finalY_dataset[[q]], n_tree, m_feature, min_leaf, final_test[[q]])
  }
  Common_cell=NULL
  Final_test=NULL
  Drug_sensitivity_cell_test2=finalY_test_cell
  for (q in 1:length(Serial)){
    D=NULL
    D=Drug_sensitivity_cell_test2
    for (check in 1:length(Serial[[q]])){
      D=intersect(D,Cell[Serial[[q]][check]][[1]])
    }
    Common_cell[[q]]=D
    Drug_sensitivity_cell_test2=setdiff(Drug_sensitivity_cell_test2,D)
    
    F_test=matrix(rep(0,length(Common_cell[[q]])*Variable_number),ncol=Variable_number)
    match_ind=NULL
    W=Serial[[q]]
    for (RR in 1:Variable_number){
      for (R in 1:length(W)){
        match_ind[[R]]=match(Common_cell[[q]],Common_cell_test[[W[R]]])
        F_test[,RR]=F_test[,RR]+BSP_coeff[[q]][R,RR]*matrix(Y_hat_test[[W[R]]][match_ind[[R]],RR],ncol=1)
      }
    }
    
    Final_test[[q]]=F_test
  }
  Final_result=NULL
  for (q in 1:length(Serial)){
    if (length(Common_cell[[q]])>0){
      Final_result=rbind(Final_result,cbind(matrix(Common_cell[[q]],ncol=1),matrix(Final_test[[q]],ncol=Variable_number)))
    }
  }
  match_ind=match(finalY_test_cell,Final_result[,1])
  Final_prediction=Final_result[match_ind,]
  return(Final_prediction)
}
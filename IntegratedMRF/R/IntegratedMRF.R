#' Weights for combination of predictions from different data subtypes using Least Square Regression based on various error estimation techniques 
#' 
#' Calculates combination weights for different subtypes of dataset combinations to generate integrated Random Forest (RF) or Multivariate 
#' Random Forest (MRF) model based on different error estimate options 
#' of Bootstrap, Re-substitution, 0.632 Bootstrap or Leave one out.   
#'  
#' @param finalX List of Matrices where each matrix represent a specific data subtype (such as genomic characterizations for 
#' drug sensitivity prediction). Each subtype can have different types of features. For example, if there are three subtypes containing
#'  100, 200 and 250 features respectively,  finalX will be a list containing 3 matrices of sizes M x 100, M x 200 and M x 250 
#'  where M is the number of Samples. 
#' @param finalY_train A M x T matrix of output features for training samples, where M is number of samples and T is the number of output features. 
#' The dataset is assumed to contain no missing values. If there are missing values, an imputation method should be applied before using the function. 
#' A function 'Imputation' is included within the package.
#' @param Cell It contains a list of samples (the samples can be represented either numerically by indices or by names) for each data subtype. 
#' For the example of 3 data subtypes, it will be a list containing 3 arrays where each array contains the sample information for each data subtype.
#' @param finalY_train_cell Sample names of output features for training samples
#' @param n_tree Number of trees in the forest
#' @param m_feature Number of randomly selected features considered for a split in each regression tree node.
#' @param min_leaf Minimum number of samples in the leaf node 
#' @param Serial Consists of a  list of all combinations of different subtypes of a dataset (except for the case with no dataset being selected). 
#' For example, if a 
#' dataset has 3 subtypes, then Serial is a list of size 2^3-1=7.  The ordering of the seven sets will be [1 2 3], [1 2], [1 3], [2 3], [1], [2], [3]
#' @param Confidence_Level Confidence level for calculation of confidence interval (User Defined)
#' @return 
#' List with the following components: 
#' \item{BSP_coeff}{Combination weights using Bootstrap Error Model, where index are in list format. 
#' If number of genomic characterizations or subtypes of dataset is 5, then there will be 2^5-1=31 list of weights}
#' \item{Resub_coeff}{Combination weights using Resubstituition Error Model, where index are in list format. 
#' If number of genomic characterizations or subtypes of dataset is 5, then there will be 2^5-1=31 list of weights}
#' \item{BSP632_coeff}{Combination weights using 0.632Bootstrap Error Model, where index are in list format. 
#' If number of genomic characterizations or subtypes of dataset is 5, then there will be 2^5-1=31 list of weights}
#' \item{LOO_coeff}{Combination weights using Leave-One-Out Error Model, where index are in list format. 
#' If number of genomic characterizations or subtypes of dataset is 5, then there will be 2^5-1=31 list of weights} 
#' \item{Error}{Matrix of Mean Absolute Error, Mean Square Error and correlation between actual and predicted responses for integrated model based 
#' on Bootstrap, Re-substitution, 0.632Bootstrap and Leave-one-out error estimation sampling techniques for the integrated model 
#' containing all the data subtypes}
#' \item{Confidence Interval}{Low and High confidence interval for a user defined confidence level for the drug using Jackknife-After-Bootstrap Approach in a list}
#' \item{BSP_error_all_mae}{Bootstrap Mean Absolute Errors (MAE) for all combinations of the dataset subtypes. Size C x R, where C is the number of 
#' combinations and R is the number of output responses. C is in decreasing order, which means first value is combination of all subtypes 
#' and next ones are in decreasing order. For example, if a dataset has 3 subtypes, then C is equal to 2^3-1=7.  The ordering of C is the combination of
#' subtypes [1 2 3], [1 2], [1 3], [2 3], [1], [2], [3] }
#' \item{Resub_error_all_mae}{Re-substituition Mean Absolute Errors (MAE) for all combinations of the dataset subtypes. Size C x R, where C is the number of 
#' combinations and R is the number of output responses. C is in decreasing order, which means first value is combination of all subtypes 
#' and next ones are in decreasing order. For example, if a dataset has 3 subtypes, then C is equal to 2^3-1=7.  The ordering of C is the combination of
#' subtypes [1 2 3], [1 2], [1 3], [2 3], [1], [2], [3] }
#' \item{BSP632_error_all_mae}{0.632Bootstrap Mean Absolute Errors (MAE) for all combinations of the dataset subtypes. Size C x R, where C is the number of 
#' combinations and R is the number of output responses. C is in decreasing order, which means first value is combination of all subtypes 
#' and next ones are in decreasing order. For example, if a dataset has 3 subtypes, then C is equal to 2^3-1=7.  The ordering of C is the combination of
#' subtypes [1 2 3], [1 2], [1 3], [2 3], [1], [2], [3] }
#' \item{LOO_error_all_mae}{Leave One Out Mean Absolute Errors (MAE) for all combinations of the dataset subtypes. Size C x R, where C is the number of 
#' combinations and R is the number of output responses. C is in decreasing order, which means first value is combination of all subtypes 
#' and next ones are in decreasing order. For example, if a dataset has 3 subtypes, then C is equal to 2^3-1=7.  The ordering of C is the combination of
#' subtypes [1 2 3], [1 2], [1 3], [2 3], [1], [2], [3] }
#' The function also returns figures of different error estimation in .tiff format 
#' @details
#' The function takes all the subtypes of dataset in matrix format and its corresponding sample information.
#' For the calculation purpose, we have taken the data of the samples that are common in all the subtypes and output training responses.
#' For example,  let a dataset has 3 sub-types with different number of samples and features, while indices of samples in subtype 1, 2, 3  and output feature matrix
#' is 1:10, 3:15, 5:16 and 5:11 respectively. So, features of sample index 5:10 (common to all subtypes and output feature matrix) of all subtypes and output feature 
#' matrix will be separated and considered for all calculations. 
#' 
#' For M x N dataset, N number of bootstrap sampling sets are considered. For each bootstrap sampling set and each subtype, a Random Forest (RF) 
#' or, Multivariate Random Forest (MRF) model is generated, which is used for calculating the prediction performance for out-of-bag samples.  
#' The prediction performance for each subtype of the dataset is based on the averaging over different bootstrap training sets. 
#' The combination weights (regression coefficients) for each combination of subtypes are generated using least Square Regression from the 
#' individual subtype predictions and used later to calculate mean absolute error, mean square error and correlation coefficient between 
#' predicted and actual values.
#' 
#' For re-substitution error estimation with M cell lines, 1 model is generated for each subtype of dataset, which is then used to 
#' calculate errors and combination weights for different data subtype combinations.   
#' 
#' For 0.632 Bootstrap error estimation, prediction of bootstrap and re-substitution error estimation is combined using 
#' 0.632xBootstrap Error + 0.368xRe-substitution Error. 
#' These prediction results are then used to compute the errors and combination weights for different data subtype combinations.
#' 
#' Confidence Interval has been calculated using Jackkniffe-After-Bootstrap Approach and prediction result of bootstrap error estimation.
#' 
#' For leave-one-out error estimation using M cell lines, M models are generated for each subtype of dataset, which are then used to 
#' calculate the errors and combination weights for different data subtype combinations.
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
#' ## Combination Index using Different Error Estimation Method
#' #Result=Combination(finalX,finalY_train,Cell,finalY_train_cell,n_tree,m_feature,min_leaf,Serial)
#' @importFrom
#' grDevices dev.off tiff
#' @references
#' Wan, Qian, and Ranadip Pal. "An ensemble based top performing approach for 
#' NCI-DREAM drug sensitivity prediction challenge." PloS one 9.6 (2014): e101183.
#' @export

Combination <- function(finalX,finalY_train,Cell,finalY_train_cell,n_tree,m_feature,min_leaf,Serial,Confidence_Level){
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
  ptm1=proc.time()
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
  Store_Jack=rep( list(NULL), Variable_number )
  for (RR in 1:Variable_number){
    Store_Jack[[RR]]=rep( list(NULL), length(Serial))
    for (q in 1:length(Serial)){
      Store_Jack[[RR]][[q]]=rep( list(NULL), nrow(finalY) )
    }
  }
  BSP_error_alll_mae=rep(list(NULL), Variable_number)
  BSP_error_alll_mse=rep(list(NULL), Variable_number)
  BSP_error_alll_corr=rep(list(NULL), Variable_number)
  for (S in 1:Variable_number){
    BSP_error_alll_mae[[S]]=matrix(rep(0,length(Serial)*N),ncol=N)
    BSP_error_alll_mse[[S]]=matrix(rep(0,length(Serial)*N),ncol=N)
    BSP_error_alll_corr[[S]]=matrix(rep(0,length(Serial)*N),ncol=N)
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
    
    for (RR in 1:Variable_number){
      for (S in 1:length(Serial)){
        final_genome_BSP1=NULL
        W=Serial[[S]]
        for (q in 1:length(Serial[[S]])){
          final_genome_BSP1=cbind(final_genome_BSP1,matrix(finalY_pred[[W[q]]][,RR],ncol=1))
        }
        BSP_temp_coeff1=matrix(limSolve::lsei(A=final_genome_BSP1, B=matrix(finalY_bsp_pred[,RR],ncol=1), E=rep(1,dim(final_genome_BSP1)[2]), F=1)$X, ncol=1)
        BSP_error_alll_mae[[RR]][S,FF]=mean(abs(final_genome_BSP1%*%BSP_temp_coeff1-finalY_bsp_pred[,RR]))
        BSP_error_alll_mse[[RR]][S,FF]=mean((final_genome_BSP1%*%BSP_temp_coeff1-finalY_bsp_pred[,RR])^2)
        BSP_error_alll_corr[[RR]][S,FF]=stats::cor(final_genome_BSP1%*%BSP_temp_coeff1,finalY_bsp_pred[,RR])
        Err_bsp_temp=abs(final_genome_BSP1%*%BSP_temp_coeff1-finalY_bsp_pred[,RR])
        for (R in 1:length(Index_pred)){
          Store_Jack[[RR]][[S]][[Index_pred[R]]]=rbind(Store_Jack[[RR]][[S]][[Index_pred[R]]],Err_bsp_temp[R,])
        }
      }
    }
  }
  Error_BSP_Jack=rep( list(NULL), Variable_number )
  for (RR in 1:Variable_number){
    Error_BSP_Jack[[RR]]=matrix(rep(0,length(Serial)*nrow(finalY)),ncol=nrow(finalY))
    for (S in 1:length(Serial)){
      for (FFF in 1:nrow(finalY)){
        Error_BSP_Jack[[RR]][S,FFF]=mean(Store_Jack[[RR]][[S]][[FFF]])
      }
    }
  }
  BSP_error_all_mae=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  BSP_error_all_mse=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  BSP_error_all_corr=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  BSP_errors=NULL
  for (RR in 1:Variable_number){
    for (S in 1:length(Serial)){
      BSP_error_all_mae[S,RR]=mean(BSP_error_alll_mae[[RR]][S,])
      BSP_error_all_mse[S,RR]=mean(BSP_error_alll_mse[[RR]][S,])
      BSP_error_all_corr[S,RR]=mean(BSP_error_alll_corr[[RR]][S,])
    }
    BSP_errors=cbind( BSP_errors,rbind(BSP_error_all_mae[1,RR],BSP_error_all_mse[1,RR],BSP_error_all_corr[1,RR]))
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
  final_BSP=matrix(rep(0,length(finalY)),ncol=Variable_number)
  for (RR in 1:Variable_number){
    Store2=rep( list(NULL), nrow(finalY) )
    for (FF in 1:N){
      bootsam_FF=bootsam[,FF]
      Index_FF=unique(bootsam_FF)
      Index_pred=setdiff(Index, Index_FF)
      
      finalY_bsp1=matrix(finalY[bootsam_FF,RR],ncol=1)
      finalY_bsp1_pred=matrix(finalY[Index_pred,RR],ncol=1)
      
      BSP_temp_coeff=matrix(limSolve::lsei(A=final_genome_BSP[[RR]][bootsam_FF,], B=finalY_bsp1, E=rep(1,dim(final_genome_BSP[[RR]])[2]), F=1)$X, ncol=1)
      final_BSP_index=final_genome_BSP[[RR]][Index_pred,]%*%BSP_temp_coeff
      
      for (R in 1:length(Index_pred)){
        Store2[[Index_pred[R]]]=c(Store2[[Index_pred[R]]],final_BSP_index[R,])
      }
    }
    for (FF in 1:nrow(finalY)){
      final_BSP[FF,RR]=mean(Store2[[FF]])
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
  ptm2=proc.time()-ptm1
  message("Elapsed Time for Bootstrap Error Estimation is ", ptm2[[3]])
  ####################### Resubstitution Error######################
  ptm1 <- proc.time()
  finalY_pred_resub=NULL
  final_genome_resub=rep(list(NULL), Variable_number)
  for (q in 1:length(Cell)){
    finalY_pred_resub[[q]]=build_forest_predict(final[[q]], finalY, n_tree, m_feature, min_leaf, final[[q]])
    for (RR in 1:Variable_number){
      final_genome_resub[[RR]]=cbind(final_genome_resub[[RR]],finalY_pred_resub[[q]][,RR])
    }
  }
  Resub_errors=NULL
  Resub_error_mae=NULL
  Resub_error_mse=NULL
  Resub_corr=NULL
  final_resub=matrix(rep(0,length(finalY)),ncol=Variable_number)
  for (RR in 1:Variable_number){
    error=error_calculation(final_genome_resub[[RR]],matrix(finalY[,RR],ncol=1))
    Resub_error_mae[RR]=error[[2]]
    Resub_error_mse[RR]=error[[3]]
    Resub_corr[RR]=error[[4]]
    final_resub[,RR]=error[[1]]
  }
  Resub_errors=rbind(Resub_error_mae,Resub_error_mse,Resub_corr)
  
  # Resub_coeff=matrix(lsei(A=final_genome_resub, B=finalY, E=rep(1,dim(final_genome_resub)[2]), F=1)$X, ncol=1)
  Resub_coeff=rep(list(NULL), length(Serial))
  Y_hat_resub=NULL
  for (S in 1:length(Serial)){
    Resub_coeff[[S]]=matrix(rep(0,length(Serial[[S]])*Variable_number),ncol=Variable_number)
  }
  Error_resub_Jack=rep( list(NULL), Variable_number )
  Resub_error_all_mae=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  for (RR in 1:Variable_number){
    Error_resub_Jack[[RR]]=matrix(rep(0,length(Serial)*nrow(finalY)),ncol=nrow(finalY))
    for (q in 1:length(Cell)){
      Y_hat_resub[[q]]=final_genome_resub[[RR]][,q]
    }
    for (S in 1:length(Serial)){
      final_genome_Resub1=NULL
      W=Serial[[S]]
      for (q in 1:length(Serial[[S]])){
        final_genome_Resub1=cbind(final_genome_Resub1,matrix(Y_hat_resub[[W[q]]],ncol=1))
      }
      Resub_coeff[[S]][,RR]=matrix(limSolve::lsei(A=final_genome_Resub1, B=finalY[,RR], E=rep(1,dim(final_genome_Resub1)[2]), F=1)$X, ncol=1)
      Resub_error_all_mae[S,RR]=mean(abs(final_genome_Resub1%*%Resub_coeff[[S]][,RR]-finalY[,RR]))
      Error_resub_Jack[[RR]][S,]=abs(final_genome_Resub1%*%Resub_coeff[[S]][,RR]-finalY[,RR])
    }
  }
  ptm2=proc.time()-ptm1
  message("Elapsed Time for Resubstitution Error Estimation is ", ptm2[[3]])
  ####################### 0.632BSP Error######################
  ptm1 <- proc.time()
  BSP632_error_mae=0.632*BSP_error_all_mae[1,]+0.368*Resub_error_mae
  BSP632_error_mse=0.632*BSP_error_all_mse[1,]+0.368*Resub_error_mse
  BSP632_corr=0.632*BSP_error_all_corr[1,]+0.368*Resub_corr
  BSP632_errors=rbind(BSP632_error_mae,BSP632_error_mse,BSP632_corr)
  
  BSP632_coeff=rep(list(NULL), length(Serial))
  Y_hat_BSP632=NULL
  for (S in 1:length(Serial)){
    BSP632_coeff[[S]]=matrix(rep(0,length(Serial[[S]])*Variable_number),ncol=Variable_number)
  }
  BSP632_error_all_mae=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  for (RR in 1:Variable_number){
    for (q in 1:length(Cell)){
      Y_hat_BSP632[[q]]=0.632*final_genome_BSP[[RR]][,q]+0.368*final_genome_resub[[RR]][,q]
    }
    for (S in 1:length(Serial)){
      final_genome_BSP6321=NULL
      W=Serial[[S]]
      for (q in 1:length(Serial[[S]])){
        final_genome_BSP6321=cbind(final_genome_BSP6321,matrix(Y_hat_BSP632[[W[q]]],ncol=1))
      }
      BSP632_coeff[[S]][,RR]=matrix(limSolve::lsei(A=final_genome_BSP6321, B=finalY[,RR], E=rep(1,dim(final_genome_BSP6321)[2]), F=1)$X, ncol=1)
      BSP632_error_all_mae[S,RR]=0.632*BSP_error_all_mae[S,RR]+0.368*Resub_error_all_mae[S,RR]#mean(abs(final_genome_BSP6321%*%BSP632_coeff[[S]][,RR]-finalY[,RR]))
    }
  }
  ptm2=proc.time()-ptm1
  message("Elapsed Time for 0.632 Bootstrap Error Estimation is ", ptm2[[3]])
  ################################ 0.632BSP error Jackkniffe Confidence Interval ####################
  ptm1 <- proc.time()
  Low_confidence=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  High_confidence=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  mean_D=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  for (RR in 1:Variable_number){  
    for (S in 1:length(Serial)){
      err=0.632*Error_BSP_Jack[[RR]][S,]+0.368*Error_resub_Jack[[RR]][S,]
      #err=abs((0.632*final_BSP[,RR]+0.368*final_resub[,RR])-finalY[,RR])
      n=length(err)
      s=sqrt(((n-1)/n)*sum((err-mean(err))^2))
      alpha=1-(Confidence_Level/100)
      # z_alpha=sqrt(2)*erfinv(2*alpha/2-1);
      z_alpha=stats::qnorm((1+2*alpha/2-1)/2)
      # z_1alpha=sqrt(2)*erfinv(2*(1-alpha/2)-1);
      z_1alpha=stats::qnorm((1+2*(1-alpha/2)-1)/2)
      mean_D[S,RR]=mean(err)
      
      Low_confidence[S,RR]=max(0,(mean(err)-(s*abs(z_alpha))))# Z_a/2 error function or, quantile function
      High_confidence[S,RR]=mean(err)+(s*(z_1alpha))
    }
  }
  ptm2=proc.time()-ptm1
  message("Elapsed Time for Jackknife Confidence Interval Calculation is ", ptm2[[3]])
  ################################ Leave One out error #######################################
  ptm1 <- proc.time()
  #   if (nrow(finalY)<50){
  #     N_LOO=floor(0.75*nrow(finalY))
  #   }else if (nrow(finalY)>=50 && nrow(finalY)<101){
  #     N_LOO=floor(nrow(finalY)/2)
  #   }else if (nrow(finalY)>=101){
  #     N_LOO=floor(nrow(finalY)/5)
  #   }
  Y_hat_LOO=NULL
  for (q in 1:length(Cell)){
    Y_hat_LOO[[q]]=matrix(rep(0,length(finalY)),ncol=Variable_number)
  }
  
  for (q in 1:length(Cell)){
    for (FF in 1:nrow(finalY)){
      Index2=1:nrow(finalY)
      Index2=setdiff(Index2,FF)
      Y1=matrix(finalY[Index2,],ncol=Variable_number)
      ##
      X1=final[[q]][Index2,]
      Xt=matrix(final[[q]][FF,],nrow=1)
      Y_hat_LOO[[q]][FF,] = build_forest_predict(X1,Y1, n_tree, m_feature, min_leaf, Xt)
    }
  }
  LOO_errors=NULL
  LOO_error_mae=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  LOO_error_mse=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  LOO_corr=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  
  LOO_coeff=rep(list(NULL), length(Serial))
  for (S in 1:length(Serial)){
    LOO_coeff[[S]]=matrix(rep(0,length(Serial[[S]])*Variable_number),ncol=Variable_number)
  }
  for (RR in 1:Variable_number){
    for (S in 1:length(Serial)){
      final_genome_LOO=NULL
      W=Serial[[S]]
      for (q in 1:length(Serial[[S]])){
        final_genome_LOO=cbind(final_genome_LOO,matrix(Y_hat_LOO[[W[q]]][,RR],ncol=1))
      }
      LOO_coeff[[S]][,RR]=matrix(limSolve::lsei(A=final_genome_LOO, B=finalY[,RR], E=rep(1,dim(final_genome_LOO)[2]), F=1)$X, ncol=1)
      LOO_error_mae[S,RR]=mean(abs(final_genome_LOO%*%LOO_coeff[[S]][,RR]-finalY[,RR]))
      LOO_error_mse[S,RR]=mean((final_genome_LOO%*%LOO_coeff[[S]][,RR]-finalY[,RR])^2)
      LOO_corr[S,RR]=stats::cor(final_genome_LOO%*%LOO_coeff[[S]][,RR],finalY[,RR])
    }
  }
  LOO_errors=rbind(LOO_error_mae[1,],LOO_error_mse[1,],LOO_corr[1,])
  ptm2=proc.time()-ptm1
  message("Elapsed Time for Leave-one-out Error Estimation is ", ptm2[[3]])
  
  Result=NULL
  Result[[1]]=BSP_coeff
  Result[[2]]=Resub_coeff
  Result[[3]]=BSP632_coeff
  Result[[4]]=LOO_coeff
  Result[[5]]=rbind(BSP_errors,Resub_errors,BSP632_errors,LOO_errors)
  Result[[6]]=rep( list(NULL), 2 )
  Result[[6]][[1]]=Low_confidence
  Result[[6]][[2]]=High_confidence
  Result[[7]]=BSP_error_all_mae
  Result[[8]]=Resub_error_all_mae
  Result[[9]]=BSP632_error_all_mae
  Result[[10]]=LOO_error_mae
  #### Figures ######
  Error_Matrix1=rbind(Result[[7]],Result[[8]],Result[[9]],Result[[10]])
  Error_Matrix2=matrix(c(rep("BSP_error",length(Serial)),rep("Resub_error",length(Serial)),rep("0.632 BSP_error",length(Serial)),rep("LOO_error",length(Serial))),ncol=1)
  
  SS=matrix(rep(0,length(Serial)*length(Cell)),ncol=length(Cell))
  for (i in 1:length(Serial)){
    for (j in 1:length(Cell)){
      if (length(which(j==Serial[[i]]))==1){
        SS[i,j]=1
      }
    }
  }
  SSC=transform(paste0(SS[,1],SS[,2]))
  if (length(Cell)>2) {
    for (k in 3:length(Cell)){
      SSC=transform(paste0(SSC[,1],SS[,k]))
    }
  }
  Error_type=NULL
  Error_Matrix3=cbind(rep(matrix(SSC[,1],ncol=1),4))
  Error_Matrix=cbind(Error_Matrix2,Error_Matrix3,Error_Matrix1)
  colnames(Error_Matrix) <- c("Error_type","Combination",rep("Response",ncol(finalY_train)))
  for (i in 1:ncol(finalY_train)){
    Error_Matrix[,i+2]=matrix(format(round(as.numeric(Error_Matrix[,i+2]), 4), nsmall = 4),ncol=1)
  }
  #  library(ggplot2)
  plot_list = list()
  for (i in 1:ncol(finalY_train)){
    qq=ggplot2::ggplot(data.frame(Error_Matrix), ggplot2::aes(x=Combination,y=Error_Matrix[,(i+2)], fill=Error_type))+ ggplot2::geom_bar( stat="identity",position="dodge")+ ggplot2::labs(title = "Different Error Estimation for Integrated Models")+ ggplot2::xlab("Combination Index of Integrated Models")+ggplot2::ylab("Mean Absolute Error")+ ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust=0.5, hjust = 0))
    plot_list[[i]] = qq  
  }
  for (i in 1:ncol(finalY_train)){
    tiff(paste("Error_Estimation_Response_",i,".tiff",sep=""), width = 7, height = 10, units = 'in', res = 700)
    print(plot_list[[i]])
    dev.off()
  }
  
  return(Result)
}
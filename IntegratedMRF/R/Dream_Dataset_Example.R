#' NCI-Dream Drug Sensitivity Prediction Challenge Dataset Implementation using IntegratedPredictionUsingRandomForest Package
#' @details
#' This is a practical implementation of Combination and CombPredict function. The dataset has been used in NCI-Dream Drug Sensitivity 
#' Prediction Challenge. 2 genomic characterization Gene Expression and Methylation is considered here to show the Integrated prediction.
#' If the number of trees in the forest is 20, then it takes around 13-15 minutes to run the whole code.
#' @examples
#' ## Download the data from "https://www.synapse.org/#!Synapse:syn2785783"
#' ## Uncomment all the lines with # to run the code. 
#' ## Join the line with "+" with previous line
#' ## The line with ## have been used to comment on the code. 
#' # library(IntegratedPredictionUsingRandomForest)
#' # Response=c(2,3,10,15) 
#' ## Output Response Number 2, 3, 10 and 15 of the dataset used to run the package, 
#' ## user can specify the response number, he wants to work with.
#' # n_tree=20  ## Number of trees in the forest
#' # m_feature=5 ## Number of random features for a split in each regression tree node
#' # min_leaf=1 ## Minimum number of samples in the leaf node
#' # Confidence_Level=80 ##Confidence level for calculation of confidence interval (User Defined)
#' 
#' ## Take Gene Expression and Methylation subtype dataset to create Random Forest Model and Prediction
#' # Cell=NULL ## A list of samples for each data subtype.
#' # Expression=NULL ## List of Matrices where each matrix represent a specific data subtype
#' # finalX=NULL
#' # library(openxlsx)
#' ## Gene Expression
#' # Expression[[1]]=read.xlsx("Dream7_Gene_Expression.xlsx")
#' # Cell[[1]]=colnames(Expression[[1]], do.NULL = TRUE, prefix = "col")[-1]
#' # finalX[[1]]=matrix(as.numeric(t(Expression[[1]])[-1,]),nrow=length(Cell[[1]]))
#' ## Methylation
#' # Expression[[2]]=read.xlsx("Dream7_Methylation.xlsx")
#' # Cell[[2]]=colnames(Expression[[2]], do.NULL = TRUE, prefix = "col")[-1]
#' # finalX[[2]]=matrix(as.numeric(t(Expression[[2]])[-1,]),nrow=length(Cell[[2]]))
#'
#' ## Drug Sensitivity for Modeling 
#' # Drug_Sen_train <- read.xlsx("Dream7_Drug_Sensitivity_Train.xlsx",
#' # +sheet = 1,  startRow = 1, colNames = TRUE)
#' # finalY_train=NULL
#' # for (i in 1:length(Response)){
#' #  Training_response=matrix(Drug_Sen_train[,Response[i]+1],ncol=1)
#' #  Training_response=Imputation(Training_response) ## Imputation for missing values
#' #  finalY_train=cbind(finalY_train,matrix(Training_response,ncol=1)/9.8330) 
#' ## 9.8330 is used for data normalization 
#' # }
#' # finalY_train_cell=Drug_Sen_train[,1] ## Training Sample Index
#' ## Create a list of all combinations of different subtypes of a dataset. 
#' ## Here, for 2 subtypes of dataset,
#' ## then Serial is a list of size 2^2-1=3. The ordering of the three sets are [1 2], [1], [2]
#' # Serial=NULL
#' # library(caTools)
#' # for (p in length(Cell):1){
#' #   nk=combs(1:length(Cell),p)
#' #   sk=length(Serial)
#' #   for (q in 1:dim(nk)[1]){
#' #     Serial[[sk+q]]=nk[q, ]
#' #   }
#' # }
#' ## Combination function gives result of different errors and its' corresponding combination weights.
#' # Result=Combination(finalX,finalY_train,Cell,finalY_train_cell,n_tree,m_feature,
#' # +min_leaf,Serial,Confidence_Level)
#' ## Drug Sensitivity for Prediction
#' # Drug_Sen_test <- read.xlsx("Dream7_Drug_Sensitivity_Test.xlsx", 
#' # +sheet = 1,  startRow = 1, colNames = TRUE)
#' # finalY_test=NULL
#' # for (i in 1:length(Response)){
#' #  Test_response=matrix(Drug_Sen_test[,Response[i]+1],ncol=1)
#' #  Test_response=Imputation(Test_response) ## Imputation for missing values
#' #  finalY_test=cbind(finalY_test,matrix(Test_response,ncol=1)/9.8330) 
#' ## 9.8330 is used for data normalization 
#' # }
#' # finalY_test_cell=Drug_Sen_test[,1] ## Testing Sample Index
#' # Weight=Result[[4]] ##Combination weight using Leave-one-out error is 
#' ## considered for weight of Integrated Prediction
#' ## Gives Integrated Prediction
#' # CombPredict(finalX,finalY_train,Cell,finalY_train_cell,finalY_test_cell,n_tree
#' # +,m_feature,min_leaf,Serial,Weight)
Dream_Dataset_Example <- function(){
  
}
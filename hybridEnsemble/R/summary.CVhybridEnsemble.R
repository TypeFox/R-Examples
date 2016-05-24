#' Summarize the performance of the cross-validated Hybrid Ensemble
#'
#' This function produces summary results per combination method.
#' 
#' @param object An object of class CVhybridEnsemble
#' @param name Name of the dataset. Default is blank.
#' @param stat 'median' or 'IQR' (inter quartile range) of the performance measure used in the CVhybridEnsemble object
#' @param LateX TRUE or FALSE. If true LateX code is printed to the screen. Otherwise a data frame.
#' @param toppart TRUE or FALSE. For the LateX table. Should the top part of the table be printed. Useful for concatenating multiple runs of the \code{summary} function (see examples).
#' @param bottompart TRUE or FALSE. For the LateX table. Should the bottom part of the table be printed. Useful for concatenating multiple runs of the \code{summary} function (see examples).
#' @param ... Not used
#' @details In the output: 'RBGA' (Genetic Algorithm), 'DEOPT' (Differential Evolution), 'GENSA' (Generalized Simulated Annealing), 'MALSCHAINS' (Memetic Algorithm), 'PSOPTIM' (Particle Swarm), 'SOMA' (Self Organizing Migrating Algorithm), 'TABU' (Tabue Search), 'LHNNLS' (Lawson-Hanson Non-negative least squares), 'GINNLS' (Goldfarb-Idnani Non-negative least squares), 'NNloglik' (Non-negative binomial likelihood), 'MEAN' (Simple Mean), 'SB' (Single Best), 'AUTHORITY' (Authority Based method). SB names denote the single best for all cross-validation runs: RF= Random Forest, SV= Bagged Support Vector Machines, KF= Kernel Factory, AB=AdaBoost, LR=Bagged Logistic Regression, NN=Bagged Neural Networks, RoF= Rotation Forest, KN= K-Nearest Neighbors.
#' @examples
#' 
#' 
#' data(Credit)
#' 
#' \dontrun{
#' CVhE <- CVhybridEnsemble(x=Credit[1:200,names(Credit) != 'Response'],
#'                     y=Credit$Response[1:200],
#'                     verbose=TRUE,
#'                     RF.ntree=50,
#'                     KF.rp=1,
#'                     AB.iter=50,
#'                     NN.size=5,
#'                     NN.decay=0,
#'                     SV.gamma = 2^-15,
#'                     SV.cost = 2^-5,
#'                     SV.degree=2,
#'                     SV.kernel='radial')
#' 
#' summary(object=CVhE,stat='median')
#' summary(object=CVhE,stat='IQR')
#' 
#' #LaTeX table
#' #This code example shows how toppart and bottompart can be convenient if you want 
#' #to concatenate multiple datasets (here six time the same dataset).
#' #Paste the output of this code in your LateX document:
#' cat(
#'  summary(object=CVhE ,name="Credit", LateX=TRUE, toppart=TRUE),
#'  summary(object=CVhE ,name="Credit", LateX=TRUE),
#'  summary(object=CVhE, name="Credit", LateX=TRUE),
#'  summary(object=CVhE ,name="Credit", LateX=TRUE),
#'  summary(object=CVhE ,name="Credit", LateX=TRUE),
#'  summary(object=CVhE ,name="Credit", LateX=TRUE, bottompart=TRUE) )
#' 
#' }
#' 
#' @references Ballings, M., Vercamer, D., Van den Poel, D., Hybrid Ensemble: Many Ensembles is Better Than One, Forthcoming.
#' @seealso \code{\link{hybridEnsemble}}, \code{\link{predict.hybridEnsemble}}, \code{\link{importance.hybridEnsemble}}, \code{\link{CVhybridEnsemble}}, \code{\link{plot.CVhybridEnsemble}}
#' @author Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
#' @method summary CVhybridEnsemble
summary.CVhybridEnsemble <- function(object,name='', stat="median", LateX=FALSE, toppart=FALSE, bottompart=FALSE,... ) {
  
  
  stat <- match.arg(stat,c('median','IQR'))
  
  dat <- numeric()
  lab <- character()
  
  SBnames <- character()
  for (i in 1:length(unlist(object$SB$SBname)))  {
   
    if (i<length(unlist(object$SB$SBname))) {
      SBnames[i] <- paste(substr(unlist(object$SB$SBname)[i],1,2),",",sep="") 
    } else {
      SBnames[i] <- paste(substr(unlist(object$SB$SBname)[i],1,2)) 
    }
  }
  SBnames <- paste(SBnames,collapse="")
  
  if (LateX==FALSE) {
    
    for (i in 1:(length(object)-1)) {dat[i] <- round(object[[i]][[stat]],4)  }
    
    for (i in 1:(length(object)-1)) {lab[i] <- paste(substr(names(object)[i],1,4))}
    
    result <- data.frame(t(dat),SBnames, row.names=name)
    colnames(result) <- c(lab,"SB names")
    result  
    
  }else {
  


  for (i in 1:(length(object)-1)) {if (i==length(object)) { dat[i] <- paste("\\small{", round(object[[i]][[stat]],4),"}",sep="") } else {  dat[i] <- paste("\\small{",round(object[[i]][[stat]],4),"} &",sep="")  }  }

  for (i in (1:length(object)-1)) {if (i==length(object)) { lab[i] <- paste("\\small{", substr(names(object)[i],1,4) , "}", sep="") } else {  lab[i] <- paste("\\small{",substr(names(object)[i],1,4),"} &",sep="")  }  }
  
  if (tolower(object$eval.measure)=='auc'){
    measure <- 'AUC'
  } else if (tolower(object$eval.measure)=='sens'){ 
    measure <- 'Sensitivity'
  } else if (tolower(object$eval.measure)=='spec'){ 
    measure <- 'Specificity'
  }

  dat[length(dat)] <-  paste(dat[length(dat)],paste("\\scriptsize{", paste(SBnames,collapse=""),"}",sep=''),sep=' ')
  
  cap <- paste(toupper(substring(stat, 1,1)), substring(stat, 2),sep="", collapse=" ")
  
  if (toppart==TRUE && bottompart==TRUE) {
  
    cat(
  "\\begin{table}[!htbp]
  \\caption{",cap,measure,"of five times twofold cross validation}
  \\centering
  \\begin{tabular}{",rep("l|",length(object)),"p{1.5cm}","}
  \\hline
  \\hline
  \\small{Dataset}  & ", lab,"\\\\
  \\hline \n",name," & ", dat ,"\\\\   
  \\hline
  \\hline
  \\end{tabular}
  \\\\[0.2cm]
  \\scriptsize{S:SVM, R:Random Forest, K:Kernel Factory, L:Logit, A:AdaBoost, N:Neural Network}
  \\label{table:",cap,"}
  \\end{table}" )
    
  } else if (toppart==TRUE && bottompart==FALSE) {
  
    cat(
      "\\begin{table}[!htbp]
  \\caption{",cap," AUC of five times twofold cross validation}
  \\centering
  \\begin{tabular}{",rep("l|",length(object)),"p{1.5cm}","}
  \\hline
  \\hline
  \\small{Dataset}  & ", lab,"\\\\
  \\hline \n",name," & ", dat ,"\\\\   " )
    
  } else if (toppart==FALSE && bottompart==TRUE) {
  
  cat(
  "\n",name," & ", dat ,"\\\\   
  \\hline
  \\hline
  \\end{tabular}
  \\\\[0.2cm]
  \\scriptsize{S:SVM, R:Random Forest, K:Kernel Factory, L:Logit, A:AdaBoost, N:Neural Network}
  \\label{table:",cap,"}
  \\end{table}" )
  
  }else {
  
  cat(
  "\n",name," & ", dat ,"\\\\")
    
  }

  }
  
}
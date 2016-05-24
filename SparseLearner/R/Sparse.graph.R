#' Graphic Modeling Using LASSO-Type Sparse Learning Algorithm.
#'
#' This function builds a gaussian or binary graph based on the bootstrap ranking LASSO regression method.
#' 
#' @param x input matrix. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables. Gaussian or binary data is supported.
#' @param graph.type the type of gaussian or binary graph. Defaults to gaussian.
#' @param B the number of external loop for intersection operation. Defaults to 5.
#' @param Boots the number of internal loop for bootstrap sampling. Defaults to 100.
#' @param edge.rule the rule indicating whether the AND-rule or the OR-rule should be used to define the edges in the graph. Defaults to AND. 
#' @param kfold the number of folds of cross validation - default is 10. Although kfold can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is kfold=3.
#' @param plot logical. Should the resulting graph be plotted? Defaults to TRUE.
#' @param seed the seed for random sampling, with the default value 0123.
#' @export
#' @import glmnet
#' @import SiZer
#' @import qgraph
#' @import SIS
#' @references
#' [1] Guo, P., Zeng, F., Hu, X., Zhang, D., Zhu, S., Deng, Y., & Hao, Y. (2015). Improved Variable 
#' Selection Algorithm Using a LASSO-Type Penalty, with an Application to Assessing Hepatitis B 
#' Infection Relevant Factors in Community Residents. PLoS One, 27;10(7):e0134151.
#' [2] Strobl, R., Grill, E., Mansmann, U. (2012). Graphical modeling of binary data using the LASSO: 
#' a simulation study. BMC Medical Research Methodology, 12:16. 
#' [3] Meinshausen, N., Buehlmann, P. (2006). High-dimensional graphs and variable selection with the Lasso. 
#' Ann Stat, 34:1436-1462.
#' @examples
#' # Example 1: Gene network estimation using the bootstrap ranking LASSO method.
#' # Gaussian graph with OR-rule.
#' library(SIS)
#' data(leukemia.train)
#' # Genes screened by the LASSO algorithm as candidates for graphical modeling.
#' x <- as.matrix(leukemia.train[, -7130])
#' y <- as.numeric(leukemia.train[, 7130])
#' set.seed(0123)
#' cvfit <- cv.glmnet(x=x, y=y, type.measure="deviance", nfolds=3, family="binomial")
#' model.final <- cvfit$glmnet.fit
#' nzero <- as.matrix(coef(model.final, s=cvfit$lambda.min))
#' # To reduce the running time, only half of significant genes are shown.
#' var_nz <- sort(abs(nzero[nzero[,1]!=0, ][-1]), decreasing=TRUE)
#' var_nz <- names(var_nz[1:(length(var_nz)/2)])
#' sub_data <- leukemia.train[, c(var_nz, "V7130")]
#' # Gene expression data subset from patients with acute myeloid leukemia.
#' subset_1 <- subset(sub_data, sub_data$V7130==1)
#' subset_1 <- as.matrix(subset_1[, -dim(subset_1)[2]])
#' # The parameters of B and Boots in the following example are set as small values to
#' # reduce the running time, however the default values are proposed.
#' Sparse.graph.fit1 <- Sparse.graph(subset_1, graph.type=c("gaussian"), 
#'                                    B=2, Boots=1, edge.rule=c("OR"))
#' # Give out the adjacency matrix of variables.
#' Sparse.graph.fit1$adj.matrix

#' # Example 2: Gaussian graph with AND-rule.
#' # The parameters of B and Boots in the following example are set as small values to
#' # reduce the running time, however the default values are proposed.
#' Sparse.graph.fit2 <- Sparse.graph(subset_1, graph.type=c("gaussian"), 
#'                         B=2, Boots=1, edge.rule=c("OR"), plot=FALSE)
#' # Give out the adjacency matrix of variables.
#' Sparse.graph.fit2$adj.matrix
#' # Plot the graph based on the adjacency matrix of variables using the qgraph package.
#' library(qgraph)
#' qgraph(Sparse.graph.fit2$adj.matrix, directed=FALSE, color="blue", 
#'         negCol="red", edge.labels=TRUE, layout="circle")
Sparse.graph=function(x, graph.type=c("gaussian"), B=5, Boots=100, edge.rule=c("AND"), kfold=10, plot=TRUE, seed=0123){
   RM <- NULL; Neighmatrix <- NULL; ii <- NULL; Mymatrix <-NULL
   RM <- matrix(0, ncol(x), ncol(x))
   colnames(RM) <- colnames(x)
   rownames(RM) <- colnames(x)
   Neighmatrix <- RM
   set.seed(seed)

 if (graph.type==c("gaussian")) {
  for(ii in 1: ncol(x)){
    y0 <- NULL; x0 <- NULL; varx <- NULL; varx0 <- NULL; nvar <- NULL; rowx <- NULL; n <- NULL;
    Nboot <- NULL; s <- NULL; BoostrapX <- NULL; BoostrapY <- NULL; cvfit <- NULL; model.final <- NULL;
    nzero <- NULL; var_nz <- NULL; Impor <- NULL; score <- NULL; index <- NULL; model.pwl <- NULL;
    break_p <- NULL; Impor0 <- NULL; Myresult <- NULL; var_nonzero <- NULL; RecordM2 <- NULL;
    RecordM <- NULL; var_list_all <- NULL; Adjmatrix <- NULL; Var.selected <- NULL; Var.coef <- NULL;
    ij <- NULL; i <- NULL; i0 <- NULL; j0 <- NULL; i1<-NULL; j1<-NULL

    y0=as.numeric(x[,ii])
    x0=as.matrix(x[,-ii])
  #########
    varx <- colnames(x0)
    varx0 <- varx
    nvar <- length(varx)
    rowx <- nrow(x0)
    n <- length(y0)
    var_list_all <- vector("list", B)
    if (rowx!=n){
      stop("The number of rows in predictor matrix is not equal to the length of response vector!")
        }
    if (nvar==1){
      stop("The x matrix has only one variable!")
        }
    RecordM2 <- matrix(0, nrow=B, ncol=nvar)
    rownames(RecordM2) <- paste("B", 1:B, sep="")
    colnames(RecordM2) <- varx
  for (ij in 1:B){
    Nboot <- n
    RecordM <- matrix(0, nrow=Boots, ncol=nvar)
    rownames(RecordM) <- paste("BSample", 1:Boots, sep="")
    colnames(RecordM) <- varx 
    for(i in 1:Boots){
        s <- sample(Nboot, replace=TRUE)
        BoostrapX <- x0[s, ]
        BoostrapY <- y0[s]
        cvfit <- cv.glmnet(BoostrapX, BoostrapY, nfolds=kfold, family="gaussian")
        model.final <- cvfit$glmnet.fit
        nzero <- as.matrix(coef(model.final, s=cvfit$lambda.min))
        var_nz <- names(nzero[nzero[,1]!=0,])
        for(i0 in 1:dim(nzero)[1]){
          for(j0 in 1:nvar){
            if (dimnames(nzero)[[1]][i0]==colnames(RecordM)[j0]){RecordM[i, j0] <- nzero[i0, ]}
                            }
                                 }
        cat("Step 1-Current iteration: ", i, "\n")
                      }
    Impor <- abs(apply(RecordM, 2, mean))
    score <- sort(Impor)
    index <- seq(1:length(score))
    model.pwl <- piecewise.linear(index, score, middle=1, CI=FALSE, bootstrap.samples=1000, sig.level=0.05)
    break_p <- round(model.pwl$change.point)
    var_list_all[[ij]] <- names(score[break_p:length(index)])
    if (length(var_list_all[[ij]])!=0){
      varx0 <- intersect(var_list_all[[ij]], varx0)}
    Impor0 <- apply(RecordM, 2, mean)[var_list_all[[ij]]]
    for(gp in 1:length(Impor0)){
      RecordM2[ij, colnames(RecordM2)==names(Impor0)[gp]]=Impor0[gp]
                                }
    cat("Step 2-Current iteration: ", ij, "\n")
                   }
    if(length(varx0)==1){
      RecordM2 <- as.matrix(RecordM2[, varx0], ncol=1)
      colnames(RecordM2) <- varx0
                         } else {
      RecordM2 <- RecordM2[, varx0]
                                }
    Var.selected <- colnames(RecordM2)
    Var.coef <- apply(RecordM2, 2, mean)
    Myresult <- list(Var.selected=Var.selected, Var.coef=Var.coef)
  #########
        var_nonzero <- Myresult$Var.selected
        if (length(var_nonzero)!=0){
          for(i0 in 1:length(var_nonzero)){
            for(j0 in 1:ncol(RM)){
              if (var_nonzero[i0]==colnames(RM)[j0]) RM[ii,j0]=Myresult$Var.coef[i0]
              else RM[ii,j0]=RM[ii,j0]
                                 }
                                           }
                                    }
        cat("Iteration = ", ii, "\n")
                      }
    if (edge.rule==c("AND")) {
                             for(i1 in 1:nrow(RM)){
                               for(j1 in 1:ncol(RM)){
                                 if(RM[i1,j1]!=0 & RM[j1,i1]!=0) Neighmatrix[i1,j1] <- (RM[i1,j1]+RM[j1,i1])/2
                                 else Neighmatrix[i1,j1] <- Neighmatrix[i1,j1]
                                                     }
                                                  }
    Neighmatrix[lower.tri(Neighmatrix)] <- 0
                            }
    if (edge.rule==c("OR")) {
                             for(i1 in 1:nrow(RM)){
                               for(j1 in 1:ncol(RM)){
                                 if(RM[i1,j1]!=0 | RM[j1,i1]!=0) Neighmatrix[i1,j1]=RM[i1,j1]+RM[j1,i1]
                                 else Neighmatrix[i1,j1]=Neighmatrix[i1,j1]
                                                     }
                                                   }
    Neighmatrix[lower.tri(Neighmatrix)] <- 0
                           }
                              }

 if (graph.type==c("binary")) {
  for(ii in 1: ncol(x)){
    y0 <- NULL; x0 <- NULL; varx <- NULL; varx0 <- NULL; nvar <- NULL; rowx <- NULL; n <- NULL;
    Nboot <- NULL; s <- NULL; BoostrapX <- NULL; BoostrapY <- NULL; cvfit <- NULL; model.final <- NULL;
    nzero <- NULL; var_nz <- NULL; Impor <- NULL; score <- NULL; index <- NULL; model.pwl <- NULL;
    break_p <- NULL; Impor0 <- NULL; Myresult <- NULL; var_nonzero <- NULL; RecordM2 <- NULL;
    RecordM <- NULL; var_list_all <- NULL; Adjmatrix <- NULL; Var.selected <- NULL; Var.coef <- NULL;
    ij <- NULL; i <- NULL; i0 <- NULL; j0 <- NULL; i1<-NULL; j1<-NULL;

    y0 <- as.numeric(x[, ii])
    x0 <- as.matrix(x[, -ii])
  #########
    varx <- colnames(x0)
    varx0 <- varx
    nvar <- length(varx)
    rowx <- nrow(x0)
    n <- length(y0)
    var_list_all <- vector("list",B)
    if (rowx!=n){
      stop("The number of rows in predictor matrix is not equal to the length of response vector!")
                }
    if (nvar==1){
      stop("The x matrix has only one variable!")
                 }
    RecordM2 <- matrix(0, nrow=B, ncol=nvar)
    rownames(RecordM2) <-paste("B", 1:B, sep="")
    colnames(RecordM2) <-varx
  for (ij in 1:B){
    Nboot <- n
    RecordM <- matrix(0, nrow=Boots, ncol=nvar)
    rownames(RecordM) <- paste("BSample", 1:Boots, sep="")
    colnames(RecordM) <- varx 
    for(i in 1:Boots){
      repeat{ 
        s <- sample(Nboot, replace=TRUE)  
        if(length(table(y0[s])) >= 2 & length(table(y0[-s])) >= 2)
        break
             }
      BoostrapX <- x0[s, ]
      BoostrapY <- y0[s]
      cvfit <- cv.glmnet(BoostrapX, BoostrapY, type.measure="deviance", nfolds=kfold, family="binomial")
      model.final <- cvfit$glmnet.fit
      nzero <- as.matrix(coef(model.final, s=cvfit$lambda.min))
      var_nz <- names(nzero[nzero[,1]!=0, ])
      for(i0 in 1:dim(nzero)[1]){
        for(j0 in 1:nvar){
           if (dimnames(nzero)[[1]][i0]==colnames(RecordM)[j0]){RecordM[i, j0] <- nzero[i0, ]}
                          }
                                 }
      cat("Step 1-Current iteration: ", i, "\n")
                      }
    Impor <- abs(apply(RecordM, 2, mean))
    score <- sort(Impor)
    index <- seq(1:length(score))
    model.pwl <- piecewise.linear(x=index, y=score, middle=1, CI=FALSE, bootstrap.samples=1000, sig.level=0.05)
    break_p <- round(model.pwl$change.point)
    var_list_all[[ij]] <- names(score[break_p:length(index)])
    if (length(var_list_all[[ij]])!=0){
      varx0 <- intersect(var_list_all[[ij]], varx0)}
    Impor0 <- apply(RecordM, 2, mean)[var_list_all[[ij]]]
    for(gp in 1:length(Impor0)){
      RecordM2[ij, colnames(RecordM2)==names(Impor0)[gp]]=Impor0[gp]
                                }
    cat("Step 2-Current iteration: ", ij, "\n")
                     }
    if(length(varx0)==1){
      RecordM2 <- as.matrix(RecordM2[, varx0], ncol=1)
      colnames(RecordM2) <- varx0
                         } else {
      RecordM2 <- RecordM2[, varx0]
                                }
    Var.selected=colnames(RecordM2)
    Var.coef=apply(RecordM2, 2, mean)
    Myresult <- list(Var.selected=Var.selected, Var.coef=Var.coef)
  #########
        var_nonzero <- Myresult$Var.selected
        if (length(var_nonzero)!=0){
          for(i0 in 1:length(var_nonzero)){
            for(j0 in 1:ncol(RM)){
              if (var_nonzero[i0]==colnames(RM)[j0]) RM[ii,j0] <- Myresult$Var.coef[i0]
              else RM[ii,j0] <- RM[ii,j0]
                                 }
                                           }
                                    }
        cat("Iteration = ", ii, "\n")
                      }
    if (edge.rule==c("AND")) {
                             for(i1 in 1:nrow(RM)){
                               for(j1 in 1:ncol(RM)){
                                 if(RM[i1,j1]!=0 & RM[j1,i1]!=0) Neighmatrix[i1,j1] <- (RM[i1,j1]+RM[j1,i1])/2
                                 else Neighmatrix[i1,j1] <- Neighmatrix[i1,j1]
                                                     }
                                                  }
    Neighmatrix[lower.tri(Neighmatrix)] <- 0
                            }
    if (edge.rule==c("OR")) {
                             for(i1 in 1:nrow(RM)){
                               for(j1 in 1:ncol(RM)){
                                 if(RM[i1,j1]!=0 | RM[j1,i1]!=0) Neighmatrix[i1,j1] <- RM[i1,j1]+RM[j1,i1]
                                 else Neighmatrix[i1,j1] <- Neighmatrix[i1,j1]
                                                     }
                                                   }
    Neighmatrix[lower.tri(Neighmatrix)] <- 0
                           }
                               }
   Neighmatrix <- as.matrix(Neighmatrix)
  if(plot){
      qgraph(Neighmatrix, edge.color=NA, directed=FALSE, node.width=1, node.height=1, color="blue", borders=FALSE, 
      vTrans=100, label.cex=10, label.color="black", labels=colnames(Neighmatrix), label.prop=.9, posCol="black", 
      negCol="red", edge.width=1, edge.labels=TRUE, layout="circle", edge.label.position=0.5)
          }
   Mymatrix <- list(adj.matrix=Neighmatrix, graph.type=graph.type, B=B, Boots=Boots, edge.rule=edge.rule)
   return(Mymatrix)
}

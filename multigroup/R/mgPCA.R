#' @title Multigroup Principal Component Analysis
#' 
#' @description 
#' Multigroup PCA algorithm (NIPALS for Multigroup PCA)
#' 
#' @param Data a numeric matrix or data frame
#' @param Group a vector of factors associated with group structure
#' @param ncomp number of components, if NULL number of components is equal to 2
#' @param Scale scaling variables, by defalt is FALSE. By default data are centered within groups
#' @param graph should loading and component be plotted
#' @return list with the following results:
#' @return \item{Data}{     Original data}
#' @return \item{Con.Data}{     Concatenated centered data}
#' @return \item{split.Data}{     Group centered data}
#' @return \item{Group}{      Group as a factor vector}
#' @return \item{loadings.group}{     Loadings associated with each group}
#' @return \item{score.group}{      Scores associated with each group}
#' @return \item{loadings.common}{      Matrix of common loadings}
#' @return \item{score.Global}{     Global scores}
#' @return \item{cumper.inertigroup}{     Cumulative percentage of group components inertia}
#' @return \item{cumper.inertiglobal}{      Cumulative percentage of global component inertia}
#' @return \item{noncumper.inertiglobal}{     Percentage of global component inertia}
#' @return \item{lambda}{     The specific variances of groups}
#' @return \item{exp.var}{      Percentages of total variance recovered associated with each dimension }
#' @return \item{Similarity.Common.Group.load}{Cumulative similarity between group and common loadings}
#' @return \item{Similarity.noncum.Common.Group.load}{ NonCumulative  similarity between group and common loadings}
#' @seealso \code{\link{BGC}}, \code{\link{FCPCA}}, \code{\link{DCCSWA}}, \code{\link{DSTATIS}}, \code{\link{DGPA}}, \code{\link{summarize}}, \code{\link{TBWvariance}}, \code{\link{loadingsplot}}, \code{\link{scoreplot}}, \code{\link{iris}}  
#' @export
#' 
#' 
#' @references A. Eslami, E. M. Qannari, A. Kohler and S. Bougeard (2013). General overview
#'  of methods of analysis of multi-group datasets,
#'  \emph{Revue des Nouvelles Technologies de l'Information}, 25, 108-123.
#' 
#' @references A. Eslami, E. M. Qannari, A. Kohler and S. Bougeard (2013). Analyses
#' factorielles de donnces structurces en groupes d'individus,
#' \emph{Journal de la Societe Francaise de Statistique}, 154(3), 44-57.
#'    
#'       
#' @examples
#' Data = iris[,-5]
#' Group = iris[,5]
#' res.mgPCA = mgPCA (Data, Group)
#' barplot(res.mgPCA$noncumper.inertiglobal)
#' #----------------
#' #Similarity index: group loadings are compared to the common structure (first dimension)
#' Xzero = rep(0, 3)
#' MIN = min(res.mgPCA$Similarity.noncum.Common.Group.load[[1]][-1, 1])-0.0005
#' XLAB = paste("Dim1, %",res.mgPCA$noncumper.inertiglobal[1])
#' plot(Xzero, res.mgPCA$Similarity.noncum.Common.Group.load[[1]][-1, 1], pch=15, ylim=c(MIN, 1), 
#' main="Similarity between groups and common structure", xlab=XLAB, ylab="", xaxt="n")
#' abline(v=0)
#' abline(h=seq(MIN, 1, by=0.05), col="black", lty=3)
#' XX=res.mgPCA$Similarity.noncum.Common.Group.load[[1]][-1, 1, drop=FALSE]
#' text(Xzero, XX, labels=rownames(XX), pos=4)
#' #----------------
#' # Similarity index: group loadings are compared to the common structure (dimensions 1 and 2)
#' XX1=res.mgPCA$Similarity.noncum.Common.Group.load[[1]][-1, 1]
#' XX2=res.mgPCA$Similarity.noncum.Common.Group.load[[2]][-1, 1]
#' simil <- cbind(XX1, XX2)
#' YLAB = paste("Dim1, %",res.mgPCA$noncumper.inertiglobal[2])
#' plot(simil, xlab=XLAB, ylab=YLAB, main="Similarity between groups and common structure", pch=20)
#' text(simil, labels=rownames(simil), cex=1, font.lab=1, pos=3)
#' #------------------
#' loadingsplot(res.mgPCA, axes=c(1,2), INERTIE=res.mgPCA$noncumper.inertiglobal)
#' scoreplot(res.mgPCA, axes=c(1,2))
mgPCA <- function(Data, Group, ncomp=NULL, Scale=FALSE, graph=FALSE){


  #=========================================================================
  #                             1. Checking the inputs
  #=========================================================================
  check(Data, Group)
  
  
  #=========================================================================
  #                              2. preparing Data
  #=========================================================================
  if (class(Data) == 'data.frame') {
    Data=as.matrix(Data)
  }
  if(is.null(ncomp)) {ncomp=2}
  if(is.null(colnames(Data))) {
    colnames(Data) = paste('V', 1:ncol(Data), sep='')
  }
  Group = as.factor(Group)
  
  
  
  rownames(Data) = Group                 #---- rownames of data=groups
  M = length(levels(Group))              #----number of groups: M
  P = dim(Data)[2]                       #----number of variables: P
  n = as.vector(table(Group))            #----number of individuals in each group
  N = sum(n)                             #----number of individuals
  split.Data = split(Data,Group)         #----split Data to M parts 
  
  # centering and scaling if TRUE
  for(m in 1:M){  
    split.Data[[m]] = matrix(split.Data[[m]], nrow=n[m])
    split.Data[[m]] = scale(split.Data[[m]], center=TRUE, scale=Scale)
  }
  
  # concatinated dataset by row as groups
  Con.Data = split.Data[[1]]  
  for(m in 2:M) {
    Con.Data = rbind(Con.Data, split.Data[[m]])
  }
  rownames(Con.Data) = Group
  colnames(Con.Data) = colnames(Data)
  
  # Variance-covariance matrix for each group
  cov.Group = vector("list", M)
  for(m in 1:M){    
    cov.Group[[m]] = t(split.Data[[m]]) %*% split.Data[[m]] / n[m]
  }
  #==========================================================================
  #  				                      3. Outputs
  #==========================================================================
  res <- list(
    Data       = Data,
    Con.Data   = Con.Data,
    split.Data = split.Data,
    Group=Group)
  
  
  
  res$loadings.group = vector("list", M)
  res$score.group = vector("list", M)
  for(m in 1:M){
    res$score.group[[m]]=matrix(0, nrow=n[m], ncol=ncomp)
    #rownames(res$score.group[[m]])=Group
    colnames(res$score.group[[m]]) = paste("Dim", 1:ncomp, sep="")
    
    
    res$loadings.group[[m]] = matrix(0, nrow=P, ncol=ncomp)
    rownames(res$loadings.group[[m]]) = colnames(Data)
    colnames(res$loadings.group[[m]]) = paste("Dim", 1:ncomp, sep="")
  }
  
  res$loadings.common = matrix(0, nrow=P, ncol=ncomp)  
  rownames(res$loadings.common) = colnames(Data)
  colnames(res$loadings.common) = paste("Dim", 1:ncomp, sep="")
  
  
  res$score.Global = matrix(0, nrow=N, ncol=ncomp)     
  rownames(res$score.Global) = Group
  colnames(res$score.Global) = paste("Dim", 1:ncomp, sep="")
  
  
  res$Group.weight=matrix(0, nrow=M, ncol=ncomp) 
  rownames(res$Group.weight)=levels(Group)  
  colnames(res$Group.weight)=paste("Dim", 1:ncomp, sep="")
  
  res$cumper.inertigroup=matrix(0, ncol=ncomp,nrow=M)
  rownames(res$cumper.inertigroup)=levels(Group)
  colnames(res$cumper.inertigroup)=paste("Dim", 1:ncomp, sep="")
  
  res$cumper.inertiglobal=matrix(0, ncol=ncomp, nrow=1)
  colnames(res$cumper.inertiglobal)=paste("Dim", 1:ncomp, sep="")
  
  res$noncumper.inertiglobal = matrix(0, ncol=ncomp, nrow=1)
  colnames(res$noncumper.inertiglobal) = paste("Dim", 1:ncomp, sep="")
  #=========================================================================
  #                          4. NIPALS for multi_group CPCA
  #=========================================================================
  eps = 1e-7
  
  for(h in 1:ncomp){ 		          #------ for iteration loops
    threshold = 1.0;
    iNIPALS = 0;
    critt = 0  
    
    #=======================  4.1 initial value of vector of common loadings
    w.common = rnorm(P)       
    w.common = normv(w.common)
    W=matrix(0, ncol=P, nrow=M)     #------ combined group loadings
      
    
    #======================= 4.2 starting iteration
    while (threshold > eps) {
      
      t.global = matrix(0, ncol=1)  # concatinated components
      
      
      # 4.2.1 group analysis
      for (m in 1:M){     
        t.group = (split.Data[[m]]) %*% w.common
        t.group = normv(t.group)
        res$score.group[[m]][,h] = t.group
        t.global=rbind(t.global,t.group)
        
        w.group=t(split.Data[[m]])%*%t.group             # group loadings
        res$loadings.group[[m]][,h]= normv(w.group)      #w.group
        W[m,] = w.group                                # combined group loadings
      }
      
      res$score.Global[,h]= Con.Data %*% w.common   
      
      # 4.2.2 new value of common loadings (w.common)
      Group.weight = W %*% w.common   # global loading weights, W is projected to w
      Group.weight = normv(Group.weight)
      res$Group.weight[,h] = Group.weight
      w.common = t(t(Group.weight) %*% W )
      w.common = normv(w.common)
      
      
      # 4.2.3  criterion  
      iNIPALS=iNIPALS+1
      Crit=0
      for(m in 1:M){
        Crit= Crit+(t(res$score.group[[m]][,h])%*% 
                          (split.Data[[m]]) %*% w.common %*% t(w.common) 
                        %*% t(split.Data[[m]]) %*%
                          res$score.group[[m]][,h] )
      }
      critt[iNIPALS]=Crit
      
      if (iNIPALS>1){
        threshold=critt[iNIPALS]-critt[(iNIPALS-1)]
      }
      res$loadings.common[,h]=w.common
           
    } # end of iteration
    
    #=========================================================================
    #                               5. Explained variance 
    #=========================================================================
    for(m in 1:M){
      project.group=res$score.group[[m]][,1:h] %*% ginv( t(res$score.group[[m]][,1:h]) %*% res$score.group[[m]][,1:h]) %*%   t(res$score.group[[m]][,1:h])
      projet.data.group=project.group %*% res$split.Data[[m]]
      res$cumper.inertigroup[m,h]=100*sum(diag(t(projet.data.group) %*% projet.data.group))/ sum(diag(t(res$split.Data[[m]]) %*% res$split.Data[[m]]))
    }
 
    project.global=res$score.Global[,1:h] %*% ginv( t(res$score.Global[,1:h]) %*% res$score.Global[,1:h]) %*%   t(res$score.Global[,1:h])
    projet.data.global=project.global%*% res$Con.Data
    res$cumper.inertiglobal[h]=100*sum(diag(t(projet.data.global) %*% projet.data.global))/ sum(diag(t(res$Con.Data) %*% res$Con.Data))
    
    projectn.global=res$score.Global[,h] %*% ginv( t(res$score.Global[,h]) %*% res$score.Global[,h]) %*%   t(res$score.Global[,h])
    projetn.data.global=projectn.global %*% res$Con.Data
    res$noncumper.inertiglobal[h] = round(100*sum(diag(t(projetn.data.global) %*% projetn.data.global))/ sum(diag(t(res$Con.Data) %*%res$Con.Data)),1)
    #=========================================================================
    #                                   Deflation 
    #=========================================================================
    Con.Data= deflation(Con.Data, res$score.Global[,h])
    rownames(Con.Data) = Group
    split.Data = split(Con.Data,Group)  
    
    for(m in 1:M){  
      split.Data[[m]] = matrix(split.Data[[m]],ncol=P)
      colnames(split.Data[[m]]) = colnames(Con.Data)
    }
    
  }  #------ END OF DIMENSION
  
  lambda = matrix(0, nrow=M, ncol=ncomp)
  for(m in 1:M){
    lambda[m,]=round(diag(t(res$loadings.common) %*% cov.Group[[m]] %*% res$loadings.common),3)
  }
  res$lambda = lambda
  rownames(res$lambda) = levels(Group)
  colnames(res$lambda) = paste("Dim", 1:ncomp, sep="")

  exp.var = matrix(0,M,ncomp)
  for(m in 1:M){
    exp.var[m,] = 100 * lambda[m,]/ sum(diag(cov.Group[[m]]))
  }
  res$exp.var = exp.var
  rownames(res$exp.var) = levels(Group)
  colnames(res$exp.var) = paste("Dim", 1:ncomp, sep="")
  ##----------------------------------------------------------------------------
  ##  		       Similarity between partial and common loadings 
  ##----------------------------------------------------------------------------
  loadings_matrices = list()
  loadings_matrices[[1]] = res$loadings.common
  for(m in 2:(M+1)){
    loadings_matrices[[m]] = res$loadings.group[[(m-1)]]
  }
  NAMES = c("Commonload", levels(Group))
  
  Similarity =similarity_function(loadings_matrices=loadings_matrices, NAMES)
  Similarity_noncum = similarity_noncum(loadings_matrices=loadings_matrices, NAMES)

  res$Similarity.Common.Group.load          = Similarity
  res$Similarity.noncum.Common.Group.load   = Similarity_noncum
  #--------------------------------------
  if(graph) {plot.mg(res)}
  
  # add class
  class(res) = c("mgpca", "mg")
  return(res)
}

#' @S3method print mgpca
print.mgpca <- function(x, ...)
{
  cat("\nMulti-group Principal Component Analysis\n")
  cat(rep("-",43), sep="")
  cat("\n$score.group       ", "score groups")
  cat("\n$loadings.common   ", "common loadings")
  cat("\n$Data              ", "Data set")
  cat("\n")
  invisible(x)
}
  

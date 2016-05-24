#' @title Dual Common Component and Specific Weights Analysis
#' 
#' @description 
#' Dual Common Component and Specific Weights Analysis: to find 
#' common structure among variables of different groups
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
#' @return \item{loadings.common}{      Matrix of common loadings} 
#' @return \item{saliences}{      Each group having a specific contribution to the determination of this common
#' space, namely the salience, for each dimension under study} 
#' @return \item{lambda}{     The specific variances of groups}
#' @return \item{exp.var}{      Percentages of total variance recovered associated with each dimension}
#' @seealso \code{\link{mgPCA}}, \code{\link{FCPCA}}, \code{\link{BGC}}, \code{\link{DSTATIS}}, \code{\link{DGPA}}, \code{\link{summarize}}, \code{\link{TBWvariance}}, \code{\link{loadingsplot}}, \code{\link{scoreplot}}, \code{\link{iris}}  
#' @export
#' @references E. M. Qannari, P. Courcoux, and E. Vigneau (2001). Common components and specific weights analysis performed 
#'  on preference data. \emph{Food Quality and Preference}, 12(5-7), 365-368.
#'  
#'  
#' @references A. Eslami (2013). Multivariate data analysis of multi-group datasets: application to biology. University  of Rennes I.
#'  
#'  
#'    
#' @examples
#' Data = iris[,-5]
#' Group = iris[,5]
#' res.DCCSWA = DCCSWA(Data, Group, graph=TRUE)
#' loadingsplot(res.DCCSWA, axes=c(1,2))
#' scoreplot(res.DCCSWA, axes=c(1,2)) 
DCCSWA <- function(Data, Group, ncomp=NULL, Scale=FALSE, graph=FALSE){

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
  
  
  
  rownames(Data)=Group                 #---- rownames of data=groups
  M=length(levels(Group))              #----number of groups: M
  P=dim(Data)[2]                       #----number of variables: P
  n=as.vector(table(Group))            #----number of individuals in each group
  N=sum(n)                             #----number of individuals
  split.Data=split(Data,Group)         #----split Data to M parts 
  
   
  # centering and scaling if TRUE
  for(m in 1:M){  
    split.Data[[m]]=matrix(split.Data[[m]],nrow=n[m])
    split.Data[[m]]<-scale(split.Data[[m]],center=TRUE,scale=Scale)
  }

    # concatinated dataset by row as groups
  Con.Data = split.Data[[1]]  
  for(m in 2:M) {
    Con.Data = rbind(Con.Data, split.Data[[m]])
  }
  rownames(Con.Data) = Group
  colnames(Con.Data) = colnames(Data)

  # Variance-covariance matrix for each group
  cov.Group=vector("list", M)
  for(m in 1:M){    
    cov.Group[[m]] = t(split.Data[[m]]) %*% split.Data[[m]] / n[m]
  }
  
  # Total variance of all dataset sum(trace(Vm*Vm)
  Itot = 0      
  for (m in 1:M) {
    Itot = Itot + sum(cov.Group[[m]]^2)
  }
  
  #==========================================================================
  #    			                      Outputs
  #==========================================================================
  saliences = matrix(0,M,ncomp)               #----saliences
  A         = matrix(0,nrow=P,ncol=ncomp)     #----common loadings matrix
  explained = matrix(0,ncomp) 
  res <- list(
    Data       = Data,
    Con.Data   = Con.Data,
    split.Data = split.Data,
    Group=Group)
  
  #============================================================================
  #      		                      3. Iterative algorithm
  #   computation of matrix of common loading (A) and saliences 
  #============================================================================
   for (h in 1:ncomp)  {
    
    # Initialization: Choose initial values for salience=1
    salience= c(rep(1,M))      
    threshold = 1e-10
    deltafit = 1000000;
    previousfit = Itot;
    
    while(deltafit > threshold){
      
      
      # Compute compromise varaiance-covariance matrix
      Vc = 0
      for(m in 1:M){
        Vc = salience[m] * cov.Group[[m]] + Vc
      }
      
      
      #Set a to the eigenvector of Vc associated with the largest eigenvalue
      a = eigen(Vc)$vectors[,1]
      
      #Update new saliences 
      for(m in 1:M){
        salience[m] = t(a) %*% cov.Group[[m]] %*% a
      }
      
      # criterion
      def=0
      for(m in 1:M){
        def = sum(cov.Group[[m]] - salience[m] * as.matrix(a) %*% t(as.matrix(a)))^2 + def
      }
      deltafit = previousfit-def
      previousfit = def
    } # end of iteration  
    
    explained[h,1] = 100 * sum(salience ^2) / Itot
    
    saliences[,h] = salience;
    A[,h] = a;
    
    #==========================================================================
    #                               Deflation
    #==========================================================================
    
    aux = diag(1,P)- matrix(a,ncol=1) %*% matrix(a,nrow=1);
    for (m in 1:M)   {
      split.Data[[m]] = split.Data[[m]] %*% aux
      cov.Group[[m]]  = t(split.Data[[m]]) %*% split.Data[[m]]
    }
  }
  
   
  res$loadings.common=A
  rownames(res$loadings.common)=colnames(Data)
  colnames(res$loadings.common)=paste("Dim", 1:ncomp, sep="")
  
  
  res$saliences= saliences  
  rownames(res$saliences) = levels(Group)
  colnames(res$saliences) = paste("Dim", 1:ncomp, sep="")
  
  
  expl = matrix(0, nrow=ncomp, ncol=2)
  expl[,1] = c(explained)
  expl[,2] = cumsum(expl[,1]) 
  res$exp.variance = expl
  rownames(res$exp.variance)  = paste("Dim", 1:ncomp, sep="")
  colnames(res$exp.variance)  = c("%Total Var expl", "Cumul % total Var")

  #============================================================================
  #                                Explained variance 
  #============================================================================
  
  Y<-scale(Data,center=TRUE,scale=FALSE)
  Datam=split(Y, Group)
  for(m in 1:M){
    Datam[[m]] = matrix(Datam[[m]], nrow=n[m])
    Datam[[m]] = scale(Datam[[m]], center=TRUE, scale=FALSE)
  }

  for(m in 1:M) {
    cov.Group[[m]] = t(Datam[[m]]) %*% Datam[[m]] / n[m]
  }

  # variance of each loading
  # lambda = t(common loading)*(t(Xm)* Xm) * common loading
  lambda = matrix(0, nrow=M, ncol=ncomp)
  
  for(m in 1:M){
    lambda[m,]=round(diag(t(A) %*%  (cov.Group[[m]]) %*% A),3)
  }
  res$lambda = lambda
  rownames(res$lambda) = levels(Group)
  colnames(res$lambda) = paste("Dim", 1:ncomp, sep="")
 
  # exp.var 
  exp.var = matrix(0,M,ncomp)
  for(m in 1:M){
    exp.var[m,] = 100 * lambda[m,]/ sum(diag(cov.Group[[m]]))
  }
  res$exp.var = exp.var
  rownames(res$exp.var) = levels(Group)
  colnames(res$exp.var) = paste("Dim", 1:ncomp, sep="")
  
  #=================================================
  if(graph) {plot.mg(res)}

 # add class
  class(res) = c("DCCSWA", "mg")
  return(res)
}


#' @S3method print DCCSWA
print.DCCSWA <- function(x, ...)
{
  cat("\nDual Common Component and Specific Weights Analysis\n")
  cat(rep("-",43), sep="")
  cat("\n$loadings.common   ", "common loadings")
  cat("\n$Data              ", "Data set")
  cat("\n")
  invisible(x)
}
#' @title Dual Generalized Procrustes Analysis
#' 
#' @description 
#' Dual Generalized Procrustes Analysis to study multigroup data
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
#' @return \item{lambda}{     The specific variances of groups}
#' @return \item{exp.var}{      Percentages of total variance recovered associated with each dimension }
#' @seealso \code{\link{mgPCA}}, \code{\link{FCPCA}}, \code{\link{DCCSWA}}, \code{\link{DSTATIS}}, \code{\link{BGC}}, \code{\link{summarize}}, \code{\link{TBWvariance}}, \code{\link{loadingsplot}}, \code{\link{scoreplot}}, \code{\link{iris}}  
#' @export
#' @references J. Gower (1975). Generalized procrustes analysis. \emph{Psychometrika}, 40(1), 3-51.
#' 
#' @references A. Eslami, E. M. Qannari, A. Kohler and S. Bougeard (2013). General overview
#'  of methods of analysis of multi-group datasets,
#'  \emph{Revue des Nouvelles Technologies de l'Information}, 25, 108-123.
#' 
#'  @references A. Eslami, E. M. Qannari, A. Kohler and S. Bougeard (2013). Analyses
#' factorielles de donnees structurees en groupes d'individus,
#' \emph{Journal de la Societe Francaise de Statistique}, 154(3), 44-57.
#'   
#'    
#'      
#' @examples
#' Data = iris[,-5]
#' Group = iris[,5]
#' res.DGPA = DGPA(Data, Group, graph=TRUE)
#' loadingsplot(res.DGPA, axes=c(1,2))
#' scoreplot(res.DGPA, axes=c(1,2)) 
DGPA <- function(Data, Group, ncomp=NULL, Scale=FALSE, graph=FALSE){
    
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
  if(is.null(ncomp)) {ncomp=2}  # or H=qr(g_data_group)$rank
  if(is.null(colnames(Data))) {
    colnames(Data) = paste('V', 1:ncol(Data), sep='')
  }
  Group = as.factor(Group)
  
  rownames(Data) = Group                 #---- rownames of data=groups
  M = length(levels(Group))              #----number of groups: M
  P = dim(Data)[2]                       #----number of variables: P
  n = as.vector(table(Group))            #----number of individuals in each group
  N = sum(n)                             #----number of individuals
  nmax=max(n)
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
  
  # Transposed dataset and divided by sqrt(n_m)
  tab=vector("list", M)    
  for (m in 1:M){
    tab[[m]] = t(split.Data[[m]]) / sqrt(n[m])
    nc = ncol(tab[[m]])  
    if (nc < nmax)  
    {tab[[m]] = cbind(tab[[m]], matrix(0,ncol=(nmax-nc), nrow=P))}
  }
  
  # Total variance of all dataset sum(trace(Vm*Vm) 
  Itot=0      
  for (m in 1:M) {
    Itot = Itot +   sum(diag(cov.Group[[m]]))
  }
  #============================================================================
  #    			                       Outputs
  #============================================================================
  res <- list(
    Data       = Data,
    Con.Data   = Con.Data,
    split.Data = split.Data,
    Group=Group)
  
  #============================================================================
  #    			                     Procrustes function
  #                   for finding the rotaion matrices and distances
  #============================================================================
  gpa <- function(A,B){
    A <- as.matrix(A)
    B <- as.matrix(B)  
    v.sing <- svd(t(A)%*%B)            # svd(G)=uDv'
    H <- v.sing$u %*% t(v.sing$v)      # rotation          
    d <- sum((A %*% H - B)^2)   #distance|AH-B|=trac((AH-B)'AH-B)=sum((AH-B)^2)
    list(H=H, d=d)
  }
  #============================================================================
  #      		                    Iterative algorithm 
  #============================================================================
  # finding: C  compromis:   geometrical centroid
  # of the transformed matrices c=sum(Xm* Hm) /M
  # by minimizing the sum(norm(Xm* Hm- C)^2)
  # (Xm* Hm) is considered as transformed matrices
  # the solution is achived by applying an iterative algorithm
  # first:  define the inical centroid C
  # second: finding the similarity taransformation matrices (Hm)
  # finaly: after the calculation of (Xm* Hm), iterative updating of the 
  # centroid C 
  # until the global convergence 
  
  # define the inical centroid C
  Cold=tab[[1]]
  threshold = 1e-10
  fit = 10;
  previousfit = Itot;
  
  
  # iterative algorithm
  while (fit > threshold) {  
    C = matrix(0,P,nmax)
    for (m in 1:M) {
      Hm = gpa(tab[[m]], Cold)$H
      C = C+ tab[[m]] %*% Hm
    }
    
    C <- C/M
    
    def = 0
    for (m in 1:M) {
      d  <- gpa(tab[[m]], C)$d
      def<- def + d
    }
    
    
    #------
    fit = previousfit - def
    previousfit = def    
    Cold <- C
  }   
  
  # finding common loadings #C<-svd(C)$u   
  C <- svd(C %*% t(C))$u[, 1:ncomp] 
  res$loadings.common = C
  rownames(res$loadings.common) = colnames(Data)
  colnames(res$loadings.common) = paste("Dim", 1:ncomp, sep="")
  
  
  # variance of each loading # lambda = t(common loading)*(t(Xm)* Xm) * common loading
  lambda = matrix(0, nrow=M, ncol=ncomp)
  for(m in 1:M){
    lambda[m,] = round(diag(t(C) %*% cov.Group[[m]] %*% C),3)
  }
  res$lambda = lambda
  rownames(res$lambda) = levels(Group)
  colnames(res$lambda) = paste("Dim", 1:ncomp, sep="")
  
  
  #
  exp.var = matrix(0,M,ncomp)
  for(m in 1:M){
    exp.var[m,] = 100 * lambda[m,]/ sum(diag(cov.Group[[m]]))
  }
  res$exp.var = exp.var
  rownames(res$exp.var) = levels(Group)
  colnames(res$exp.var) = paste("Dim", 1:ncomp, sep="")
  
  
  #============================================================================
  if(graph) {plot.mg(res)}
  
  # add class
  class(res) = c("DGPA", "mg")
  return(res)
}


#' @S3method print DGPA
print.DGPA <- function(x, ...)
{
  cat("\nDual Generalized Procrustes Analysis\n")
  cat(rep("-",43), sep="")
  cat("\n$loadings.common   ", "common loadings")
  cat("\n$Data              ", "Data set")
  invisible(x)
}

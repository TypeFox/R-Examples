#'Point estimate of the partition using a modified Binder loss function
#'
#'Get a point estimate of the partition using a modified Binder loss function
#'for Gaussian components
#'
#'@details
#'Note that he current implementation only allows Gaussian components.
#'
#'The modified Binder loss function takes into account the distance between
#'mixture components using #'the  Bhattacharyya distance.
#'
#'@param c a list of vector of length \code{n}. \code{c[[j]][i]} is
#'the cluster allocation of observation \code{i=1...n} at iteration
#'\code{j=1...N}.
#'
#'@param Mu is a list of length \code{n} composed of \code{p x l} matrices.
#'Where \code{l} is the maximum number of components per partition.
#'
#'@param Sigma is list of length \code{n} composed of arrays containing a maximum of
#'\code{l} \code{p x p} covariance matrices.
#'
#'@param lambda is a nonnegative tunning parameter allowing further control over the distance
#'function. Default is 0.
#'
#'@param  a nonnegative constant seen as the unit cost for pairwise misclassification. Default is 1.
#'
#'@param  b nonnegative constant seen as the unit cost for the other kind of pairwise misclassification.
#'Default is 1.
#'
#'@param logposterior vector of logposterior correponding to each
#'partition from \code{c} used to break ties when minimizing the cost function
#'
#'@return a \code{list}:
#'  \itemize{
#'      \item{\code{c_est}:}{ a vector of length \code{n}. Point estimate of the partition}
#'      \item{\code{cost}:}{ a vector of length \code{N}. \code{cost[j]} is the cost
#'      associated to partition \code{c[[j]]}}
#'      \item{\code{similarity}:}{  matrix of size \code{n x n}. Similarity matrix
#'      (see \link{similarityMat})}
#'      \item{\code{opt_ind}:}{ the index of the optimal partition
#'      among the MCMC iterations.}
#'  }
#'
#'
#'@author Chariff Alkhassim
#'
#'
#'@references JW Lau, PJ Green, Bayesian Model-Based Clustering Procedures,
#'\emph{Journal of Computational and Graphical Statistics}, 16(3):526-558, 2007.
#'
#'DA Binder, Bayesian cluster analysis, \emph{Biometrika} 65(1):31-38, 1978.
#'
#'@seealso \code{\link{similarityMat}} \code{\link{similarityMatC}}
#'\code{\link{similarityMat_nocostC}}
#'
#'@export

cluster_est_Mbinder_norm <- function(c, Mu, Sigma, lambda = 0, a = 1, b = a, logposterior){
  # Computes a weighted Bhattacharyya distance between two Gaussians.
  Bhattacharyya <- function(Mu,Sigma,lambda,c){
    theta_grid <- NULL
    uni<-unique(c)
    theta_mat<-matrix(0,max(uni),max(uni))
    for (i in uni)
    {
      theta_grid<- rbind(theta_grid, cbind(i,uni))
    }
    theta_grid<-theta_grid[-which(theta_grid[,1]==theta_grid[,2]),]
    theta_grid<-unique(t(apply(theta_grid,1,FUN=sort)))
    # Number of iterations is the binomial coefficients n=2, k=total number of centers/partition.
    for (i in 1:nrow(theta_grid)){
      mu1<-Mu[,theta_grid[i,1]]
      mu2<-Mu[,theta_grid[i,2]]
      Sigma1<-Sigma[,,theta_grid[i,1]]
      Sigma2<-Sigma[,,theta_grid[i,2]]
      S_bar<-(Sigma1+Sigma2)/2
      if (!det(S_bar)){
        stop("check if each group has its center") # If this occurs then error is in DPMGibbsN
      }
      temp<-((((t(mu1-mu2)%*%solve(S_bar)%*%(mu1-mu2))/8)+
                .5*log(det(S_bar)/sqrt(det(Sigma1)*det(Sigma2)))))
      theta_mat[theta_grid[i,1],theta_grid[i,2]]<-(1-exp(-lambda*temp))
    }
    return(theta_mat)
  }

  len_c<-length(c)

  NuMat_res<-0
  for (i in 1:len_c){
    NuMat_res<-NuMat_res+
      NuMatParC(c[[i]], Bhattacharyya(Mu[[i]],Sigma[[i]],lambda,c[[i]]))$NuMatParC
  }
  NuMat_res<-NuMat_res/len_c
  similarityMat_res<-similarityMat_nocostC(sapply(c, "["))$similarity

  cost<-numeric(len_c)
  for (i in 1:len_c){
    hatc_i_equals_hatc_j<-vclust2mcoclustC(c[[i]])$Coclust
    hatc_i_equals_hatc_j[lower.tri(hatc_i_equals_hatc_j,diag=TRUE)] <- NA
    cost[i]<-sum((b*NuMat_res-a*similarityMat_res)[which(hatc_i_equals_hatc_j==1,arr.ind=TRUE)],na.rm=TRUE)
  }
  opt_ind<-which(cost==min(cost))
  if(length(opt_ind)>1){   # in case of ties
    opt_ind <- opt_ind[which.max(logposterior[opt_ind])]
  }
  c_est<-c[[opt_ind]]
  return(list("opt_ind"=opt_ind,"c_est"=c_est,"cost"=cost))
}

# Below is an implementation in R which doesn't call any C function

#
# cluster_est_MBinderN<-function(c,Mu,Sigma,lambda,a=1,b=1,logposterior)
# {
#   Bhattacharyya<-function(Mu,Sigma,coefs,lambda,c)
#   {
#     theta_grid <- NULL
#     uni<-unique(c)
#     for (i in uni)
#     {
#       theta_grid<- rbind(theta_grid, cbind(i,uni))
#     }
#     theta_grid<-theta_grid[-which(theta_grid[,1]==theta_grid[,2]),]
#     theta_grid<-unique(t(apply(theta_grid,1,FUN=sort)))
#     d<-numeric(nrow(theta_grid))
#     for (i in 1:nrow(theta_grid))
#     {
#       mu1<-Mu[,theta_grid[i,1]]
#       mu2<-Mu[,theta_grid[i,2]]
#       Sigma1<-Sigma[,,theta_grid[i,1]]
#       Sigma2<-Sigma[,,theta_grid[i,2]]
#       S_bar<-(Sigma1+Sigma2)/2
#       if (!det(S_bar)){
#         stop("check if each group has its center") # If this occurs then error is in DPMGibbsN
#       }
#       temp<-((((t(mu1-mu2)%*%solve(S_bar)%*%(mu1-mu2))/8)+
#                 .5*log(det(S_bar)/sqrt(det(Sigma1)*det(Sigma2)))))
#       d[i]<-(1-exp(-lambda*temp))
#     }
#     theta_dist<-cbind(theta_grid,d)
#     dist_coef<-numeric(nrow(coefs))
#
#     # Heavy part of the algorithm . Indeed, each coeficients of the complementary
#     # of a given partition must be assigned its distance weight
#
#     rowmatch <- function(A,B) {
#       f <- function(...) paste(..., sep=":")
#       a <- do.call("f", as.data.frame(A))
#       b <- do.call("f", as.data.frame(B))
#       match(b, a)
#     }
#     index<-rowmatch(theta_dist[,1:2],
#                     t(apply(matrix(c(c[coefs[,1]],c[coefs[,2]]),ncol=2),
#                             1,FUN=sort)))
#     return(theta_dist[index,3])
#   }
#
#   NuMat <- function(c,mu,Sigma,lambda){
#     vclust2mcoclust2 <-function(v,mu,Sigma,lambda){
#       co_clust_dif<-(sapply(v, FUN=function(x){x!=v}))*1
#       co_clust_dif[lower.tri(co_clust_dif,diag=TRUE)] <- NA
#       Ci_dif_Cj<-which(co_clust_dif==1,arr.ind=TRUE)
#       co_clust_dif[Ci_dif_Cj]<-co_clust_dif[Ci_dif_Cj]*Bhattacharyya(mu,Sigma,Ci_dif_Cj,lambda,v)
#       return(co_clust_dif)}
#     temp<-vclust2mcoclust2(c[[1]],mu[[1]],Sigma[[1]],lambda)
#     for (i in 2:length(c))
#     {
#       temp<-temp+vclust2mcoclust2(c[[i]],mu[[i]],Sigma[[i]],lambda)
#     }
#     return(temp/length(c))
#   }
#
#   cost<-numeric(length(c))
#   NuMat_res<-NuMat(c,Mu,Sigma,lambda)
#   similarityMat_res<-similarityMat(c)
#   vclust2mcoclust <- function(v)
#   {
#     m <- sapply(v, FUN=function(x){x==v})
#     #m[lower.tri(m,diag=TRUE)]<-NA
#     return(m*1)
#   }
#   for (i in 1:length(c))
#   {
#     hatc_i_equals_hatc_j<-vclust2mcoclust(c[[i]])
#     hatc_i_equals_hatc_j[lower.tri(hatc_i_equals_hatc_j,diag=TRUE)] <- NA
#     cost[i]<-sum((b*NuMat_res-a*similarityMat_res)[which(hatc_i_equals_hatc_j==1,arr.ind=TRUE)],na.rm=TRUE)
#   }
#
#   opt_ind<-which(cost==min(cost))
#   if(length(opt_ind)>1){   # in case of ties
#     opt_ind <- opt_ind[which.max(logposterior[opt_ind])]
#   }
#   c_est<-c[[opt_ind]]
#   return(list("opt_ind"=opt_ind,"c_est"=c_est,"cost"=cost))
# }




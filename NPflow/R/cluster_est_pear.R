#'Gets a point estimate of the partition using posterior expected adjusted
#'Rand index (PEAR)
#'
#'@param c a list of vector of length \code{n}. \code{c[[j]][i]} is
#'the cluster allocation of observation \code{i=1...n} at iteration
#'\code{j=1...N}.
#'
#'
#'@return a \code{list}:
#'  \itemize{
#'      \item{\code{c_est}:}{ a vector of length \code{n}. Point estimate of the partition}
#'      \item{\code{pear}:}{ a vector of length \code{N}. \code{pear[j]} is the
#'      posterior expected adjusted Rand index associated to partition \code{c[[j]]}}
#'      \item{\code{similarity}:}{  matrix of size \code{n x n}. Similarity matrix
#'      (see \link{similarityMat})}
#'      \item{\code{opt_ind}:}{ the index of the optimal partition
#'      among the MCMC iterations.}
#'  }
#'
#'
#'@author Chariff Alkhassim
#'
#'@export
#'
#'@references A. Fritsch, K. Ickstadt. Improved Criteria for Clustering Based
#' on the Posterior Similarity Matrix, in Bayesian Analysis, Vol.4 : p.367-392
#' (2009)
#'
#'@seealso \code{\link{similarityMat}} \code{\link{similarityMatC}}
#'

cluster_est_pear <- function(c)
{
  cat("Estimating posterior similarity matrix...\n(this may take some time, complexity in O(n^2))\n")
  cmat <- sapply(c, "[")
  tempC <- similarityMatC(cmat)
  cat("DONE!\n")

  similarity <- tempC$similarity

  Nclust<-length(c)
  pear<-rep(NA,Nclust)
  PEARclust<--Inf
  clust<-NULL
  opt_ind<-NULL
  n<-length(c[[1]])
  comb<-choose(n,2)
  indic<-seq_len(n)
  for (i in 1:Nclust)
  {
    s<-split(indic,c[[i]])
    clust_ind<-rbind()
    for (t in 1:length(s))
    {
      len_st<-length(s[[t]])
      if (len_st==1)
      {
        clust_ind<-rbind(clust_ind,rep(s[[t]],2))
      }
      else
      {
        row<-rep(min(s[[t]]),(len_st-1))
        col<-s[[t]][2:len_st]
        clust_ind<-rbind(clust_ind,cbind(row,col))
      }
    }
    quant<-nrow(clust_ind)
    Index<-sum(similarity[clust_ind])
    ExpectedIndex<-Index/comb
    MaximumIndex<-.5*(quant+sum(similarity,na.rm=TRUE))-ExpectedIndex
    PEAR<-(Index-ExpectedIndex)/(MaximumIndex-ExpectedIndex)
    pear[i]<-PEAR
    if (PEAR>PEARclust)
    {
      PEARclust<-PEAR
      clust<-c[[i]]
      opt_ind<-i
    }

    return(list("c_est"=clust, "pear"=pear, "similarity"=similarity,
                "opt_ind"=opt_ind))
  }
}




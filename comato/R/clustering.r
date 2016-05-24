#' Non-randomness of data
#' 
#' \code{Hopkins.index} calculates the Hopkins index that can be used as an indicator of the non-randomness of data prior to clustering.
#' @param data A numeric matrix.
#' @return The Hopkins index as a numeric value 
#' @seealso The index is described in, e.g.: Han, Jiawei; Kamber, Micheline (2010): Data mining. Concepts and techniques. 2nd ed., Amsterdam: Elsevier/Morgan Kaufmann (The Morgan Kaufmann series in data management systems).
#' @examples
#' \dontrun{
#' #Random data generation, 10 dimensions, 500 observations, 2 clusters, 
#' #Multivariate-Bernoulli distributed
#' require("gtools")
#' data = c()
#' p = 0.0
#' for (i in 1:2)
#' {
#' temp = c()
#' for (j in 1:10)
#' temp = cbind(temp, rbinom(250, 1, p+(i-1)*0.5+(0.025*i)*j))  
#' data=rbind(data, temp)
#' }
#' data = data[permute(1:500),]
#' 
#' Hopkins.index(data)
#' }
#'@export
Hopkins.index <- function(data)
{
  n = round(0.05*dim(data)[1])
  
  sample = round(runif(n, 0, dim(data)[1]))
  W = c()
  for (s in sample)
    W = c(W, min(rowSums(abs(data[-s,] - matrix(rep(data[s,], dim(data)[1]-1), nrow=dim(data)[1]-1, byrow=T)))))
  
  points = matrix(round(runif(n*dim(data)[2], min=min(data), max=max(data))), nrow = n, ncol=dim(data)[2])
  U = c()
  for (p in 1:n)
    U = c(U, min(rowSums(abs(data - matrix(rep(points[p,], dim(data)[1]), nrow=dim(data)[1], byrow=T)))))
  
  return(sum(U)/(sum(U) + sum(W)))
  
}


MBM.em <- function(data, M, eps)
{
 
  n = dim(data)[1]    #Number of observations
  l = dim(data)[2]    #Length of an observation
  
  
  theta = matrix(data=runif(M*l), nrow=M, ncol=l) #Random start values
  alpha = rep(1/M, M)                             #All components equally mixed
  prob = matrix(1/M, ncol=M, nrow=n)              #all components equally probable
  
  iter=0
  diff = eps+1
  while (diff > eps)  #Repeat as long as values change significantly
  {   
    #E-Step              
    prob = t(apply(data, 1, FUN=function(val) {
      x = apply(t(theta)*val + (1-t(theta))*(1-val), 2, prod)
      alpha*x/sum(alpha*x)
    }))
    if (M == 1)   
      prob = t(prob)
    
    #M-Step
    new_alpha= colSums(prob)/n
    
    new_theta = (t(apply(prob, 2, FUN=function(val) {    
      colSums(data*val)
    }))/new_alpha)/n

    diff = max(max(abs(alpha-new_alpha)), max(abs(theta-new_theta))) #Compute change in values
  
    alpha = new_alpha   #Prepare next iteration
    theta = new_theta    
    iter = iter+1
  }
  return(list(alpha, theta, prob, iter))
}


#' MBMM clustering
#' 
#' \code{MBM.cluster} calculates a model based clustering using multivariate Bernoulli-mixtures as probabilistic model of the data.
#' The quality of the clustering is judged using the AIC criterion.
#' @param data A numeric matrix.
#' @param min The minimal number of components that is tested.
#' @param max The maximal number of components that is tested.
#' @return A list with 3 elements. The first element is the minimal AIC value for each tested number of components.
#' The second element is a vector of all AIC values. The third is the actual clustering as returned by the EM algorithm using
#' the optimal number of components according to AIC. The element is again a list that contains the mixture coefficients, the actual
#' parameters of the mutlivariate Benroulli distributions, the probability matrix of each observation (i.e. row if  \code{data}) 
#' and component and the number of iterations that the EM algorithm needed to converge.
#' @examples
#' #Random data generation, 100 observations, 5 dimensions, dependencies within the dimensions
#' data = cbind(round(runif(100)), round(runif(100)), round(runif(100)))
#' data = cbind(data, data[,2], 1-data[,3])
#' 
#' #Noisy data:
#' s = round(runif(2, 1, 100))
#' data[s, c(4,5)] = 1 - data[s, c(4,5)]
#' 
#' #MBMM Clustering
#' res = MBM.cluster(data, 1,8)
#'@export 
MBM.cluster <- function(data, min=1, max=10)
{
  AICc = Inf
  em_res = list()
  for (m in min:max)
  {    
    cat("Testing model with ")
    cat(m)
    cat(" components\n")
    k=(m-1 + dim(data)[2]*m)
    if (k+1 > dim(data)[1])
    {
      cat("Not enough samples!")
      break
    }
    res = MBM.em(data, m, 0.0001)
    temp = 0
    for (i in 1:dim(data)[1])
    {
      val=c()      
      for (j in 1:m)
        val = c(val, prod(c(res[[2]][j, which(data[i,] == 1)], 1-res[[2]][j, which(data[i,] == 0)])))
      temp = temp + log(sum(res[[1]]*val))      
    }
    temp = -2*temp + 2*k*dim(data)[1]/(dim(data)[1] - k - 1)
    if (temp < min(AICc))
      em_res = res
    AICc = c(AICc, temp)
  }
  return(list(min(AICc), AICc[-1], em_res))
}


#' Similarity based clustering
#' 
#' \code{PAM.cluster} calculates a clustering using the PAM algorithm (k-medoids). The quality of the clustering is judged using the 
#' G1 index.
#' @param data A numeric matrix.
#' @param min The minimal number of components that is tested. Must be at least 2.
#' @param max The maximal number of components that is tested.
#' @param metric If empty, data will be treated as a distance matrix. Otherwise, the value will be passed to the call of \code{dist} 
#' to compute the distance matrix from \code{data}
#' @return A list with 3 elements. The first element contains the optimal number of components according to the G1 index. 
#' The second element contains a vector of the G1 values. The thrid element contains the clustering itself, 
#' i.e. the return value of \code{PAM}.
#' @examples
#' \dontrun{
#' #Random data generation, 10 dimensions, 500 observations, 2 clusters
#' require("gtools")
#' data = c()
#' p = 0.0
#' for (i in 1:2)
#' {
#' temp = c()
#' for (j in 1:10)
#' temp = cbind(temp, rbinom(250, 1, p+(i-1)*0.5+(0.025*i)*j))  
#' data=rbind(data, temp)
#' }
#' data = data[permute(1:500),]
#' 
#' PAM.cluster(data)
#' }
#'@export 
PAM.cluster <- function(data, min=2, max=10, metric="manhattan")
{
  G1vals = c()
  clustering = NULL
  comps = 0
  for (m in min:max)
  {    
    cat("Testing ")
    cat(m)
    cat(" clusters.\n")
    if (metric != "")
      data = dist(data, method=metric)
    res = pam(data, m, diss=T)  
    g1 = index.G1(data, res$clustering)
    if (length(G1vals) == 0 || g1 > max(G1vals))
    {
      clustering = res
      comps = m
    }   
    G1vals = c(G1vals, g1)
  }
  return(list(comps, G1vals, clustering))
}


#' Clustering maps of a conceptmaps object
#' 
#' \code{clustering} is a convenience function that implements two frequently used ways of clustering conceptmaps directly.
#' The first is clustering using the MBMM algorithm and the concept matrix, the second is clustering using the PAM algorithm and 
#' the graph similarity matrix.
#' @param maps A conceptmaps object.
#' @param method Either "PAM" or "MBMM", indicating which algorithm should be used.
#' @param min The minimal number of components that is tested. For the PAM algorithm, 1 is not allowed.
#' @param max The maximal number of components that is tested.
#' @return The return value of either \code{\link{MBM.cluster}} or \code{\link{PAM.cluster}}, depending on the value of \code{method}.
#' @examples
#' \dontrun{
#' #Assuming that there are concept maps in folder "~/maps"
#' cms = read.folder.tgf("~/maps")
#' 
#' clustering(cms, method="MBMM")
#' }
#'@export 
clustering <- function(maps, method=c("MBMM", "PAM"), min=1, max=10)
{
  if (method == "MBMM") 
  {
    mat = landscape(maps, result="matrix", mode="concept.vector")
    return(MBM.cluster(mat, min, max))
  }
  if (method == "PAM") 
  {
    mat = landscape(maps, result="matrix", mode="graph.sim")    
    return(PAM.cluster(mat, max(min, 2), max, metric=""))
  }
}

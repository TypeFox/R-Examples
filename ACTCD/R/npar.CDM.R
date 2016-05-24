npar.CDM <-
function(Y, Q, cluster.method = c("HACA","Kmeans"), Kmeans.centers = NULL, 
                Kmeans.itermax = 10, Kmeans.nstart = 1, HACA.link = 
                c("complete", "ward", "single","average", "mcquitty", 
                "median", "centroid"),HACA.cut = NULL,label.method = 
                c("2b","2a","1","3"),perm=NULL)
{
  s1 <- Sys.time()
  
  cluster.method <- match.arg(cluster.method)
  HACA.link <- match.arg(HACA.link)
  label.method <- match.arg(label.method)
  
  #-----------------------Input check--------------------------#
  input.check(Y = Y,Q = Q,cluster.method = cluster.method,label.method = label.method,perm=perm)
  
  #-----------------------Basic variables----------------------#
  #N:number of examinees
  #J:number of items
  #K:number of attributes
  #M:number of ideal attribute patterns, which is equal to 2^K
  N <- dim(Y)[1]
  J <- dim(Y)[2]
  K <- dim(Q)[2]
  M <- 2^K
  #------------------------------------------------------------#
  
  
  #A:alpha matrix
  #E:eta matrix
  A <- alpha (K)
  E <- eta (K,J,Q)
  
   #====================cluster analysis================================#
  cd.cluster.object <- cd.cluster(Y, Q, cluster.method, Kmeans.centers, 
                                Kmeans.itermax, Kmeans.nstart, HACA.link, 
                                HACA.cut)
  
  
  #===============================Labelling=============================#
  labels <- labeling(Y, Q, cd.cluster.object, method = label.method,perm=perm)
  s2 <- Sys.time()
  output <- list(att.pattern = labels$att.pattern, att.class = labels$att.class,
                 att.dist = labels$att.dist, alpha = A, eta = E, 
                 cluster.size = cd.cluster.object$size, 
                 cluster.class = cd.cluster.object$class, 
                 starting.time = s1, end.time = s2, cluster.method = cluster.method,
                 label.method = label.method)
  class(output) <- "npar.CDM"
  return(output)
}

#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# sugm(): The user interface for sugm()                                            #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Dec 2nd 2013                                                               #
# Version: 1.1.0                                                                   #
#----------------------------------------------------------------------------------#

sugm <- function(data,
                 lambda = NULL,
                 nlambda = NULL,
                 lambda.min.ratio = NULL,
                 rho = NULL,
                 method = "tiger",
                 sym = "or",
                 shrink = NULL,
                 prec = 1e-4,
                 max.ite = 1e4,
                 standardize = FALSE,
                 perturb = TRUE,
                 verbose = TRUE)
{
  if(verbose) {
    cat("High-deimensional Sparse Undirected Graphical Models.\n")
  }
  if(method!="clime" && method!="tiger") {
    cat("\"method\" must be either \"clime\" or \"tiger\" \n")
    return(NULL)
  }
  n = nrow(data)
  d = ncol(data)
  if(method == "tiger" && d<3){
    cat("d>=3 is required for \"tiger\" \n")
    cat("More on help(sugm) \n")
    return(NULL)
  }
  if(method == "clime" && d<2){
    cat("d>=2 is required for \"clime\" \n")
    cat("More on help(sugm) \n")
    return(NULL)
  }
  maxdf = max(n,d)
  est = list()
  est$cov.input = isSymmetric(data)
  if(est$cov.input)
  {
    if(verbose) {
      cat("The input is identified as the covriance matrix.\n")
    }
    if(sum(is.na(data))>0){
      cat("The input has missing values for covariance/correlation input.\n")
      return(NULL)
    }
    if(method=="tiger") {
      cat("The input for \"tiger\" cannot be covriance matrix.\n")
      return(NULL)
    }
    if(standardize){
      S0 = data
      diag.cov=diag(S0)
      diag.cov.invsq = diag(1/sqrt(diag.cov))
      S = diag.cov.invsq%*%S0%*%diag.cov.invsq
    }
    else{
      S0 = data
      S = S0
    }
  }
  if(!est$cov.input)
  {
    X0=data
    if(method=="tiger" && sum(is.na(X0))>0){
      cat("The input for \"tiger\" has missing values.\n")
      return(NULL)
    }
    X1 = X0 - matrix(rep(colMeans(X0),n), nrow=n, byrow=TRUE)
    S0 = crossprod(X1)/(n-1)
    diag.cov=diag(S0)
    diag.cov.invsq = diag(1/sqrt(diag.cov))
    if(method=="tiger") {
      data = X1%*%diag.cov.invsq
      S = S0
    }else{
      if(standardize){
#         S = diag.cov.invsq%*%S0%*%diag.cov.invsq
        S = cor(X0,use="pairwise.complete.obs")
        diag.cov=diag(S)
        diag.cov.invsq = diag(1/sqrt(diag.cov))
      }else{
#         S = S0
        S = cov(X0,use="pairwise.complete.obs")
      }
    }
  }
  
  if(!is.null(lambda)) nlambda = length(lambda)
  if(is.null(lambda))
  {
    if(method == "tiger") {
      if(is.null(nlambda)){
        nlambda = 5
      }
      if(is.null(lambda.min.ratio))
        lambda.min.ratio = 0.4
      lambda.max= pi*sqrt(log(d)/n)
      lambda = seq(lambda.max,lambda.min.ratio*lambda.max,length=nlambda)
    }
    else {
      if(is.null(nlambda))
        nlambda = 5
      if(is.null(lambda.min.ratio))
        lambda.min.ratio = 0.4
      lambda.max.tmp1 = min(max(S-diag(diag(S))),-min(S-diag(diag(S))))
      lambda.max.tmp2 = max(max(S-diag(diag(S))),-min(S-diag(diag(S))))
      if(lambda.max.tmp1==0)
        lambda.max = lambda.max.tmp2
      else
        lambda.max = lambda.max.tmp1
      lambda.min = lambda.min.ratio*lambda.max
      lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
      rm(lambda.max,lambda.min,lambda.min.ratio)
      gc()
    }
  }
  
  est$lambda = lambda
  est$nlambda = nlambda
  
  if(is.null(rho))
    rho = 1
  begt = Sys.time()
  if(method == "clime"){
    if (is.logical(perturb)) {
      if (perturb) { 
        #eigvals = eigen(S, only.values=T)$values
        #perturb = max(max(max(eigvals) - d*min(eigvals), 0)/(d-1), 1/n)
        perturb = 1/sqrt(n)
      } else {
        perturb = 0
      }
    }
    S = S + diag(d)*perturb
    if(method == "clime"){
      if(is.null(shrink)) shrink = 0#1.5
      if(is.null(max.ite)) max.ite=1e4
      re.sugm = sugm.clime.ladm.scr(S, lambda, nlambda, n, d, maxdf, rho, shrink, prec, max.ite, verbose)
      
      if(standardize){
        for(i in 1:nlambda){
          re.sugm$icov1[[i]] = diag.cov.invsq%*%re.sugm$icov1[[i]]%*%diag.cov.invsq
          icov.i = re.sugm$icov1[[i]]
          re.sugm$icov[[i]] = icov.i*(abs(icov.i)<=abs(t(icov.i)))+t(icov.i)*(abs(t(icov.i))<abs(icov.i))
        }
      }else{
        for(i in 1:nlambda){
          icov.i = re.sugm$icov1[[i]]
          re.sugm$icov[[i]] = icov.i*(abs(icov.i)<=abs(t(icov.i)))+t(icov.i)*(abs(t(icov.i))<abs(icov.i))
        }
      }
    }
  }
  
  if(method == "tiger"){
    if(is.null(shrink)) shrink=0
    if(is.null(max.ite)) max.ite=1e4
    
    re.sugm = sugm.tiger.ladm.scr(data, n, d, maxdf, rho, lambda, shrink, prec, max.ite, verbose)
    
    for(i in 1:nlambda){
      re.sugm$icov1[[i]] = diag.cov.invsq%*%re.sugm$icov1[[i]]%*%diag.cov.invsq
      icov.i = re.sugm$icov1[[i]]
      re.sugm$icov[[i]] = icov.i*(abs(icov.i)<=abs(t(icov.i)))+t(icov.i)*(abs(t(icov.i))<abs(icov.i))
    }
  }
  est$ite = re.sugm$ite
  runt = Sys.time()-begt
  
  for(j in 1:d) {
    if(re.sugm$col.cnz[j+1]>re.sugm$col.cnz[j])
    {
      idx.tmp = (re.sugm$col.cnz[j]+1):re.sugm$col.cnz[j+1]
      ord = order(re.sugm$row.idx[idx.tmp])
      re.sugm$row.idx[idx.tmp] = re.sugm$row.idx[ord + re.sugm$col.cnz[j]]
      re.sugm$x[idx.tmp] = re.sugm$x[ord + re.sugm$col.cnz[j]]
    }
  }
  G = new("dgCMatrix", Dim = as.integer(c(d*nlambda,d)), x = as.vector(re.sugm$x[1:re.sugm$col.cnz[d+1]]),
          p = as.integer(re.sugm$col.cnz), i = as.integer(re.sugm$row.idx[1:re.sugm$col.cnz[d+1]]))
  
  est$beta = list()
  est$path = list()
  est$df = matrix(0,d,nlambda)  
  est$sparsity = rep(0,nlambda)  
  for(i in 1:nlambda) {
    est$beta[[i]] = G[((i-1)*d+1):(i*d),]
    est$path[[i]] = abs(est$beta[[i]])
    est$df[,i] = apply(sign(est$path[[i]]),2,sum)
    
    if(sym == "or")
      est$path[[i]] = sign(est$path[[i]] + t(est$path[[i]]))
    if(sym == "and")
      est$path[[i]] = sign(est$path[[i]] * t(est$path[[i]]))
    est$sparsity[i] = sum(est$path[[i]])/d/(d-1)
  }
  rm(G)
  est$runtime = runt
  est$icov = re.sugm$icov
  est$icov1 = re.sugm$icov1
  est$sigma2 = S
  est$sigma = S0
  est$data = data
  est$method = method
  est$sym = sym
  est$verbose = verbose
  est$standardize = standardize
  est$perturb = perturb
  class(est) = "sugm"
  gc()
  return(est)
}

print.sugm <- function(x, ...)
{  
  cat("\n sugm options summary: \n")
  cat(x$nlambda, " lambdas used:\n")
  print(signif(x$lambda,digits=3))
  cat("Method=", x$method, "\n")
  cat("Path length:",x$nlambda,"\n")
  cat("Graph dimension:",ncol(x$data),"\n")
  cat("Sparsity level:",min(x$sparsity),"----->",max(x$sparsity),"\n")
}

plot.sugm = function(x, align = FALSE, ...){
  gcinfo(FALSE)
  
  if(x$nlambda == 1)	par(mfrow = c(1, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
  if(x$nlambda == 2)	par(mfrow = c(1, 3), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
  if(x$nlambda >= 3)	par(mfrow = c(1, 4), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
  
  if(x$nlambda <= 3)	z.final = 1:x$nlambda
  
  if(x$nlambda >=4){
    z.max = max(x$sparsity)
    z.min = min(x$sparsity)
    z = z.max - z.min
    z.unique = unique(c(which(x$sparsity>=(z.min + 0.03*z))[1],which(x$sparsity>=(z.min + 0.07*z))[1],which(x$sparsity>=(z.min + 0.15*z))[1]))
    
    
    if(length(z.unique) == 1){
      if(z.unique<(x$nlambda-1))	z.final = c(z.unique,z.unique+1,z.unique+2)
      if(z.unique==(x$nlambda-1)) z.final = c(z.unique-1,z.unique,z.unique+1)
      if(z.unique==x$nlambda) 	z.final = c(z.unique-2,z.unique-1,z.unique)
    }
    
    if(length(z.unique) == 2){
      if(diff(z.unique)==1){
        if(z.unique[2]<x$nlambda) z.final = c(z.unique,z.unique[2]+1) 
        if(z.unique[2]==x$nlambda) z.final = c(z.unique[1]-1,z.unique)
      }
      if(diff(z.unique)>1) z.final = c(z.unique[1],z.unique[1]+1,z.unique[2])
    }
    
    if(length(z.unique) == 3) z.final = z.unique
    
    rm(z.max,z.min,z,z.unique)
    gc()
    
  }
  plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Sparsity vs. Regularization")
  
  lines(x$lambda[z.final],x$sparsity[z.final],type = "p")
  
  if(align){
    layout.grid = layout.fruchterman.reingold(graph.adjacency(as.matrix(x$path[[z.final[length(z.final)]]]), mode="undirected", diag=FALSE))
    for(i in z.final){
      g = graph.adjacency(as.matrix(x$path[[i]]), mode="undirected", diag=FALSE)
      plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
      rm(g)
      gc()
    }
    rm(layout.grid)
  }
  if(!align){
    for(i in z.final){
      g = graph.adjacency(as.matrix(x$path[[i]]), mode="undirected", diag=FALSE)
      layout.grid = layout.fruchterman.reingold(g)
      plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
      rm(g,layout.grid)
      gc()
    }
  }
  if(align) cat("Three plotted graphs are aligned according to the third graph\n")
}

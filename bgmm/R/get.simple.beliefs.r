get.simple.beliefs <- function(id, k=length(unique(id)), b.min=0.025) {
  kk = length(unique(id))
  others = b.min
  alpha = 1-(k-1)* b.min
  bel = matrix(others, length(id), k)
  colnames(bel) = if(k==kk) unique(id) else c(unique(id), paste("Class", (kk+1):k, sep="."))
  for (i in seq_along(id)) 
     bel[i,id[i]] = alpha
  bel
}

# only for 1d data
init.model.params.HC <- function(dat, k) {
  mu = matrix(sapply(1:k, function(i) mean(dat[dat < quantile(dat, i/k) & dat > quantile(dat, (i - 1)/k)])), k, 1)
  cvar = array(sapply(1:k, function(i) var(dat[dat < quantile(dat, i/k) & dat > quantile(dat, (i - 1)/k)])), c(k,1,1))
  pi = rep(1/k, k)
#  mu = matrix(tapply(sort(dat),
#          ceiling(seq_along(dat)*k/length(dat)),
#          mean), k, 1)
#  cvar = array(tapply(sort(dat),
#          ceiling(seq_along(dat)*k/length(dat)),
#          var), c(k,1,1))
#  pi = rep(1/k, k)
  list(pi = pi, mu = mu, cvar = cvar)
}

init.model.params.knowns <- function(knowns, class, k, d) {
 mu = matrix(0,k,d)
 pro = numeric(k)
 cvar = array(0, c(k, d, d))
 labs = unique(class)
 for (i in seq_along(labs)) {
     mu[i,] = unlist(colMeans(knowns[class==labs[i],,drop=FALSE]))
     cvar[i,,] = cov(knowns[class==labs[i],,drop=FALSE])
     pro[i] = mean(class==labs[i])
  } 
  list(pi = pro, mu = mu, cvar = cvar)
}

init.model.params <- function(X=NULL, knowns=NULL, class=NULL, k=length(unique(class)), method = "all", B=P, P=NULL) {
 d=1
 if (is.null(dim(X)) & !is.null(X))                 {X = as.matrix(X); d=1}
 if (is.null(dim(knowns)) & !is.null(knowns))       {knowns = as.matrix(knowns); d = 1}
 if (!is.null(knowns)) d = ncol(knowns)
 if (!is.null(X)) d = ncol(X)
 if (is.null(class) & !is.null(B))             class = map(B)
 if (is.null(k))   length(unique(class))
 
 if (class(class) == "factor")   class = as.numeric(class)
 if (method == "knowns") {
     model.params = init.model.params.knowns(knowns, class, k, d)
	 # if there is not enough labeled cases just cast an error
	 nmis = sum(is.na(model.params$mu))
	 if (nmis > 0) 		stop("If method='knowns' then labeled cases for all components should be given.")
 }
 if (method == "all") {
    kX = rbind(X, knowns)
    mu = matrix(0, k, d)
    pro = numeric(k)
    cvar = array(0, c(k, d, d))
    if (d == 1) {
      model.params = init.model.params.HC(c(knowns,X), k)
    } else if (k == 1) {
      mu[1,] = unlist(colMeans(kX))
      cvar[1,,] = cov(kX)
      pro = 1
      model.params = list(pi = pro, mu = mu, cvar = cvar)
    } else {
      kres = kmeans(kX, k, nstart=10) 
      labs = 1:k
      for (i in 1:k) {
         mu[i,] = unlist(colMeans(kX[kres$cluster==labs[i],,drop=F]))
         cvar[i,,] = cov(kX[kres$cluster==labs[i],,drop=F])
         pro[i] = mean(kres$cluster==labs[i])
      } 
      model.params = list(pi = pro, mu = mu, cvar = cvar)
    }
    if (!is.null(knowns) & !is.null(class)) {
       clParams     = clusterAssigmentMeans(knowns, class, model.params$mu)
       model.params = permute.model.params(model.params, clParams$li)
    }
 }
 
 model.params
}


permute.model.params <-function(model.params, li) {
  new.model.params = model.params
  for (i in 1:nrow(li)) {
      new.model.params$pi[li[i,2]]    = model.params$pi[li[i,1]]
      new.model.params$mu[li[i,2],]   = model.params$mu[li[i,1],]
      new.model.params$cvar[li[i,2],,]= model.params$cvar[li[i,1],,]
  }
  new.model.params
}

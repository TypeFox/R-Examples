setClass(
         'browntree',
         contains='ouchtree',
         representation=representation(
           call='call',
           nchar='integer',
           data='list',
           theta='list',
           sigma='numeric',
           loglik='numeric'
           )
         )

brown <- function (data, tree) {
  
  if (!is(tree,'ouchtree'))
    stop(sQuote("tree")," must be an object of class ",sQuote("ouchtree"))
  
  if (is.data.frame(data)) {
    nm <- rownames(data)
    data <- lapply(as.list(data),function(x){names(x)<-nm;x})
  }
  if (is.numeric(data)) {
    nm <- deparse(substitute(data))[1]
    data <- list(data)
    names(data) <- nm
  }
  if (is.list(data)) {
    if (
        any(sapply(data,class)!='numeric') ||
        any(sapply(data,length)!=tree@nnodes)
        )
      stop(sQuote("data")," vector(s) must be numeric, with one entry per node of the tree")
    if (any(sapply(data,function(x)(is.null(names(x)))||(!setequal(names(x),tree@nodes)))))
      stop(sQuote("data"), " vector names (or data-frame row names) must match node names of ", sQuote("tree"))
    for (xx in data) {
      no.dats <- which(is.na(xx[tree@nodes[tree@term]]))
      if (length(no.dats)>0)
        stop("missing data on terminal node(s): ",paste(tree@nodes[tree@term[no.dats]],collapse=','))
    }
  } else
  stop(sQuote("data")," must be either a single numeric data set or a list of numeric data sets")

  nterm <- tree@nterm
  nchar <- length(data)
  if (is.null(names(data))) names(data) <- paste('char',1:nchar,sep='')

  data <- lapply(data,function(x)x[tree@nodes])
  
  dat <- lapply(data,function(y)y[tree@term])

  w <- matrix(data=1,nrow=nterm,ncol=1)
  b <- tree@branch.times
  
  sols <- lapply(dat,function(x)glssoln(w,x,b))
  theta <- lapply(sols,function(x)x$coeff) # optima
  e <- sapply(sols,function(x)x$residuals) # residuals
  q <- t(e)%*%solve(b,e)
  v <- q/nterm
  dev <- nchar*nterm*(1+log(2*pi))+nchar*log(det(b))+nterm*log(det(v))
  sigma <- t(chol(v))
  
  new(
      'browntree',
      as(tree,'ouchtree'),
      call=match.call(),
      nchar=nchar,
      data=as.list(data),
      theta=theta,
      sigma=sigma[lower.tri(sigma,diag=T)],
      loglik=-0.5*dev
      )
}

brown.deviate <- function(n = 1, object) {
  sigma <- matrix(0,object@nchar,object@nchar)
  sigma[lower.tri(sigma,diag=T)] <- object@sigma
  v <- kronecker(sigma%*%t(sigma),object@branch.times)
  m <- sapply(object@theta,rep,object@nterm)
  X <- array(
             data=NA,
             dim=c(object@nnodes,object@nchar,n),
             dimnames=list(
               object@nodes,
               names(object@data),
               paste('rep',seq(n),sep='.')
               )
             )
  X[object@term,,] <- array(
                            rmvnorm(n=n,mean=m,var=v),
                            dim=c(object@nterm,object@nchar,n)
                            )
  apply(X,3,as.data.frame)
}

setMethod(
          'simulate',
          'browntree',
          function (object, nsim = 1, seed = NULL, ...) {
            if (!is.null(seed)) {
              if (!exists('.Random.seed',envir=.GlobalEnv)) runif(1)
              save.seed <- get('.Random.seed',envir=.GlobalEnv)
              set.seed(seed)
            }
            X <- brown.deviate(n=nsim,object)
            if (!is.null(seed)) {
              assign('.Random.seed',save.seed,envir=.GlobalEnv)
            }
            X
          }
          )

setMethod(
          'update',
          'browntree',
          function (object, data) {
            if (missing(data)) data <- object@data
            brown(
                  data=data,
                  tree=object
                  )
          }
          )

bootstrap <- function (object, nboot = 200, seed = NULL, ...) {
  stop("function ",sQuote("bootstrap")," is undefined for objects of class ",sQuote(class(object)))
}
setGeneric('bootstrap')  

setMethod(
          "bootstrap",
          "browntree",
          function (object, nboot = 200, seed = NULL, ...) {
            simdata <- simulate(object,nsim=nboot,seed=seed)
            results <- vector(mode='list',length=nboot)
            toshow <- c("sigma.squared","theta","loglik","aic","aic.c","sic","dof")
            for (b in seq_len(nboot)) {
              results[[b]] <- summary(update(object,data=simdata[[b]],...))
            }
            as.data.frame(t(sapply(results,function(x)unlist(x[toshow]))))
          }
          )

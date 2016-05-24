.ergmm.add.fixed<-function(model, X, mean, var, coef.names=NULL, where=c("append","prepend")){
  where <- match.arg(where)
  
  if(length(dim(X))==2) dim(X) <- c(dim(X),1)
  if(!is.null(coef.names)) dimnames(X) <- list(NULL, NULL, coef.names)

  p <- dim(X)[3]
  model[["p"]]<-model[["p"]]+p
  model[["coef.names"]]<-switch(where,
                                append = c(model[["coef.names"]], dimnames(X)[[3]]),
                                prepend = c(dimnames(X)[[3]], model[["coef.names"]])
                                )
                                
  model[["prior"]][["beta.mean"]]<-switch(where,
                                          append = c(model[["prior"]][["beta.mean"]], rep(mean,length.out=p)),
                                          prepend = c(rep(mean,length.out=p), model[["prior"]][["beta.mean"]])
                                          )
  
  model[["prior"]][["beta.var"]]<-switch(where,
                                         append = c(model[["prior"]][["beta.var"]], rep(var,length.out=p)),
                                         prepend = c(rep(var,length.out=p), model[["prior"]][["beta.var"]])
                                         )
                                         

  Yg <- model[["Yg"]]


  Xtmp <- list()
  for(i in seq_len(p)){
    xm <- seldrop(X[,,i,drop=FALSE],3)
    # If the network is undirected, symmetrize:
    if(!is.directed(Yg)) xm <- xm + t(xm)
    
    # If the network is bipartite and the matrix has b1*b2 dimensions,
    # it needs to be augment to (b1+b2)*(b1+b2):
    if(is.bipartite(Yg)
       && all(dim(xm)==c(Yg%n%"bipartite", network.size(Yg)-Yg%n%"bipartite")))
      xm <- bipartite.augment(xm)
    
    Xtmp <- c(Xtmp,list(xm))
  }

  model[["X"]]<-switch(where,
                       append = c(model[["X"]],Xtmp),
                       prepend = c(Xtmp, model[["X"]])
                       )

  model
}

.import.ergm.term <-function(model, term.index, term.name, ..., mean=0, var=9){
  Yg<-model[["Yg"]]
  f <- ~Yg
  
  # this is a stupid hack to pick out one term from the formula; is there a better way?
  f[[3]] <- as.list(attr(terms(model$formula), 'variables'))[[term.index+2]]
  #f[[3]] <- as.call(c(list(term.name), ...))
  if(!is.dyad.independent(f)) warning("Term `", term.name, "` induces dyadic dependence. Likelihood will be effectively replaced by pseudolikelihood.", call.=FALSE)
  if(has.loops(Yg)) warning("Imported ergm term `", term.name, "` will set its dyadic covariate for self-loops, X[i,i,k], to 0. Use `loopfactor` and `loopcov` to model self-loops.", call.=FALSE)
  
  X <- ergmMPLE(f, output="array")$predictor
  
  X[c(is.na(X))] <- 0

  .ergmm.add.fixed(model, X, mean, var)
}

InitErgmm.Intercept<-InitErgmm.intercept<-InitErgmm.1<-function(model, mean=0, var=9){
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 1:4))
    stop(paste("`edges` model term expected between 0 and 2 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)

  .ergmm.add.fixed(model,
                   matrix(1,network.size(model[["Yg"]]),network.size(model[["Yg"]])),
                   mean, var,
                   "(Intercept)")
}

InitErgmm.loops<-function (model, mean=0, var=9){
  if(!has.loops(model[["Yg"]]))
    stop("Self-loop term is  meaningless in a network without self-loops", call.=FALSE)

  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 1:3))
    stop(paste("`loops` model term expected between 0 and 2 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)

  .ergmm.add.fixed(model, diag(1,network.size(model[["Yg"]]),network.size(model[["Yg"]])), mean, var, "loops")
}

InitErgmm.loopcov <- function (model, attrname, mean=0, var=9){
  if(!has.loops(model[["Yg"]]))
    stop("Self-loop covariates are meaningless in a network without self-loops", call.=FALSE)

  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:4))
    stop(paste("loopcov() model term expected between 1 and 3 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)
  n<-network.size(model[["Yg"]])
  x<-model[["Yg"]] %v% attrname

  xm<-diag(x,n,n)
  cn<-paste("loopcov",attrname,sep=".")

  .ergmm.add.fixed(model, xm, mean, var, cn)
}

InitErgmm.loopfactor <- function (model, attrname, base=1, mean=0, var=9){
  if(!has.loops(model[["Yg"]]))
    stop("Self-loop covariates are meaningless in a network without self-loops", call.=FALSE)

  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:4))
    stop(paste("loopfactor() model term expected between 1 and 3 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)

  n<-network.size(model[["Yg"]])

  if(length(attrname)==1)
    x<-model[["Yg"]] %v% attrname
  else
    do.call(paste, c(lapply(attrname, function(a) get.vertex.attribute(model[["Yg"]], a)), sep="."))

  ls<-sort(unique(x))
  if(NVL(base,0)!=0){
    ls <- ls[-base]
    if(length(ls)==0){
      warning("`loopfactor` term deleted because contributes no statistics.")
      return(model)
    }
  }
  
  mean<-rep(mean,length.out=length(ls))
  var<-rep(var,length.out=length(ls))
  for(li in seq_along(ls)){
    l<-ls[[li]]
    xm<-diag(x==l,n,n)
    cn<-paste('loopfactor',paste(attrname,collapse="."),l,sep=".")
    model <- .ergmm.add.fixed(model, xm, mean[li], var[li], cn)
  }
  model
}

InitErgmm.latentcov<-function (model, x, attrname=NULL,
                               mean=0, var=9) 
{
  if(!has.loops(model[["Yg"]])) warning("Term `latentcov` is deprecated for networks without self-loops. Use `edgecov` from package `ergm` instead.")
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
    stop(paste("latentcov() model term expected between 1 and 4 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)
  #Coerce x to an adjacency matrix
  if(is.network(x)){
    xm<-as.matrix.network(x,matrix.type="adjacency",attrname)
    cn<-if(!is.null(attrname)) attrname else paste("network",length(model[["X"]])+1)
  }else if(is.character(x)){
    xm<-as.matrix.network(model[["Yg"]],matrix.type="adjacency",x)
    cn<-x
  }else{
    xm<-as.matrix(x)
    cn<-if(!is.null(attrname)) attrname else 
				paste("latentcov", as.character(sys.call(0)[[3]]),	sep = ".")
  }

  .ergmm.add.fixed(model, xm, mean, var, cn)
}

InitErgmm.sendercov<-function (model, attrname, force.factor=FALSE, mean=0, var=9) 
{
  if(!has.loops(model[["Yg"]])) warning("Term `sendercov` is deprecated for networks without self-loops. Use `nodeocov`, `nodecov`, `nodeofactor`, or `nodefactor` from package `ergm` instead.")

  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
    stop(paste("sendercov() model term expected between 1 and 4 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)
  if (!is.directed(model[["Yg"]]))
    stop("Sender covariates are not allowed with an undirected network; use 'socialitycov'", call.=FALSE)

  n<-network.size(model[["Yg"]])
  x<-model[["Yg"]] %v% attrname
  if(is.numeric(x) && ! force.factor){
    # Numeric covariate.
    xm<-matrix(x,n,n,byrow=FALSE)
    cn<-paste("sendercov",attrname,sep=".")
    model <- .ergmm.add.fixed(model,xm,mean,var,cn)
  }else{
    # Factor covariate.
    x<-as.factor(x)
    ls<-levels(x)
    mean<-rep(mean,length.out=length(ls))
    var<-rep(var,length.out=length(ls))
    for(li in 2:length(ls)){
      l<-ls[[li]]
      xm<-matrix(x==l,n,n,byrow=FALSE)
      cn<-paste('sendercov',attrname,l,sep=".")
      model <- .ergmm.add.fixed(model, xm, mean[li], var[li], cn)
    }
  }
  model
}
InitErgmm.receivercov<-function (model, attrname, force.factor=FALSE, mean=0, var=9) 
{
  if(!has.loops(model[["Yg"]])) warning("Term `receivercov` is deprecated for networks without self-loops. Use `nodeicov`, `nodecov`, `nodeifactor`, or `nodefactor` from package `ergm` instead.")

  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
    stop(paste("receivercov() model term expected between 1 and 3 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)
  if (!is.directed(model[["Yg"]]))
    stop("Receiver covariates are not allowed with an undirected network; use 'socialitycov'", call.=FALSE)

  n<-network.size(model[["Yg"]])
  x<-model[["Yg"]] %v% attrname
  if(is.numeric(x) && ! force.factor){
    # Numeric covariate.
    xm<-matrix(x,n,n,byrow=TRUE)
    cn<-paste("receivercov",attrname,sep=".")
    model <- .ergmm.add.fixed(model,xm,mean,var,cn)
  }else{
    # Factor covariate.
    x<-as.factor(x)
    ls<-levels(x)
    mean<-rep(mean,length.out=length(ls))
    var<-rep(var,length.out=length(ls))
    for(li in 2:length(ls)){
      l<-ls[[li]]
      xm<-matrix(x==l,n,n,byrow=TRUE)
      cn<-paste('receivercov',attrname,l,sep=".")
      model <- .ergmm.add.fixed(model, xm, mean[li], var[li], cn)
    }
  }
  model
}

InitErgmm.socialitycov<-function (model, attrname, force.factor=FALSE, mean=0, var=9) 
{
  if(!has.loops(model[["Yg"]])) warning("Term `receivercov` is deprecated for networks without self-loops. Use `nodeicov`, `nodecov`, `nodeifactor`, or `nodefactor` from package `ergm` instead.")

  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
    stop(paste("socialitycov() model term expected between 1 and 4 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)

  n<-network.size(model[["Yg"]])
  x<-model[["Yg"]] %v% attrname
  if(is.numeric(x) && ! force.factor){
    # Numeric covariate.
    xm<-matrix(x,n,n,byrow=FALSE)+matrix(x,n,n,byrow=TRUE)
    cn<-paste("socialitycov",attrname,sep=".")
    model <- .ergmm.add.fixed(model,xm,mean,var,cn)
  }else{
    # Factor covariate.
    x<-as.factor(x)
    ls<-levels(x)
    mean<-rep(mean,length.out=length(ls))
    var<-rep(var,length.out=length(ls))
    for(li in 2:length(ls)){
      l<-ls[[li]]
      xm<-matrix(x==l,n,n,byrow=FALSE)+matrix(x==l,n,n,byrow=TRUE)
      cn<-paste('socialitycov',attrname,l,sep=".")
      model <- .ergmm.add.fixed(model, xm, mean[li], var[li], cn)
    }
  }
  model
}

# ped.set is a list of pedigrees
# Va utiliser l'environnement .mem pour la mémoisation du cluster 
# et celle des pedigree si on tombe sur quelqu'un qui veut pas de
# calcul parallèle
.mem <- new.env()
ClusterEnv <- new.env()

Likelihood <- function(ped.set, modele, theta, n.cores = getOption("mc.cores", 2L), optim.alloc = TRUE, sum.likelihoods = TRUE, PSOCK = Sys.info()["sysname"] != "Linux")
{
  if(n.cores == 1)
  {
    lik <- vector("list", length(ped.set))
    for(i in seq_along(ped.set))
    {
      a <- Elston(ped.set[[i]], modele, theta, .mem);
      lik[[i]] <- log(a$result);
      .mem <- a$mem;
    }
    if(sum.likelihoods)
      return(Reduce("+", lik))
    else
      return(Reduce(cbind,lik))
  }
  key <- paste("cluster",digest(list("cluster", ped.set, modele)),sep="")
  if(exists(key, envir = .mem)) 
  {
    # cat("reusing cluster\n")
    x <- get(key, envir = .mem)
    Likelihood.next.runs(theta, x, sum.likelihoods)
  }
  else
  {
    # cat("building cluster\n")
    if(PSOCK) {
      cl <- makePSOCKcluster(n.cores)
      clusterExport(cl, c("likelihood.0", "likelihood.00", "likelihood.0.vec", "likelihood.00.vec", "Elston", "ClusterEnv"), environment(likelihood.0))
    } else
      cl <- makeForkCluster(n.cores)

    x <- Likelihood.first.run(ped.set, modele, theta, cl, optim.alloc, sum.likelihoods, PSOCK)
    mem.sv(key, list(cl=cl, o=x$o), .mem)
    return(x$likelihood)
  }
}


es.stopCluster <- function(verbose=TRUE)
{
  for(key in ls(pattern="cluster.*", envir = .mem))  
  {
    x <- get(key, envir = .mem)
    if(verbose) cat("stopping one cluster with",length(x$cl),"nodes\n")
    rm(list=key, envir = .mem)
    stopCluster(x$cl)
  }
}

Likelihood.first.run <- function(ped.set, modele, theta, cl, optim.alloc, sum.likelihoods, PSOCK)
{
  n.cores <- length(cl)
  # préparation découpage des données
  n.fam <- length(ped.set)
  X <- vector("list", n.cores)
  Indices <- vector("list", n.cores)
  i1 <- 1
  for(i in seq_along(cl))
  {
    i2 <- round(i*n.fam/n.cores)
    Indices[[i]] <- i1:i2
    i1 <- i2+1
  
    X[[i]] <- ped.set[ Indices[[i]] ]
  }

  # initialise PED et MODELE
  clusterApply(cl, X, function(x) assign( "PED", x, ClusterEnv) )# PED <<- X[[i]] )
  clusterApply(cl, rep(list(modele), n.cores), function(x) assign("MODELE", x, ClusterEnv)) # MODELE <<- x)
  # initialise MEM 
  clusterApply( cl, seq_len(n.cores), function(i) assign("MEMO", replicate( length(get("PED", envir=ClusterEnv)) , new.env()), ClusterEnv)) #  MEMO <<- replicate(length(PED),new.env()))

  if(!optim.alloc)   
  {
    # calcule
    if(sum.likelihoods)
      return( list( likelihood = Reduce("+", clusterApply(cl, rep(list(theta), n.cores), likelihood.0)), o = seq_along(ped.set) ))
    else
      return( list( likelihood = Reduce(cbind, clusterApply(cl, rep(list(theta), n.cores), likelihood.0.vec)),  o = seq_along(ped.set) ))
  }

  # calcule
  if(sum.likelihoods)
    R <- clusterApply(cl, rep(list(theta), n.cores), likelihood.00)
  else
    R <- clusterApply(cl, rep(list(theta), n.cores), likelihood.00.vec)

  time <- Reduce(function(x,y) c(x,y$time), R, NULL)

  # re-réparti
  Indices <- repartition(time, n.cores)
  o <- order( Reduce(c, Indices) )
  Y <- vector("list", n.cores)
  for(i in seq_along(cl))
  {
    X[[i]] <- ped.set[ Indices[[i]] ]
  }
  clusterApply(cl, X, function(x) assign( "PED", x, ClusterEnv) )# PED <<- X[[i]] )
  clusterApply( cl, seq_len(n.cores), function(i) assign("MEMO", replicate( length(get("PED", envir=ClusterEnv)) , new.env()), ClusterEnv))  #MEMO <<- replicate(length(PED),new.env()))
  if(sum.likelihoods)
    return(list( likelihood = Reduce(function(x,y) x + y$likelihood, R, 0), o = o) )


  return(list( likelihood = Reduce(function(x,y) cbind(x, y$likelihood), R,NULL), o = o ));
   
}

Likelihood.next.runs <- function(theta, x, sum.likelihoods)
{
  n.cores <- length(x$cl)
  if(sum.likelihoods)
    return( Reduce("+", clusterApply(x$cl, rep(list(theta), n.cores), likelihood.0)) )
  else
    return( Reduce(cbind, clusterApply(x$cl, rep(list(theta), n.cores), likelihood.0.vec))[,x$o] )
}

repartition <- function(x, n)
{
  Ind <- rep(list(numeric(0)), n)
  S <- rep(0,n)
  I <- order(x, decreasing=TRUE)
  x <- x[I]
  for(i in seq_along(x))
  {
    a <- x[i];
    k <- which.min(S)  # tjs remplir le moins plein...
    S[k] <- S[k] + a
    Ind[[k]] <- c(Ind[[k]], I[i])
  }
  return(Ind)
}


# --------------------------------------------------------------
# Ces fonctions utilisent les variables globales PED MEMO MODELE 
# et renvoient simplement le résultat
# Elles sont destinées à tourner sur les noeuds du cluster seulement.

#MODELE <- modele.pel
likelihood.00 <- function(theta, ped=get("PED", envir = ClusterEnv), modele = get("MODELE", envir = ClusterEnv), memo = get("MEMO", envir = ClusterEnv) ) 
{ 
  T <- numeric(length(ped))
  lik <- vector("list",length(ped))
  for(i in seq_along(ped))
  {
    T[i] <- system.time( {a <- Elston(ped[[i]], modele, theta, memo[[i]]); } )[1]
    gc() # si on ne force pas gc() après chaque calcul ça fausse l'estimation du temps nécessaire...
    lik[[i]] <- log(a$result);
    memo[[i]] <- a$mem; 
    # cat(i," : ", log(a$result), " (", T[i], ")\n", sep='')
  }
  return(list( likelihood=Reduce("+", lik) , time = T ) )
}

likelihood.0 <- function(theta, ped=get("PED", envir = ClusterEnv), modele = get("MODELE", envir = ClusterEnv), memo = get("MEMO", envir = ClusterEnv) ) 
{ 
  lik <- vector("list",length(ped))
  for(i in seq_along(ped))
  {
    a <- Elston(ped[[i]], modele, theta, memo[[i]]); 
    lik[[i]] <- log(a$result);
    memo[[i]] <- a$mem; 
    # cat(i," : ", log(a$result), "\n", sep='')
  }
  return(Reduce("+", lik))
}

likelihood.00.vec <- function(theta, ped=get("PED", envir = ClusterEnv), modele = get("MODELE", envir = ClusterEnv), memo = get("MEMO", envir = ClusterEnv) ) 
{ 
  T <- numeric(length(ped))
  lik <- vector("list",length(ped))
  for(i in seq_along(ped))
  {
    T[i] <- system.time( {a <- Elston(ped[[i]], modele, theta, memo[[i]]); } )[1]
    gc() # si on ne force pas gc() après chaque calcul ça fausse l'estimation du temps nécessaire...
    lik[[i]] <- log(a$result);
    memo[[i]] <- a$mem; 
    # cat(i," : ", log(a$result), " (", T[i], ")\n", sep='')
  }
  return(list( likelihood=Reduce(cbind,lik), time = T ) )
}

likelihood.0.vec <- function(theta, ped=get("PED", envir = ClusterEnv), modele = get("MODELE", envir = ClusterEnv), memo = get("MEMO", envir = ClusterEnv) ) 
{ 
  lik <- vector("list",length(ped))
  for(i in seq_along(ped))
  {
    a <- Elston(ped[[i]], modele, theta, memo[[i]]); 
    lik[[i]] <- log(a$result);
    memo[[i]] <- a$mem; 
    # cat(i," : ", log(a$result), "\n", sep='')
  }
  return(Reduce(cbind,lik))
}



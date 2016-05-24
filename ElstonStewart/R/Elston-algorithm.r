# split pedigree for ES

split.ped <- function(ped, mem = new.env())
{
  key <- digest(list("split.ped",ped))
  if(exists(key, envir = mem)) return(list(result= get(key, envir = mem), mem = mem))

  n.vec <- which(nb.enfants(ped) > 0 & ped$st[,"pere"] != 0);
  ped$st <- ped$st[, c('id', 'pere', 'mere', 'sexe') ]

  id.pivots <- matrix(ncol=length(n.vec), nrow=2)
  n.pivots  <- matrix(ncol=length(n.vec), nrow=2)
   
  new.id <- max(ped$st[,"id"])
  new.n  <- dim(ped$st)[1]
  for(i in seq_along(n.vec))
  {
    n <- n.vec[i]

    new.n  <- new.n  + 1;
    new.id <- new.id + 1

    id.pivots[1,i] <- ped$st[n,"id"] # la matrice id.pivot contient les ids des
    id.pivots[2,i] <- new.id         # pivots et de leur copie.
    n.pivots[1,i]  <- n;
    n.pivots[2,i]  <- new.n;

    ped$st <- rbind(ped$st, ped$st[n,]);
    ped$st[new.n,"id"] <- new.id
    ped$st[new.n,"pere"] <- ped$st[n,"pere"]    # le nouveau prend les parents
    ped$st[new.n,"mere"] <- ped$st[n,"mere"]
    ped$st[n,"pere"] <- ped$st[n,"mere"] <- 0   # le pivot garde les enfants mais n'a plus de parents
 
    ped$geno  <- c(ped$geno,  ped$geno[n])
    ped$pheno <- c(ped$pheno, ped$pheno[n])
  }
  
  # on repere les lignes qui correspondent aux familles nucleaires
  # [en fait les familles 2G]
  w <- rep(TRUE, new.n)
  f2G <- list();
  while(any(w))
  {
    b <- min( ped$st[ w , "id"] )  # on attrape un nouvel individu
    L <- 0;
    while(length(b) != L)  # on recupere sa famille nucleaire
    {
      L <- length(b);
      w.b <- is.element(ped$st[,"id"],b)
      # On ajoute les parents de tous les individus de la liste b
      b <- union(b, ped$st[w.b,"pere"])
      b <- union(b, ped$st[w.b,"mere"])
      b <- setdiff(b, 0)
      # on ajoute les enfants de tous les individus de la liste b
      b <- union(b, ped$st[ is.element(ped$st[,"pere"], b) | is.element(ped$st[,"mere"], b) , "id"]);
    }
    f2G <- c(f2G, list(w.b))
    w <- w & (!w.b); 
  } 

  # on prepare le terrain en donnant, pour chaque famille dans f2G, les id des pivots
  x <- as.vector(id.pivots)
  f2G.pivots <- vector('list', length(f2G))
  for(i in seq_along(f2G)) f2G.pivots[[i]] <- x[is.element(x, ped$st[f2G[[i]],"id"])]

  r <- list(ped = ped, id.pivots = id.pivots, n.pivots = n.pivots, f2G = f2G, f2G.pivots=f2G.pivots)
  mem.sv(key, r, mem)
  return( list(result=r, mem=mem));
}

Elston.on.splitted <- function(splitted, modele, theta, mem=new.env())
{
  key <- digest(list("bidouille",splitted, modele, theta))
  if(exists(key, envir = mem)) return(list(result= get(key, envir = mem), mem = mem))

  if(dim(splitted$n.pivots)[2] == 0) # plus de pivots
  {
    P <- 1;
    for(i in seq_along(splitted$f2G))
    {  
      w <- splitted$f2G[[i]]
      ped <- list( st = splitted$ped$st[w,,drop=FALSE], geno=splitted$ped$geno[w], pheno=splitted$ped$pheno[w] )
      lik.2g <- Likelihood.2g(ped,modele,theta,mem)
      P <- P*lik.2g$result;
      mem <- lik.2g$mem;
    }
    mem.sv(key, P, mem);  
    return( list(result=P, mem=mem) );
  }

  Re <- extract.pivot(splitted, mem);
  extract <- Re$result;
  mem <- Re$mem;

  # enumeration genotypes
  S <- 0;
  for(a in extract$geno.piv)
  {
    numerateur <- modele$proba.g(a,theta)*modele$p.pheno(extract$pheno.piv,a,theta)
    if(all(numerateur == 0)) ## pas la peine de considerer a, combi geno/pheno incompatible
      next;
    numerateur[ numerateur == 0 ] <- 1;
    extract$splitted$ped$geno[extract$n.piv] <- a
    Re <- Recall(extract$splitted, modele, theta, mem)
    P <- Re$result
    mem <- Re$mem
    for(i in seq_along(extract$F2G))
    {
      f2g <- extract$F2G[[i]]$ped
      n.piv.f2g <- extract$F2G[[i]]$n.piv
      f2g$geno[n.piv.f2g] <- a
      lik.2g <- Likelihood.2g(f2g,modele,theta,mem)
      P <- P*lik.2g$result;
      mem <- lik.2g$mem;
    }
    S <- S + P/numerateur
  }
  mem.sv(key, S, mem);
  return( list(result=S, mem=mem) );
}

extract.pivot <- function(splitted, mem = new.env())
{
  key <- digest(list("extract.pivot",splitted))
  if(exists(key, envir = mem)) return(list(result= get(key, envir = mem), mem = mem))

  # ----------------------------------------
  # une heuristique pour le choix de pivot
  # prioritaire : pivots ou le nombre de genotypes a envisager est le plus petit
  # L <- sapply(splitted$ped$geno[splitted$n.pivots[1,]], length)
  # if(any(L == 1))
  # {
  #   k <- which(L == 1)[1]
  # }
  # else # sinon : pivots pris dans la famille ou il y a le moins de pivot
  # {
  #   id.piv <- splitted$f2G.pivots[[which.min(sapply(splitted$f2G.pivots, length))]][1]
  #   k <- which( apply((splitted$id.pivots == id.piv),2,any) )
  # }
  ### k <- which.min(L)
  # ----------------------------------------
  k <- 1

  n.piv.1  <- splitted$n.pivots[1,k]
  n.piv.2  <- splitted$n.pivots[2,k]
  id.piv.1 <- splitted$id.pivots[1,k]
  id.piv.2 <- splitted$id.pivots[2,k]

  ph <- splitted$ped$pheno[[n.piv.1]];
  g <- splitted$ped$geno[[ n.piv.1 ]]

  # supression pivot...
  splitted$n.pivots  <- splitted$n.pivots[,-k,drop=FALSE]   
  splitted$id.pivots <- splitted$id.pivots[,-k,drop=FALSE]   

  F2G <- list()
  w.fix <- rep(FALSE, length(splitted$f2G))  # pour reperer dans f2G les familles qui n'ont plus de pivots
  w <- rep(FALSE, dim(splitted$ped$st)[1])      # pour reperer les individus qui appartiennent a ces familles [pour extraction]
  for(i in seq_along(splitted$f2G))
  {
    a <- splitted$f2G.pivots[[i]] 
    splitted$f2G.pivots[[i]] <- a[ a != id.piv.1 & a != id.piv.2 ];
    if( length( splitted$f2G.pivots[[i]] ) == 0 ) # plus de pivot dans cette famille
    {
      w.fix[i] <- TRUE
      w <- w | splitted$f2G[[i]];
      # extraction de la famille 2G
      ped <- list(st = splitted$ped$st[splitted$f2G[[i]],,drop=FALSE],
                  pheno = splitted$ped$pheno[splitted$f2G[[i]]], geno=splitted$ped$geno[splitted$f2G[[i]]] )
      # ou sont les pivots dans cette famille ?
      n.piv <- which( ped$st[,"id"] == id.piv.1 | ped$st[,"id"] == id.piv.2 )
      F2G <- c( F2G, list(list(ped=ped, n.piv=n.piv) ))
    }
  }

  # on efface de splitted les familles f2G qui ont ete extraites ci-dessus
  if(any(w))
  {  
    splitted$ped <- list( st = splitted$ped$st[!w,,drop=FALSE] , geno = splitted$ped$geno[!w], pheno = splitted$ped$pheno[!w] )
    splitted$f2G <- splitted$f2G[!w.fix];
    for(i in seq_along(splitted$f2G)) 
      splitted$f2G[[i]] <- splitted$f2G[[i]][!w]
    splitted$f2G.pivots <- splitted$f2G.pivots[!w.fix];
    splitted$n.pivots <-  matrix(sapply( splitted$id.pivots, function(i) which(splitted$ped$st[,"id"]==i) ), nrow=2)
  }
  n.piv = which( is.element(splitted$ped$st[,"id"], c(id.piv.1, id.piv.2)))
  r <- list(n.piv = n.piv , pheno.piv = ph, geno.piv = g, splitted = splitted, F2G = F2G)
  mem.sv(key, r, mem);
  return(list(result=r, mem=mem))
}

Elston <- function(ped, modele, theta, mem = new.env())
{
  if(class(ped) != "es.pedigree")
    stop("Argument ped is not of class es.pedigree")
  Re <- split.ped(ped, mem)
  Elston.on.splitted( Re$result, modele, theta, Re$mem)
}


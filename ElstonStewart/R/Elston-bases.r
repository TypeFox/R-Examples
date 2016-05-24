require(kinship2)
require(digest)
require(parallel)


# pour pouvoir ajouter une colonne avec le nombre d'enfants !
nb.enfants <- function(ped)
{
  n <- dim(ped$st)[1]
  e <- numeric(n)
  for(i in 1:n)
  {
    if(ped$st[i,"sexe"] == 1) e[i] <- sum( ped$st[,"pere"] == ped$st[i,"id"] )
    if(ped$st[i,"sexe"] == 2) e[i] <- sum( ped$st[,"mere"] == ped$st[i,"id"] )
  }
  return(e);
}

### modele pour un marqueur di-allélique
### (avec un placeholder pour la penetrance)
trans <- function(of, fa, mo)
{
  if(fa == 0)
  {
    if(mo == 0 & of == 0)
      return(1);
    if(mo == 1 & (of == 0 | of == 1))
      return(.5);
    if(mo == 2 & of == 1)
      return(1);
  }
  if(fa == 1)
  {
    if(mo == 0 & (of == 0 | of == 1))
      return(.5);
    if(mo == 1 & (of == 0 | of == 2))
      return(.25)
    if(mo == 1 & of == 1)
      return(.5)
    if(mo == 2 & (of == 1 | of == 2))
      return(.5)
  }
  if(fa == 2)
  {
    if(mo == 0 & of == 1)
      return(1);
    if(mo == 1 & (of == 1 | of == 2))
      return(.5);
    if(mo == 2 & of == 2)
      return(1);
  }
  return(0)
};

modele.di <- list(
   name = "diallelic",
   proba.g = function(g, theta)
    {
      if(g == 0)
        return(theta$p**2)
      if(g == 1)
        return(2*theta$p*(1-theta$p))
      if(g == 2)
        return((1-theta$p)**2)
      stop("Unknown Genotype");
    },
  trans = trans,
  p.pheno = function(x,g,theta) 1
)


# Memoization transférée dans les fonctions
# saug ça (pour l'instant ?)
mem.sv <- function(k, v, mem) assign(k, v, envir=mem)

Likelihood.2g <- function(ped, modele, theta, mem = new.env())
{
  # memoization 
  key <- digest(list("Likelihood.2g",ped,modele$name,theta))
  if(exists(key, envir = mem)) return(list(result= get(key, envir = mem), mem = mem))

  # pas (plus) d'enfants, que des individus non apparentés
  if(sum(ped$st[,"pere"]>0)==0)
    return( list( result=Likelihood.unrelated(ped, modele, theta), mem=mem) );

  eff <- effeuille.sibship(ped, mem);
  r <- eff$result
  mem <- eff$mem

  S <- 0;
 
  for(b in r$gp)
  {
    for(c in r$gm)
    {
      r$ped.x$geno[r$n.parents.x] <- c(b,c);
      lik.2g <- Recall(r$ped.x, modele, theta, mem);
      pro.x  <- lik.2g$result
      mem <- lik.2g$mem
      if(all(pro.x == 0)) next;
      P <- 1;
      for( n.1 in r$n.sib )
      { 
        S1 <- 0;
        x  <- ped$pheno[[n.1]];
        gf <- ped$geno[[n.1]];
        for(a in gf)
          S1 <- S1 + modele$p.pheno(x,a, theta)*modele$trans(a,b,c)
        P <- P*S1;
        if(all(P == 0)) break;
      } 
      S <- S + pro.x*P;
    }
  }
  mem.sv(key, S, mem)
  return( list(result = S, mem = mem) )
}  

Likelihood.unrelated <- function(ped, modele, theta)
{
  # cas du pedigree vide
  if(dim(ped$st)[1] == 0) return(1);

  P <- 1;
  for(i in 1:dim(ped$st)[1])
  {
    x <- ped$pheno[[i]];
    g <- ped$geno[[i]];
    S <- 0;
    for(a in g)
      S <- S + modele$proba.g(a,theta)*modele$p.pheno(x,a,theta);
    P <- P*S;
    if(all(P == 0)) return(P)
  }
  P;
}

effeuille.sibship <- function(ped, mem = new.env())
{
  # memoization 
  key <- digest(list("effeuille.sibship",ped))
  if(exists(key, envir = mem)) return(list(result= get(key, envir = mem), mem = mem))

  # choix d'une feuille : le premier enfant du pedigree
  # + préparation du pedigree dépouillé de cette feuille
  n   <- which(ped$st[,"pere"]>0)[1]
  p   <- ped$st[n,"pere"];
  m   <- ped$st[n,"mere"];

  n.p <- which(ped$st[,"id"] == p)
  n.m <- which(ped$st[,"id"] == m)
  gp <- ped$geno[[n.p]];
  gm <- ped$geno[[n.m]];

  # on effeuille d'un coup toute la fratrie !!
  n.sib <- which( ped$st[,"pere"] == p & ped$st[,"mere"] == m )

  ped.x <- list( st=ped$st[-n.sib,], geno=ped$geno[-n.sib], pheno=ped$pheno[-n.sib] )
  n.parents.x <- c(which(ped.x$st[,"id"] == p), which(ped.x$st[,"id"] == m))

  r <- list(gp=gp, gm=gm, n.sib=n.sib, ped.x=ped.x, n.parents.x=n.parents.x)
  mem.sv(key, r, mem)
  return( list(result = r, mem = mem) )
}


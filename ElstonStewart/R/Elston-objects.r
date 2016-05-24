
es.pedigree <- function(id, father, mother, sex, pheno, geno, famid)
{
  if(missing(famid)) 
    famid = ""
  ## data check
  n <- length(id)
  if( length(father) != n | length(mother) != n | length(sex) != n | length(pheno) != n | length(geno) != n )
    stop(paste(famid, "id, father, mother, sex, pheno, geno should have same length"))

  if(any(duplicated(id)))
    stop(paste(famid, "duplicated ids"))

  nz.father <- father[father!=0]
  if( length(setdiff(nz.father, id)) != 0 )
     stop(paste(famid, "father ids", setdiff(nz.father, id),"not found in id list"))
  if(any(sex[match(nz.father, id)] != 1))
     stop(paste(famid, "for some individuals indicated as fathers sex is not male (1)"))

  nz.mother <- mother[mother!=0]
  if( length(setdiff(nz.mother, id)) != 0 )
     stop(paste(famid, "mother ids", setdiff(nz.mother, id),"not found in id list"))
  if(any(sex[match(nz.mother, id)] != 2))
     stop(paste(famid, "for some individuals indicated as mothers sex is not female (2)"))

  if( any(xor(mother == 0, father == 0)) )
     stop(paste(famid, "individuals must have both parents in the pedigree, or be founders"))
  ## --

  ped <- list( st = cbind(id = id, pere = father, mere = mother, sexe = sex), pheno = pheno, geno = geno, famid = famid)
  ped$st <- cbind(ped$st, nb.enfants=nb.enfants(ped));
  class(ped) <- "es.pedigree"
  return(ped)
}

plot.es.pedigree <- function(x, ...)
{
  pheno <- if(hasArg("pheno")) list(...)$pheno else x$pheno 
  if(mode(pheno) == "list")
    pheno <- Reduce( function(x,y) c(x,y[1]), pheno, NULL)
  ped.k2 <- pedigree( x$st[,"id"], x$st[,"pere"], x$st[,"mere"], x$st[,"sexe"], pheno )
  plot.pedigree(ped.k2, ...)
}

print.es.pedigree <- function(x, ...)
  cat("An es.pedigree object with", dim(x$st)[1], "individuals\n")


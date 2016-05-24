#==============================================================================
# File: utils.R
#
# Author: Nathan Morris
#         Yeunjoo Song
#  
# Notes: Utility functions.
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Import hapmap data
#------------------------------------------------------------------------------
importHapmapData = function(chr, pop="CEU", ...)
{
  tmpfilename = paste("tmp_hapmap_", chr, "_", pop, ".gz", sep="")

  myurl = paste("http://hapmap.ncbi.nlm.nih.gov/downloads/phasing/2009-02_phaseIII/HapMap3_r2/",
                toupper(pop),
                "/TRIOS/hapmap3_r2_b36_fwd.consensus.qc.poly.chr", chr, "_",
                tolower(pop), ".phased.gz",
                sep = "")

  download.file(url=myurl, destfile=tmpfilename,...)

  zz = gzfile(tmpfilename, "r")
  tmp = read.table(zz, stringsAsFactors=FALSE, header=TRUE)
  tmp$chr = chr
  close(zz)

  cat("Complete... \n")

  return(tmp)
}

#------------------------------------------------------------------------------
# Print analysis progress information,
#  a property name and value with padding dots between them.
#------------------------------------------------------------------------------
.printInfoLine = function(pInfo, value, width, leadSize=8)
{
  value = as.character(value)
  if( length(value) == 0 )
    value = "NA"

  padSize   = width - nchar(value) - nchar(pInfo)
  dotPad    = paste(rep(".", padSize), collapse="")
  leadSpace = paste(rep(" ", leadSize), collapse="")

  cat(leadSpace, pInfo, dotPad, value, "\n")
}

#------------------------------------------------------------------------------
# Remove leading & trailing spaces; Input: vector, Output: vector
#------------------------------------------------------------------------------
.trimSpace = function(mystring)
{
  return(gsub("(^ +)|( +$)", "", mystring))
  #return(sub('[[:space:]]+', '',sub('[[:space:]]+$', '', mystring)))
}

#------------------------------------------------------------------------------
# Check pedigree dataframe
#  Required fields - family, id, father, mother
#  For founders - father & mother = 0
#------------------------------------------------------------------------------
.checkPed = function(peds)
{
  names = names(peds)

  # Is peds a dataframe with family, id, father, mother
  if( is.data.frame(peds) && (all(c("family", "id", "father", "mother") %in% names)) )
    return(0)  # No error
  else
    return(1)  # Error	
}

#------------------------------------------------------------------------------
# Order pedigrees
#------------------------------------------------------------------------------
.orderPedigrees = function(peds, returnFinalOrder=F)
{
  orderByFamily = order(peds$family)
  npeds = peds[orderByFamily,]
  orderWithinFamily = as.list(by(npeds[c("id","father","mother")], npeds$family, orderPed))
  orderWithinFamily = lapply(orderWithinFamily, order)
  orderWithinFamily = unlist(orderWithinFamily, use.names=FALSE)
  finalOrder = by(data.frame(ordBF=orderByFamily, ordWF=orderWithinFamily),
                  npeds$family,
                  function(x)
                  {
                    return(x$ordBF[x$ordWF])
                  })

  # contains the info on how to reorder the pedigrees in DF
  finalOrder = unlist(finalOrder, use.names=FALSE)

  if( returnFinalOrder )
    return(finalOrder)
  else
    return(peds[finalOrder,])
}

#------------------------------------------------------------------------------
# Get the position of parents in pedigrees
#------------------------------------------------------------------------------
.getParentPos = function(peds)
{
  pos = as.list(by(data.frame(absolutePos=1:nrow(peds),peds),
                   peds$family,
                   function(x)
                   {
                     # relative position of the parent within the pedigree
                     posFather = match(x$father, x$id, nomatch = NA)
                     posMother = match(x$mother, x$id, nomatch = NA)

                     # nomatch must be NA instead of 0 for this to work
                     absoluteMother = x$absolutePos[posMother]
                     absoluteFather = x$absolutePos[posFather]

                     #finding kinship requires 0 instead of NA for missing
                     posFather[is.na(posFather)] = 0
                     posMother[is.na(posMother)] = 0

                     return(list(absoluteMother=absoluteMother,
                                 absoluteFather=absoluteFather,
                                 posFather=posFather,
                                 posMother=posMother,
                                 id=as.character(x$id)))
                   }))

  return(pos)
} 

#------------------------------------------------------------------------------
# Find the missing data
#------------------------------------------------------------------------------
.findMissing = function(yk, isX)
{
  if( isX )
    return(complete.cases(yk))

  return(!is.na(yk))
}

#------------------------------------------------------------------------------
# Filter the missing entries in x and y.
#------------------------------------------------------------------------------
.filterMissing = function(yk, mk, isX)
{
  if( !isX )
    return(yk[mk, drop = FALSE])

  return(yk[mk,, drop = FALSE])
}

#------------------------------------------------------------------------------
# Filter the missing entries in vc matrix.
#------------------------------------------------------------------------------
.filterMissingVC = function(vck, mk1, mk2)
{
  return(lapply(vck, function(z) return(z[mk1,mk2, drop = FALSE])))
}

#------------------------------------------------------------------------------
# Get filtered list of data Y, X VC
#------------------------------------------------------------------------------
.getYFilteredData = function(y, x, vc, t)
{
  tmp_y = lapply(y, function(yk) return(yk[,t]))

  missingY = lapply(tmp_y, .findMissing, isX=FALSE)

  yall  = mapply(.filterMissing,   tmp_y,  missingY, FALSE,    SIMPLIFY=FALSE)
  xall  = mapply(.filterMissing,   x,      missingY, TRUE,     SIMPLIFY=FALSE)
  vcall = mapply(.filterMissingVC, vc,     missingY, missingY, SIMPLIFY=FALSE)

  return(list(y=yall,x=xall,vc=vcall, missing=missingY))
}

#------------------------------------------------------------------------------
# Get filtered list of data COV
#------------------------------------------------------------------------------
.getVCFilteredData = function(vc, missingY1, missingY2)
{
  return(mapply(.filterMissingVC, vc, missingY1, missingY2, SIMPLIFY=FALSE))
}

#------------------------------------------------------------------------------
# Calculate delta values (C and V as a vector) given theta values.
#  - equation 5 and 6 
#  - use library(MASS)
#------------------------------------------------------------------------------
.thetaToDelta = function(model, theta)
{
  B  = model@B(theta)
  L  = model@L(theta)
  Gm = model@Gm(theta)
  Gs = model@Gs(theta)
  Z  = lapply(model@Z, function(v) v(theta))
  E  = lapply(model@E, function(v) v(theta))

  return(.Call("computeDelta", B, L, Gs, Gm, Z, E))
}

#------------------------------------------------------------------------------
# A function which calculates the vector of partial derivatives of a function
#  numerically 
#------------------------------------------------------------------------------
.funDeriv = function(fun, theta, step = 1E-6)
{
  f0 = fun(theta)
  len_fun = length(f0)
  d = matrix(0, length(theta), len_fun)

  for( i in 1:length(theta) )
  {
    h     = rep(0, length(theta))
    h[i]  = step
    d[i,] = (fun(theta + h) - f0)/step
  }

  return(d)
}

#------------------------------------------------------------------------------
# A function to Null Spaces of Matrices with given tolerance
# - improvement over Null {MASS}  
#------------------------------------------------------------------------------
.strumNull = function(M, tol = 1e-07)
{
    tmp <- qr(M, tol=tol)
    set <- if(tmp$rank == 0) 1:ncol(M) else  - (1:tmp$rank)
    qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}	

#------------------------------------------------------------------------------
# A function for sampling from a Dirichlet distribution
# - improvement over Null {MASS}  
#------------------------------------------------------------------------------
.rdirichlet = function(n, alpha)
{
  k = length(alpha)
  r = matrix(0, nrow=n, ncol=k)
  for( i in 1:k )
  {
    r[,i] = rgamma(n, alpha[i], 1)
  }
  r = matrix(mapply(function(r, s) {return (r/s)}, r, rowSums(r)), ncol=k)

  return (r)
}

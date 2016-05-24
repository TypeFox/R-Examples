setMethod("initialize", "RVsharingProb", function(.Object, ... ){
  .Object <- callNextMethod()
  .Object
})

# Compute coefficients of proportionality kappa between kinship of the founders and excess sharing

# By Alexandre Bureau

setMethod("ComputeKinshipPropCoef", signature(obj="RVsharingProb"), function(obj)
{
# obj is an object returned by the function RVsharing 

# Extracting list of distance to founders of each final descendant (excluding the top founders)
desfounders=obj@desfounders
iancestors=obj@iancestors
id = obj@id
dad.id = obj@dad.id
mom.id = obj@mom.id

if (length(iancestors)>2) stop ("Not able yet to handle pedigree with more than 2 intermediate ancestors")
N = length(desfounders)-(length(iancestors)-1)

if (length(iancestors)==2)
  {
  # Append the list of distances to founders of the first intermediate ancestor to the list of distances to founders
  # of every final descendant below him, those that don't have the second intermediate ancestor as last founder in their list
  for (i in 1:N)
    if (names(desfounders[[i]])[length(desfounders[[i]])] != iancestors[2])
      {
      # Adding the distance between spouse of intermediate ancestor and final descendant to the distance between final descendant
      # and intermediate ancestor contained in desfounders[[as.character(iancestors[1])]]
      tmp = unlist(desfounders[[as.character(iancestors[1])]]) + unlist(desfounders[[i]][length(desfounders[[i]])]) 
      desfounders[[i]] = c(desfounders[[i]],tmp)
      }    
  }
  
# Vector of founders
fvec = id[dad.id==0]
nf = length(fvec)
kmat = matrix(0,N-1,N-1)

# Loop over all final descendant
for (i1 in 1:(N-1))
  {
  for (j in 1:length(desfounders[[i1]]))
    {
    for (i2 in (i1+1):N)
      {
         for (k in 1:length(desfounders[[i2]]))
          {
          fj = names(desfounders[[i1]])[j]
          fk = names(desfounders[[i2]])[k]
          # Check that the founders ancestral to each final descendant are not the same person
          if (fj != fk)
            {
            kmat[i1,i2-1] = kmat[i1,i2-1] +  1/2^(as.numeric(desfounders[[i2]][k]) + as.numeric(desfounders[[i1]][j])) 
            }
          }
      }
    }
  }
kmat
} )

setMethod("show", signature(object="RVsharingProb"), function(object)
{
cat("Probability subjects",object@carriers,"among",union(object@carriers,setdiff(names(object@desfounders),object@iancestors)),"share a rare variant: ",object@pshare,"\n")
})

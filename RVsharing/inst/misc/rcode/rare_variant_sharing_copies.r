# Computing probability that a rare variant is shared by a set of subjects in a pedigree

# By Alexandre Bureau

# Version 0.1
# 2013/02/08

# Version 0.2
# 2013/02/12
# Addition of a wrapper for pedigree objects

# Change of fdi from a logical vector to a vector of indices, to correct a problem occuring when 
# an intermediate ancestor is before the final descendants in the ID list
# Restriction of currentnonfounders to the active ones when determining whether there is an intermediate ancestor
# Removed currentfounders that are NA when computing how many there are.
# At depth 0, Include previous ancestors of only the final descendants that are active

# Version 0.3
# 2013/04/08

# Adding the names of the founders to the desfounders list of degrees between final descendants and founders
# Moving the inactivation of final descendants below an intermediate ancestor and the removal of the spouse of the intermediate ancestor
# after the update of the List of distance to founders of each final descendant to insure inclusion of the spouse in the list
# Adding desfounders, id, dad.id and mom.id to the list of returned elements for further computations with other functions

library(kinship2)

RVsharing.copies = function(id, dad.id, mom.id,a)
{
# The function RV sharing computes the probability that all final descendants in the pedigree share a rare variant
# given that a rare variant has been detected in any one of these final descendants
# For now, there can only be one lineage of intermediate ancestors with more than one child each
# Multiple mariages can only involve one of the top founders. Intermediate ancestors can have only one spouse
# All final descendants must share a common ancestor or couple of ancestors, otherwise an erroneous response may be obtained
# Correction to the extration of the kth name : names(tmp)[k]

# The function RVsharing.weighted computes the probability that all sequenced subjects inherit the variant from a single ancestor (num)
# and the probability that no sequenced subject inherited the variant from any of the founders given that a single copy of the variant 
# is present in the founders (p0), weighted by the probability that each founder was the only one to introduce the variant.

# Version 0.1 : Only a specified pair of founders are related with specified kinship coefficient
# Version 0.2 : Many founder pairs can have non-zero kinship coefficients specified in a matrix
#               Correction of an error in Version 0.1 that gave probabilities too low

N = length(id)
# Indicator vector of final descendants
fdi = which(!(id%in%dad.id | id%in%mom.id))
nfd = length(fdi)
if (nfd < 2) stop("There are fewer than 2 descendants for which to compute a rare variant sharing probability.")

# Getting the depth of each subject in the pedigree
dv = kindepth(id, dad.id, mom.id)
md = max(dv)
# Number of founders
Nf = sum(dv==0)
if (any(a<Nf | a>2*Nf)) stop("Number of distinct-by-descent gene copies must be between 0 and 1.")


# Collecting the degrees between the sequenced children, the founders and the intermediate ancestors
degvec = currentnonfounders = currentfounders = numeric(nfd)
active = rep(TRUE,nfd) 
# List of distance to founders of each final descendant
desfounders = list()
# List of list of final descendants and intermediate ancestors below each founder below each ancestor
foundersdegreedes = list()
# List of indicators of whether each descendant of a founder is an intermediate ancestor
iancestor.as.descendant = list()
# List of intermediate ancestors
iancestors = numeric(0)
# List of vectors of degrees of relatedness of final descendants below each intermediate ancestor
ancestorsdegreedes = list()
ia = 1

# Initializing the currentnonfounders vector with the final descendants at the deepest level
currentnonfounders[dv[fdi]==md] = id[dv==md]

# Loop from highest to lowest depth
if (md > 1)
{
for (lev in (md-1):1)
  {
  # Incrementing D for final descendants and intermediate ancestors
  degvec[dv[fdi]>lev&active] = degvec[dv[fdi]>lev&active] + 1
  # Listing ancestors at current depth of final descendants and intermediate ancestors
  # Loop over final descendants and intermediate ancestors with depth greater than current depth
  for (i in (1:nfd)[dv[fdi]>lev&active])
    {
    # The currentnonfounders are those from the previous level
    currentdad=dad.id[id == currentnonfounders[i]]
    currentmom=mom.id[id == currentnonfounders[i]]  
    # Identify non-founder among mom and dad
    currentnonfounders[i] = ifelse(dad.id[id == currentdad]==0,currentmom,currentdad)
    # Identify founder among mom and dad
    currentfounders[i] = ifelse(dad.id[id == currentdad]==0,currentdad,currentmom)
    }
  # Adding final descendent at the current level to the currentnonfounders for the next level
  currentnonfounders[dv[fdi]==lev] = id[fdi][dv[fdi]==lev]

  # Checking if there is an intermediate ancestor with more than one descendant at the current level
  tab.currentnonfounders = table(currentnonfounders[currentnonfounders>0&active])
  # If there is more than one, stop because it is not implemented
  if (sum(tab.currentnonfounders>1) > 1) stop ("More than one intermediate ancestor at the same level with two or more descendant")
  # If there is any, there is only one
  if (any(tab.currentnonfounders>1))
    {
    iancestors[ia] = as.numeric(names(tab.currentnonfounders)[tab.currentnonfounders>1])
    # Adding the degrees of final descendants below the current intermediate ancestor to his list
    ancestorsdegreedes[[ia]] = degvec[currentnonfounders==iancestors[ia]&active]
    # Include degrees between spouse and final descendants in list of ancestors  (assumes only one spouse)
    foundersdegreedes[[ia]] = list(degvec[currentnonfounders==iancestors[ia]&active])
    # Setting indicator of whether the descendant is the previous intermediate ancestor
    if (ia>1)
      iancestor.as.descendant[[ia]] = list(id[fdi][currentnonfounders==iancestors[ia]&active] == iancestors[ia-1])
    else iancestor.as.descendant[[ia]] = list(rep(FALSE,length(foundersdegreedes[[ia]][[1]])))
    # Include previous ancestors of final descendants if any
    if (any(as.numeric(names(desfounders)) %in% id[fdi][currentnonfounders==iancestors[ia]&active]))
      {
      ii = 1
      tmp = desfounders[as.numeric(names(desfounders)) %in% id[fdi][currentnonfounders==iancestors[ia]&active]]
      # Loop over final descendants 
      for (k in 1:length(tmp))
        {
        foundersdegreedes[[ia]][(ii+1):(ii+length(tmp[[k]]))] = tmp[[k]]
        # Setting indicator of whether the descendant of all ancestors in tmp[[k]] is the previous intermediate ancestor
        iancestor.as.descendant[[ia]][(ii+1):(ii+length(tmp[[k]]))] = list(ifelse (ia>1, as.numeric(names(tmp)[k]) == iancestors[ia-1], FALSE))
        ii = ii + length(tmp[[k]])
        # Adding spouse of intermediate ancestor to list of founders of current final descendant
        
        }
      }
    }
  # Adding the current founder ancestral to each final descendants to his list of founders
  if (length(currentfounders[currentfounders>0][!is.na(currentfounders)])>0)
    {
    for (i in (1:nfd)[dv[fdi]>lev&active])
      {
      # If there are at least i elements in desfounders
      if (length(desfounders)>=i)
        {
        desfounders[[i]][length(desfounders[[i]])+1] = degvec[i]
        # Keeping the name of the founder
        names(desfounders[[i]])[length(desfounders[[i]])] = currentfounders[i]
        }
      else 
        {
        desfounders[[i]] = list(degvec[i])
        names(desfounders[[i]])[1] = currentfounders[i]
        }
      # Assigning the ID of the subject as name
      names(desfounders)[i] = id[fdi][i]
      }
    }
  # Finishing processing the current intermediate ancestor if there is one
  if (any(tab.currentnonfounders>1))
    {
    # Turning these final descendants to inactive
    active[currentnonfounders==iancestors[ia]] = FALSE
    # Removing spouse(s) of intermediate ancestor from currentfounders 
    # Note: the spouse(s) have the same positions in the currentfounders vector as the 
    # intermediate ancestor in the currentnonfounders vector
    currentfounders = currentfounders[currentnonfounders != iancestors[ia]]
    # Adding the intermediate ancestor to the vector of subjects with a degree
    nfd = nfd + 1
    fdi[nfd] = which(id==iancestors[ia])
    degvec[nfd] = 0
    active[nfd] = TRUE
    # Adding the intermediate ancestor to the vector of currentnonfounders
    currentnonfounders[nfd] = iancestors[ia]
    # Incrementing ia
    ia = ia + 1
    }
  }
}
# Depth 0: there should be at most 2 founders common to all subjects
# We assign one of them as a dummy "intermediate" ancestor
# Incrementing D for final descendants and intermediate ancestors
degvec[active] = degvec[active] + 1
# Listing ancestors at current depth of final descendants and intermediate ancestors
# The currentnonfounders are those from the previous level
currentdads=dad.id[id %in% currentnonfounders[active]]
currentmoms=mom.id[id %in% currentnonfounders[active]]  
    
# If all subjects have the same dad, use him as last ancestor
if (all(currentdads==currentdads[1]))
  { 
  iancestors[ia] = currentdads[1]
  currentfounders = currentmoms
  }
# else if all subjects have the same mom, use her as last ancestor
else
  {
  if (all(currentmoms==currentmoms[1])) 
    {
    iancestors[ia] = currentmoms[1]
    currentfounders = currentdads
    }  
# else there is no common ancestor, and the probability of sharing is 0
  else return (0)
  }
     
    # Adding the degrees of final descendants below the current intermediate ancestor to his list
    ancestorsdegreedes[[ia]] = degvec[active]
    # Include first spouse in list of ancestors
    spousevec = unique(currentfounders)
    foundersdegreedes[[ia]]= list(degvec[currentfounders==spousevec[1]&active])
    names(foundersdegreedes[[ia]]) = spousevec[1]
    # Setting indicator of whether the descendant is the previous intermediate ancestor
    if (ia>1)
      iancestor.as.descendant[[ia]] = list(id[fdi][currentfounders==spousevec[1]&active] == iancestors[ia-1])
    else iancestor.as.descendant[[ia]] = list(rep(FALSE,length(foundersdegreedes[[ia]][[1]])))
    # Add additional spouses if any
    if(length(spousevec)>1)
      {
      for (i in 2:length(spousevec))
        {
        foundersdegreedes[[ia]][[i]] = degvec[currentfounders==spousevec[i]&active]
        names(foundersdegreedes[[ia]][[i]]) = spousevec[i]
        if (ia>1)
          iancestor.as.descendant[[ia]][[i]] = id[fdi][currentfounders==spousevec[i]&active] == iancestors[ia-1]
        else iancestor.as.descendant[[ia]][[i]] = rep(FALSE,length(foundersdegreedes[[ia]][[i]]))
        }
      }  
    # Include previous ancestors of final descendants if any
    if (any(as.numeric(names(desfounders)) %in% id[fdi][active]))
      {
      ii = length(spousevec)
      tmp = desfounders[as.numeric(names(desfounders)) %in% id[fdi][active]]
      # Loop over final descendants 
      for (k in 1:length(tmp))
        {
        foundersdegreedes[[ia]][(ii+1):(ii+length(tmp[[k]]))] = tmp[[k]]
        # Here it takes a twist to get the names to transfer to foundersdegreedes. See if can improve later.
        for (jj in 1:length(tmp[[k]]))
          {
          names(foundersdegreedes[[ia]])[ii+jj] = names(tmp[[k]][jj])
          }
        # Setting indicator of whether the descendant of all ancestors in tmp[[k]] is the previous intermediate ancestor
        iancestor.as.descendant[[ia]][(ii+1):(ii+length(tmp[[k]]))] = list(ifelse (ia>1, as.numeric(names(tmp)[k]) == iancestors[ia-1], FALSE))
        ii = ii + length(tmp[[k]])
        }
      }
      # Adding the current founder couple ancestral to each final descendants to his list of founders
      # This is not required for the sharing probability computation, but is used for kinship estimation
    # print(currentfounders)
    for (i in (1:nfd)[active])
      {
      j = 1
      # If there are at least i elements in desfounders
      if (length(desfounders)>=i)
        {
        desfounders[[i]][length(desfounders[[i]])+(1:2)] = degvec[i]
        # Keeping the name of the founder
        names(desfounders[[i]])[length(desfounders[[i]])-1] = currentfounders[j]
        names(desfounders[[i]])[length(desfounders[[i]])] = iancestors[ia]
        }
      else 
        {
        desfounders[[i]] = rep(degvec[i],2)
        names(desfounders[[i]])[1] = currentfounders[j]
        names(desfounders[[i]])[2] = iancestors[ia]
        }
      # Assigning the ID of the subject as name
      names(desfounders)[i] = id[fdi][i]
      j = j+1
      }

# Computation of numerator
num = 1
for (i in 1:ia)
  num = num * 1/2^sum(ancestorsdegreedes[[i]])
# Computation of probability of being the only one to introduce the variant for top founder or founders
#print(iancestors[ia])

pfd = 2*(a-Nf)*(2*(a-Nf)-1)/(a*Nf*(2*Nf-1)) + (a-Nf)*2*(2*Nf-a)/(a*((a-Nf)*2*(2*Nf-a) + (a-Nf)*(2*(a-Nf)-1) + 2*(2*Nf-a)*(2*Nf-a-1)))
# If there is only one spouse, then a couple of founders can transmit a variant to all final descendents
if (length(spousevec)==1) 
  {
  #print(spousevec)
  #print(relfounders==spousevec)
  pfd = pfd + 2*(a-Nf)*(2*(a-Nf)-1)/(a*Nf*(2*Nf-1)) + (a-Nf)*2*(2*Nf-a)/(a*((a-Nf)*2*(2*Nf-a) + (a-Nf)*(2*(a-Nf)-1) + 2*(2*Nf-a)*(2*Nf-a-1)))   
  }
num = num*pfd
 
# Computation of denominator
# Probability that no variant has been transmitted
p0 = 0
# Probability that no variant has been transmitted from previous intermediate ancestor
pk = 1
for (i in 1:ia)
  {
  for (j in 1:length(foundersdegreedes[[i]]))
  {
  # Weighting by the probability that names(foundersdegreedes[[i]])[j] was the only one to introduce the variant
  #print(names(foundersdegreedes[[i]])[j])
  #print( relfounders==names(foundersdegreedes[[i]])[j])
    p0 = p0 + (2*(a-Nf)*(2*(a-Nf)-1)/(a*Nf*(2*Nf-1)) + (a-Nf)*2*(2*Nf-a)/(a*((a-Nf)*2*(2*Nf-a) + (a-Nf)*(2*(a-Nf)-1) + 2*(2*Nf-a)*(2*Nf-a-1))))*prod((1-1/2^foundersdegreedes[[i]][[j]]) + ifelse(iancestor.as.descendant[[i]][[j]],(1/2^foundersdegreedes[[i]][[j]])*pk,0))
  }
  # Updating the probability for the previous intermediate ancestor, who becomes the current intermediate ancestor
  # For now, intermediate ancestors can have only one spouse, this is why we take the indicators of the first founder attached to him
  if (i<ia) pk = prod((1-1/2^ancestorsdegreedes[[i]]) + ifelse(iancestor.as.descendant[[i]][[1]],1/2^ancestorsdegreedes[[i]]*pk,0))
  }
# At the end, add the probability from the dummy "intermediate" ancestor, weighted by the probability that he was the only one 
# to introduce the variant. He is currently the only one who can have more than one spouse
# Since only one of his spouses can be the parent of the previous intermediate ancestor, sapply returns only one non-zero term.
# The summation returns in fact the value of that single non-zero term
p0 = p0 + (2*(a-Nf)*(2*(a-Nf)-1)/(a*Nf*(2*Nf-1)) + (a-Nf)*2*(2*Nf-a)/(a*((a-Nf)*2*(2*Nf-a) + (a-Nf)*(2*(a-Nf)-1) + 2*(2*Nf-a)*(2*Nf-a-1))))*prod((1-1/2^ancestorsdegreedes[[ia]]) + sum(sapply(iancestor.as.descendant[[ia]][1:length(spousevec)],function(lv,deg,pk) ifelse(lv, (1/2^deg) * pk,0), deg=ancestorsdegreedes[[i]],pk=pk)))

list(num=num,p0=p0,iancestors=iancestors,desfounders=desfounders,id=id,dad.id=dad.id,mom.id=mom.id)
}

# Computing probability of transmission of a rare variant to sets of subjects in a pedigree from a single ancestor, for 
# the setting where introduction by more than one founder is allowed (and computed separately)
# This function applies with method 2 for approximating variant sharing probabilities in presence of unknown 
# relationships between founders

# By Alexandre Bureau

# Based on function RVsharing

# Version 0.2
# 2013/05/15

# Corrects an erroneous formula in version 0.1

library(kinship2)

RVsharing.approx2 = function(id, dad.id, mom.id,PFU)
{
# The function RV sharing computes the probability that all final descendants in the pedigree share a rare variant
# given that a rare variant has been detected in any one of these final descendants
# For now, there can only be one lineage of intermediate ancestors with more than one child each
# Multiple mariages can only involve one of the top founders. Intermediate ancestors can have only one spouse
# All final descendants must share a common ancestor or couple of ancestors, otherwise an erroneous response may be obtained
# Correction to the extration of the kth name : names(tmp)[k]

# The function RVsharing.approx2 computes the probability that all sequenced subjects inherit the variant from a single ancestor (num)
# and the probability that no sequenced subject inherited the variant from any of the founders given that a single copy of the variant 
# is present in the founders (p0), scaled by the probability PFU that each founder was the only one to introduce the variant. 


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

# Collecting the degrees between the sequenced children, the founders and the intermediate ancestors
degvec = numeric(nfd)
currentnonfounders = currentfounders = character(nfd)
active = rep(TRUE,nfd) 
# List of distance to founders of each final descendant
desfounders = list()
# List of list of final descendants and intermediate ancestors below each founder below each ancestor
foundersdegreedes = list()
# List of indicators of whether each descendant of a founder is an intermediate ancestor
iancestor.as.descendant = list()
# List of intermediate ancestors
iancestors = character(0)
# List of vectors of degrees of relatedness of final descendants below each intermediate ancestor
ancestorsdegreedes = list()
lev.ia = list()
ia = lia = 1

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
    currentnonfounders[i] = ifelse(is.na(dad.id[id == currentdad]),currentmom,currentdad)
    # Identify founder among mom and dad
    currentfounders[i] = ifelse(is.na(dad.id[id == currentdad]),currentdad,currentmom)
    }
  # Adding final descendent at the current level to the currentnonfounders for the next level
  currentnonfounders[dv[fdi]==lev] = id[fdi][dv[fdi]==lev]

  # Checking if there is an intermediate ancestor with more than one descendant at the current level
  tab.currentnonfounders = table(currentnonfounders[currentnonfounders>0&active])
  # If there is more than one, stop because it is not implemented
  # if (sum(tab.currentnonfounders>1) > 1) stop ("More than one intermediate ancestor at the same level with two or more descendant")
  # If there is any intermediate ancestor with more than one descendant at the current level
  if (any(tab.currentnonfounders>1))
    {
    lev.ia[[lia]] = names(tab.currentnonfounders)[tab.currentnonfounders>1]
    iancestors = c(iancestors,lev.ia[[lia]])
    for (i in which(iancestors %in% lev.ia[[lia]]))
    {
    # Adding the degrees of final descendants below the current intermediate ancestor to his list
    ancestorsdegreedes[[i]] = degvec[currentnonfounders==iancestors[i]&active]
    # Include degrees between spouse and final descendants in list of ancestors  (assumes only one spouse)
    foundersdegreedes[[i]] = list(degvec[currentnonfounders==iancestors[i]&active])
    # Setting indicator of whether the descendant is the previous intermediate ancestor
    if (ia>1)
      {
      if (length(lev.ia[[lia]])>1) stop ("More than one intermediate ancestor at level ",lev," with two or more descendant")
      iancestor.as.descendant[[i]] = list(id[fdi][currentnonfounders==iancestors[i]&active] == iancestors[ia-1])
      }
    else iancestor.as.descendant[[i]] = list(rep(FALSE,length(foundersdegreedes[[i]][[1]])))
    # Include previous ancestors of final descendants if any
    if (any(names(desfounders) %in% id[fdi][currentnonfounders==iancestors[i]&active]))
      {
      ii = 1
      tmp = desfounders[names(desfounders) %in% id[fdi][currentnonfounders==iancestors[i]&active]]
      # Loop over final descendants 
      for (k in 1:length(tmp))
        {
        foundersdegreedes[[i]][(ii+1):(ii+length(tmp[[k]]))] = tmp[[k]]
        # Setting indicator of whether the descendant of all ancestors in tmp[[k]] is the previous intermediate ancestor
        iancestor.as.descendant[[i]][(ii+1):(ii+length(tmp[[k]]))] = list(ifelse (ia>1, names(tmp)[k] == iancestors[i-1], FALSE))
        ii = ii + length(tmp[[k]])
        # Adding spouse of intermediate ancestor to list of founders of current final descendant
        
        }
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
    for (i in which(iancestors %in% lev.ia[[lia]]))
    {
    # Turning these final descendants to inactive
    active[currentnonfounders==iancestors[i]] = FALSE
    # Removing spouse(s) of intermediate ancestor from currentfounders 
    # Note: the spouse(s) have the same positions in the currentfounders vector as the 
    # intermediate ancestor in the currentnonfounders vector
    currentfounders = currentfounders[currentnonfounders != iancestors[i]]
    # Adding the intermediate ancestor to the vector of subjects with a degree
    if (any(id==iancestors[i]))
    {
    nfd = nfd + 1
    fdi[nfd] = which(id==iancestors[i])
    degvec[nfd] = 0
    active[nfd] = TRUE
    # Adding the intermediate ancestor to the vector of currentnonfounders
    currentnonfounders[nfd] = iancestors[i]
    }
    }
    # Incrementing ia
    ia = ia + length(lev.ia[[lia]])
    lia = lia + 1
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
    # Setting indicator of whether the descendant is the previous intermediate ancestor
    if (ia>1)
      iancestor.as.descendant[[ia]] = list(id[fdi][currentfounders==spousevec[1]&active] %in% lev.ia[[lia-1]])
    else iancestor.as.descendant[[ia]] = list(rep(FALSE,length(foundersdegreedes[[ia]][[1]])))
    # Add additional spouses if any
    # Warning! This is going to work only if all previous intermediate descendents are under the same spouse
    if(length(spousevec)>1)
      {
      for (i in 2:length(spousevec))
        {
        foundersdegreedes[[ia]][[i]] = degvec[currentfounders==spousevec[i]&active]
        if (ia>1)
          iancestor.as.descendant[[ia]][[i]] = id[fdi][currentfounders==spousevec[i]&active] %in% lev.ia[[lia-1]]
        else iancestor.as.descendant[[ia]][[i]] = rep(FALSE,length(foundersdegreedes[[ia]][[i]]))
        }
      }  
    # Include previous ancestors of final descendants if any
    if (any(names(desfounders) %in% id[fdi][active]))
      {
      ii = length(spousevec)
      tmp = desfounders[names(desfounders) %in% id[fdi][active]]
      # Loop over final descendants 
      for (k in 1:length(tmp))
        {
        foundersdegreedes[[ia]][(ii+1):(ii+length(tmp[[k]]))] = tmp[[k]]
        # Setting indicator of whether the descendant of all ancestors in tmp[[k]] is a previous intermediate ancestor
        iancestor.as.descendant[[ia]][(ii+1):(ii+length(tmp[[k]]))] = list(ifelse (ia>1, names(tmp)[k] %in% iancestors[lev.ia[[lia-1]]], FALSE))
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

pfd = PFU
# If there is only one spouse, then a couple of founders can transmit a variant to all final descendents
if (length(spousevec)==1) 
  {
  #print(spousevec)
  #print(relfounders==spousevec)
  pfd = pfd + PFU   
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
    p0 = p0 + PFU*prod((1-1/2^foundersdegreedes[[i]][[j]]) + ifelse(iancestor.as.descendant[[i]][[j]],(1/2^foundersdegreedes[[i]][[j]])*pk,0))
  }
  # Updating the probability for the previous intermediate ancestor, who becomes the current intermediate ancestor
  # For now, intermediate ancestors can have only one spouse, this is why we take the indicators of the first founder attached to him
  if (i<ia) pk = prod((1-1/2^ancestorsdegreedes[[i]]) + ifelse(iancestor.as.descendant[[i]][[1]],1/2^ancestorsdegreedes[[i]]*pk,0))
  }
# At the end, add the probability from the dummy "intermediate" ancestor, weighted by the probability that he was the only one 
# to introduce the variant. He is currently the only one who can have more than one spouse
# Since only one of his spouses can be the parent of the previous intermediate ancestor, sapply returns only one non-zero term.
# The summation returns in fact the value of that single non-zero term
# This remains to be tested with >1 spouse
  tmpf = as.matrix(sapply(iancestor.as.descendant[[ia]][1:length(spousevec)],function(lv,deg,pk) ifelse(lv, (1/2^deg) * pk,0), deg=ancestorsdegreedes[[ia]],pk=pk))  
  p0 = p0 + PFU*prod((1-1/2^ancestorsdegreedes[[ia]]) + apply(tmpf,1,sum))

list(num=num,p0=p0,iancestors=iancestors,desfounders=desfounders,id=id,dad.id=dad.id,mom.id=mom.id)
}

# By Alexandre Bureau

RVsharing.fn = function(id, dad.id, mom.id,carriers)
{

N = length(id)
# vector of indices of final descendants
fdi = which(!(id%in%dad.id | id%in%mom.id))
nfd = length(fdi)
if (nfd < 2) stop("There are fewer than 2 descendants for which to compute a rare variant sharing probability.")
if (!missing(carriers))
{
	missc = setdiff(carriers,id)
	if(length(missc)>0) stop(missc," not in pedigree.")
	# Check all carriers are non-founders
  names(dad.id) = id
  if (any(is.na(dad.id[as.character(carriers)])))
    stop ("Carriers ",carriers[is.na(dad.id[as.character(carriers)])]," are founders. This is not supported by RVsharing.")
}

# Getting the depth of each subject in the pedigree
dv = kindepth(id, dad.id, mom.id)
# Number of founders
Nf = sum(dv==0)

ped_datastruct = function(fdi,dv)
{
nfd = length(fdi)
md = max(dv)
# List of distance to founders of each final descendant
desfounders = list()
# List of list of final descendants and intermediate ancestors below each founder below each ancestor
foundersdegreedes = list()
# List of indicators of whether each descendant of a founder is an intermediate ancestor
iancestor.as.descendant = list()
# List of intermediate ancestors
iancestors = character(0)
# Number of founders below each intermediate ancestor
iancestors.Nf = numeric(0)
# List of vectors of degrees of relatedness of final descendants below each intermediate ancestor
ancestorsdegreedes = list()
# List of intermediate ancestors on the current pedigree degree
lev.ia = list()
# ia: intermediate ancestors index
# lia: reverse level index from 1 at lev = md-1 to md-1 at lev=1
ia = lia = 1

# Collecting the degrees between the sequenced children, the founders and the intermediate ancestors
degvec = numeric(nfd)
currentnonfounders = currentfounders = character(nfd)
active = rep(TRUE,nfd) 

# Initializing the currentnonfounders vector with the final descendants or carriers at the deepest level
currentnonfounders[dv[fdi]==md] = id[fdi][dv[fdi]==md]

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
  # Adding final descendants at the current level to the currentnonfounders for the next level
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
    names(ancestorsdegreedes[[i]]) = id[fdi][currentnonfounders==iancestors[i]&active]
    # Include degrees between spouse and final descendants in list of ancestors  (assumes only one spouse)
    foundersdegreedes[[i]] = list(degvec[currentnonfounders==iancestors[i]&active])
    # Setting indicator of whether the descendant is the previous intermediate ancestor
    if (ia>1)
      {
      if (length(lev.ia[[lia-1]])>1) stop ("More than one subject at level ",lev," with two or more descendants (intermediate ancestors).")
      iancestor.as.descendant[[i]] = list(id[fdi][currentnonfounders==iancestors[i]&active] == iancestors[ia-1])
      # Recording number of founders below intermediate ancestor
      }
    else iancestor.as.descendant[[i]] = list(rep(FALSE,length(foundersdegreedes[[i]][[1]])))
    # Include previous ancestors of final descendants if any
    if (any(names(desfounders) %in% id[fdi][currentnonfounders==iancestors[i]&active]))
      {
      ii = 1
      tmp = desfounders[names(desfounders) %in% id[fdi][currentnonfounders==iancestors[i]&active]]
      # Loop over final descendants 
      # check whether this is correct with more than one intermediate ancestors on the same level
      for (k in 1:length(tmp))
        {
        foundersdegreedes[[i]][(ii+1):(ii+length(tmp[[k]]))] = tmp[[k]]
        # Setting indicator of whether the descendant of all ancestors in tmp[[k]] is the previous intermediate ancestor
        iancestor.as.descendant[[i]][(ii+1):(ii+length(tmp[[k]]))] = list(ifelse (ia>1, names(tmp)[k] == iancestors[ia-1], FALSE))
        ii = ii + length(tmp[[k]])
        # Adding spouse of intermediate ancestor to list of founders of current final descendant
        
        }
      }
      iancestors.Nf[i] = ifelse(any(unlist(iancestor.as.descendant[[i]])),iancestors.Nf[ia-1],0) + length(iancestor.as.descendant[[i]])
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

# The dummy intermediate ancestor has all founders below him (except himself)
iancestors.Nf[ia] = Nf - 1 
    
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
    names(ancestorsdegreedes[[ia]]) = id[fdi][active]
    # Include first spouse in list of ancestors
    spousevec = unique(currentfounders)
    foundersdegreedes[[ia]]= list(degvec[active][currentfounders==spousevec[1]])
    # Setting indicator of whether the descendant is the previous intermediate ancestor
    if (ia>1)
      iancestor.as.descendant[[ia]] = list(id[fdi[active]][currentfounders==spousevec[1]] %in% lev.ia[[lia-1]])
    else iancestor.as.descendant[[ia]] = list(rep(FALSE,length(foundersdegreedes[[ia]][[1]])))
    # Add additional spouses if any
    # Warning! This is going to work only if all previous intermediate ancestors are under the same spouse
    if(length(spousevec)>1)
      {
      for (i in 2:length(spousevec))
        {
        foundersdegreedes[[ia]][[i]] = degvec[active][currentfounders==spousevec[i]]
        if (ia>1)
          iancestor.as.descendant[[ia]][[i]] = id[fdi[active]][currentfounders==spousevec[i]] %in% lev.ia[[lia-1]]
        else iancestor.as.descendant[[ia]][[i]] = rep(FALSE,length(foundersdegreedes[[ia]][[i]]))
        }
      }  
    # Include previous ancestors of final descendants if any
    if (any(names(desfounders) %in% id[fdi[active]]))
      {
      ii = length(spousevec)
      tmp = desfounders[names(desfounders) %in% id[fdi[active]]]
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
      j = 1
    for (i in (1:nfd)[active])
      {
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
# fdi in ouput contains indices of intermediate ancestors in addition to final descendants 
list(fdi=fdi,ia=ia,lev.ia=lev.ia,iancestors=iancestors,iancestor.as.descendant=iancestor.as.descendant,desfounders=desfounders,foundersdegreedes=foundersdegreedes,ancestorsdegreedes=ancestorsdegreedes,spousevec=spousevec)
}

pl = ped_datastruct(fdi,dv)
# Computation of numerator
num = 1
for (i in 1:pl$ia)
  num = num * 1/2^sum(pl$ancestorsdegreedes[[i]])
# Computation for top founder or founders
# If there is only one spouse, then a couple of founders can transmit a variant to all final descendents
if (length(pl$spousevec)==1) num = num*2
# Division by the number of founders
num = num/Nf

# Initialisation of counter of removed subjects
ncremoved = 0

if (missing(carriers))
numo = num
else
{
  carriers0 = carriers
  noncarriers = setdiff(id[!(id%in%dad.id | id%in%mom.id)],carriers)
  if (length(noncarriers)>0)
  {
  fd.subsets = list(as.matrix(carriers))
  # Loop over number of non-carriers to include as "carrier" in the possible subset
  if (length(noncarriers)>1)
    for (k in 1:(length(noncarriers)-1))
    {
  	  tmp = combn(noncarriers,k)
  	  fd.subsets[[k+1]] = rbind(matrix(carriers,length(carriers),ncol(tmp)),tmp)
    }
  subsetp = numeric(length(fd.subsets))
  # Loop over possible subsets
  for (k in 1:length(fd.subsets))	
    {
    sn = nrow(fd.subsets[[k]])
    # If there is only one final descendant or carrier in the subset, then the probability he is a carrier is 1/Nf
    if (sn == 1) subsetp[k] = 1/Nf
    else
    {
    nsubs = ncol(fd.subsets[[k]])
    subsetkp = numeric(nsubs)
  	for (h in 1:nsubs)
  	  {
  	  # Initialize meiosis reduction counter
      meir = 0
  	  # If there is more than one intermediate ancestor
  	  if (pl$ia>1)
  	    {
		insubset = logical(sn)
  	    # Loop over intermediate ancestors
  	    i = 1
  	    iia = numeric(0)
  	    done = FALSE
  	    for (lia in 1:length(pl$lev.ia))
  	    {
  	    # At first iteration, iia is length 0, so no iancestor gets added
  	    fdsi = c(fd.subsets[[k]][,h],pl$iancestors[iia])
  	    
		  for (cr in carriers)
		    {
		    # Check that the current carrier is a final descendant or intermediate ancestor in the current subset sharing a RV
		    if (!(cr %in% intersect(fd.subsets[[k]][,h],id[pl$fdi])))
		    {
		    # check if carrier is intermediate ancestor (other than the last)
		    if (cr %in% pl$iancestors[-pl$ia])
		      {
		      	# If he has at least one descendant in the current subset sharing a RV then discard him
		      if (any(fd.subsets[[k]][,h] %in% names(pl$ancestorsdegreedes[[which(pl$iancestors==cr)]])))
		        carriers = setdiff(carriers,cr)
		        # Else do nothing, and the intermediate ancestor should be recognized as a descendant of
		        # the intermediate ancestor above him, and treated as a final descendant
		      }
		    else
		      {
		      # check if carrier is parent of a final descendant or intermediate ancestor
		      if (cr %in% union(dad.id[pl$fdi],mom.id[pl$fdi]))
		      {
		      	# discard him
		        carriers = setdiff(carriers,cr)
		        # If he or she does not have at least one child among the current subset sharing a RV then add one of his children to the carriers in his place and increment the meiosis reduction counter
		      if (!any(fd.subsets[[k]][,h] %in% id[dad.id==cr | mom.id==cr]))
		        {
		        carriers = c(carriers,id[pl$fdi][dad.id[pl$fdi]==cr | mom.id[pl$fdi]==cr][1])    	
		        meir = meir + 1      	
		        }
		      }
		      else stop("Probability computations for subsets of carriers including ",cr," cannot be performed by RVsharing.")
		    }
		  }
		}
  	    insubset = c(insubset,rep(FALSE,length(iia)))
  	    for (ian in 1:length(pl$lev.ia[[lia]]))
  	    {
		    # If at least one final descendant in the current subset is below the current intermediate ancestor i, then
		    if (any(fdsi%in%names(pl$ancestorsdegreedes[[i]])))
		    {
		        # add that intermediate ancestor to the indices used in the computation for the current subset
		    	iia = c(iia,i)
		    	# records which final descendants from the subset are below that intermediate ancestor		    	
		    	insubset[fdsi%in%names(pl$ancestorsdegreedes[[i]])] = TRUE		    	
		    }
		    # Check if all carriers in current possible subset are descendants of current intermediate ancestor 
 	  		if (all(insubset))
 	  		{
 	  			# If this is the first intermediate ancestors at this level, or all final descendants are below the current intermediate ancestor, then the following intermediate ancestors
 	  			# are not needed for the computation
 	  		   if (ian == 1 | all(fdsi%in%names(pl$ancestorsdegreedes[[i]]))) {done = TRUE; break}
			}
			i = i+1
		 }
		 if (done) break
		 }
#		 print (done)
         # If not all final descendant in the current subset are below one of the intermediate ancestors up to level 1
         # then add the last intermediate ancestor at level 0         
		 if (!done) iia = c(iia,i)
		 }
		 # Else there is only one intermediate ancestor
		 else 
		 {
		 	done = FALSE
		 	iia = 1
		    for (cr in carriers)
		      {
		      # Check that the current carrier is a final descendant in the current subset sharing a RV
		      if (!(cr %in% intersect(fd.subsets[[k]][,h],id[fdi])))
		      {
		      # check if carrier is parent of a final descendant
		      if (cr %in% union(dad.id[fdi],mom.id[fdi]))
		      {
		      	# discard him
		        carriers = setdiff(carriers,cr)
		        # If he or she does not have at least one child among the current subset sharing a RV then add one of his children to the carriers in his place and increment the meiosis reduction counter
		      if (!any(fd.subsets[[k]][,h] %in% id[dad.id==cr | mom.id==cr]))
		        {
		        carriers = c(carriers,id[fdi][dad.id[fdi]==cr | mom.id[fdi]==cr][1])    	
		        meir = meir + 1      	
		        }
		      }
		      else stop("Probability computations for subsets of carriers including ",cr," cannot be performed by RVsharing.")
		      }
		    }	
		 }
		 
		 # Compute probability of subset
		 numsub = 1
		 # Correction for replacing parent by his child
		 if (meir>0) 
		 { 
		 	numsub = numsub * 2^meir
		 	# Substitution of carriers in the current subset
		 	fd.subsets[[k]][1:length(carriers),h] = carriers
		 }
		 for (ii in iia)
		 {
		 	# Here the intermediate ancestors in iancestors are in the list used in the computation for the current subset	
		      numsub = numsub * 1/2^sum(pl$ancestorsdegreedes[[ii]][unique(c(fd.subsets[[k]][,h],pl$iancestors[iia]))],na.rm=TRUE)
		 }
		 # If all members of the current subset have the same ancestor among the spouses of the final iancestor, then
		 # done is set to true to count transmissions from that ancestor in addition to those from the final iancestor
		 if (length(unique(sapply(pl$desfounders[intersect(fd.subsets[[k]][,h],as.character(id[fdi]))],function (vec) names(vec)[length(vec)-1])))==1) done = TRUE
		 # If there is only one spouse, then a couple of founders can transmit a variant to all final descendents
		 if (done | length(pl$spousevec)==1) numsub = 2*numsub
		 # Divide by the number of founders of the highest intermediate ancestor needed
		 subsetkp[h] = numsub/Nf
		 # set carriers vector back to initial carrier vector
		 carriers = carriers0
  	     }
  	  subsetp[k] = sum(subsetkp)
  	  }
  	}
  # Add joint prob for all final descendents
  subsetp = c(subsetp,num)
  # Computation of sharing probability of observed subset
  numo = sum(subsetp*(-1)^(0:(length(subsetp)-1)))
  }
  # Else there are no non-carriers, i.e. all affected subjects are carriers
  else numo = num

# Remove children of carriers (except intermediate ancestors) from list of final descendants if there are any
# for computation of denominator
  if (!all(carriers %in% id[pl$fdi]))
  {
  fdci = fdi
  for (cr in carriers)
    if (cr %in% union(dad.id[fdi],mom.id[fdi]) & !(cr %in% pl$iancestors))
      {
      	# Remove child of carriers from list of final descendants 
      	# (he can only have one, otherwise he would be intermediate ancestor)
      	fdci = setdiff(fdci,which(dad.id==cr | mom.id==cr))
    	# Add carrier to list of final descendants
    	fdci = c(fdci,which(id==cr))
    	# increment number of children removed
    	ncremoved = ncremoved + 1
      }
  pl = ped_datastruct(fdci,dv)
  }
}

# Computation of denominator
    	
# Probability that no variant has been transmitted
p0 = 0
# Probability that no variant has been transmitted from previous intermediate ancestor
pk = 1
for (i in 1:pl$ia)
  {
  for (j in 1:length(pl$foundersdegreedes[[i]]))
    p0 = p0 + prod((1-1/2^pl$foundersdegreedes[[i]][[j]]) + ifelse(pl$iancestor.as.descendant[[i]][[j]],(1/2^pl$foundersdegreedes[[i]][[j]])*pk,0))
  if (i < pl$ia)
    {
  # The previous intermediate ancestor becomes the current intermediate ancestor
  # Updates its probability, unless he is a carrier, in which case we keep it 1
  # For now, intermediate ancestors can have only one spouse, this is why we take the indicators of the first founder attached to him
    if (missing(carriers)) compute.pk = TRUE
    else 
      {
      if (pl$iancestors[i]%in%carriers) compute.pk = FALSE 
      else compute.pk = TRUE
      }
    if (compute.pk) pk = prod((1-1/2^pl$ancestorsdegreedes[[i]]) + ifelse(pl$iancestor.as.descendant[[i]][[1]],1/2^pl$ancestorsdegreedes[[i]]*pk,0))
    else pk=0
    }
  }
# At the end, add the probability from the dummy "intermediate" ancestor. He is currently the only one who can have more than one spouse
# Since only one of his spouses can be the parent of the previous intermediate ancestor, ifelse returns only one non-zero term.
  # Debugging code
  # print (foundersdegreedes[[i]])
  # print (ancestorsdegreedes[[i]])
  # print (iancestor.as.descendant[[i]])
  # print (p0)
  tmpf = ifelse(unlist(pl$iancestor.as.descendant[[i]][1:length(pl$spousevec)]), (1/2^pl$ancestorsdegreedes[[i]]) * pk,0)
  p0 = p0 + prod((1-1/2^pl$ancestorsdegreedes[[i]]) + tmpf)

#  tmpf = as.list(sapply(pl$iancestor.as.descendant[[i]][1:length(pl$spousevec)],function(lv,deg,pk) ifelse(lv, (1/2^deg) * pk,0), deg=pl$ancestorsdegreedes[[i]],pk=pk))  
#  p0 = p0 + prod((1-1/2^pl$ancestorsdegreedes[[i]]) + sapply(tmpf,sum))
  # If children have been removed, add the probability that their other parent did not transmit them the RV
  p0 = p0 +  (1/2)*ncremoved
    # Debugging code
    # print (p0)
    
  # sharing probability
  pshare = numo/(1-p0/Nf)   
  
if (missing(carriers)) carrier.vec = as.character(id[fdi])
else carrier.vec=as.character(carriers0)
new("RVsharingProb",pshare=pshare,iancestors=pl$iancestors,desfounders=pl$desfounders,id=as.character(id),dad.id=as.character(dad.id),mom.id=as.character(mom.id),carriers=carrier.vec)
}

# Wrappers for pedigree object
# Returns only pshare
RVsharing.ped.pshare = function(ped)
{
RVsharing(ped)@pshare
} 

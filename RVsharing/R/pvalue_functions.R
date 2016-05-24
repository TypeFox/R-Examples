# Function to compute p-value when variants are shared in a subset of families
get.psubset = function(vec,not,pshare.data)
{
	# vec : vector of names of all families where a variant is seen
	# not : vector of names of families where not all affected subjects share the rare variant (RV)
	# pshare.data : data frame with at least the two following columns:
	#    pshare : vector of RV sharing probabilities
	#    ped.tocompute.vec : vector of names of the families whose sharing probability is contained in pshare. 
	#						 The names in the arguments "vec" and "not" must be found in ped.tocompute.vec.

	# check: "not" contains at least one family 
	if (length(not)==0) stop ("Vector 'not' of families not sharing the RV is empty.")
    # check: all families in not are also in vec
    if (!all(not%in%vec)) stop("Vector 'not' contains families not in 'vec'")
    
    # check: all families in vec are also in pshare.data$ped.tocompute.vec
    if (!all(vec%in%pshare.data$ped.tocompute.vec)) stop("Vector 'vec' contains families not in 'pshare.data$ped.tocompute.vec'")
        	
	# If all families share the variant, then return 1
	if (length(not)==length(vec)) return (1)
		
p.vec = pshare.data$pshare[pshare.data$ped.tocompute.vec%in%vec]
names(p.vec) = pshare.data$ped.tocompute.vec[pshare.data$ped.tocompute.vec%in%vec]
nf = length(p.vec)
nnot = sum(names(p.vec)%in%not)

# Probability of observed data
p.obs = prod(p.vec[!(names(p.vec)%in%not)],1-p.vec[as.character(not)])
#print (p.obs)

# Tail probability includes case where all families share the variant
p = prod(p.vec)
for (h in 1:nnot)
{
	comb.mat = combn(nf,h)
#	print(comb.mat)
# plus the cases where the probability with two families not sharing the variant is less extreme than the observed
# Compute probability for all pairs of families not sharing
for (i in 1:ncol(comb.mat))
  {
    ptmp = prod(p.vec[-comb.mat[,i]],1-p.vec[comb.mat[,i]])
#    print(ptmp)
    if (ptmp <= p.obs) p = p + ptmp
  }
}
p
}

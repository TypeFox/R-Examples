
flsa.general <-
function(A,L,lam1,lam2,penalize.diagonal)
{
trueA = A
if(is.matrix(A[[1]])) {p=dim(A[[1]])[1]}
if(is.vector(A[[1]])) {p=length(A[[1]])}
K = length(A)
# results matrices:
X = list()
#for(k in 1:K) {X[[k]] = matrix(NA,p,p)} 
for(k in 1:K) {X[[k]] = A[[1]]*NA} 
if(is.matrix(A[[1]])) {fusions = array(FALSE,dim=c(K,K,p,p))}  
if(is.vector(A[[1]])) {fusions = array(FALSE,dim=c(K,K,p,1))}

# get starting newc: list of matrices.  newc[[k]][i,j] gives the (2*ordermats[[k]]-K-1) adjustment for lam2 at that element.  Basically, how many rank higher vs how many rank lower?
newc = list()
for(k in 1:K)
{
	others = setdiff(1:K,k)
	others.smaller.k = 1:(k-1)
	newc[[k]] = A[[1]]*0               
	for(o in others) {newc[[k]] = newc[[k]] + (A[[o]]-A[[k]]< -1e-4) - (A[[o]]-A[[k]]>1e-4)} 
}

######### start the loop here:
for(iter in 1:(K-1))
{

# create order matrices:
ordermats = list()
for(k in 1:K)
{
	others = setdiff(1:K,k)
	others.smaller.k = 1:(k-1)
	ordermats[[k]] = A[[1]]*0   
	for(o in others) {ordermats[[k]] = ordermats[[k]] + (A[[k]]-A[[o]]>1e-4)} 
	# to deal with ties, also add a unit to ordermat[[k]] if a class with a lower k has a tie at element i,j:
	if(k>1)
	{
		for(o in others.smaller.k) {ordermats[[k]] = ordermats[[k]] + (abs(A[[o]]-A[[k]])<1e-4)} 
	}
	ordermats[[k]] = ordermats[[k]] + 1
}

# create beta.g matrices, holding the solution to Holger's "unconstrained problem" 
#  (prending we're not constraining the order of the solution to match the order of the A matrices)
betas.g = list()
for(k in 1:K)
{
	betas.g[[k]] = A[[k]] - lam2/L*newc[[k]]
}

# identify and fuse all elements for which the betas.g are out of order:
new.ordermats = list()
for(k in 1:K)
{
	others = setdiff(1:K,k)
	others.smaller.k = 1:(k-1)
	new.ordermats[[k]] = A[[1]]*0   
	for(o in others) {new.ordermats[[k]] = new.ordermats[[k]] + (betas.g[[k]]-betas.g[[o]]>1e-4)}  
	# to deal with ties, also add a unit to ordermat[[k]] if a class with a lower k has a tie at element i,j:
	if(k>1)
	{
		for(o in others.smaller.k) {new.ordermats[[k]] = new.ordermats[[k]] + (abs(betas.g[[o]]-betas.g[[k]])<1e-4)} 
	}
	new.ordermats[[k]] = new.ordermats[[k]] + 1
}

# identify neighboring fusions:  "fusions": K x K x p x p array: K x K matrices, T/F for fusions 
for(k in 1:K){
for(kp in 1:K){
	#given k,kp, declare a fusion when their ordermats entries are adjacent, and their new.ordermats entries have the opposite direction:
	fusions[k,kp,,] = fusions[k,kp,,]+           
				((ordermats[[k]]-1==ordermats[[kp]])&(new.ordermats[[k]]<new.ordermats[[kp]]))+
				((ordermats[[k]]+1==ordermats[[kp]])&(new.ordermats[[k]]>new.ordermats[[kp]]))+
				(abs(A[[k]]-A[[kp]])<1e-4)
				#(existing fusions, neighboring fusions, and ties)
	fusions = (fusions>0)*1
}}
                            

# now we've noted fusions between all entries with adjacent ordermats entries and reversed new.ordermats entries
# next: extend fusions to non-adjecent entries: if a-b and b-c, then connect a-c:
for(k in 1:K){
for(kp in 1:K){
	others = setdiff(1:K,c(k,kp))
	for(o in others)
	{
		#identify elements in o which are fused with the same element in both k and kp, then add them to the list of k-kp fusions:
		bothfused = fusions[k,o,,] & fusions[kp,o,,]    
		fusions[k,kp,,] = fusions[k,kp,,] | bothfused       
	}
}}

# now recalculate A with the new fused entries:
# to recalculate A, for each non-zero entry, identify the classes k with which it must be fused, and get their average:
for(k in 1:K)
{
	others = setdiff(1:K,k)
	#fusemean and denom: the mean value of all the trueA to be fused, and the number of values to be fused:
	fusemean = trueA[[k]]
	denom = A[[1]]*0+1   
	for(o in others)
	{
		fusemean = fusemean+fusions[k,o,,]*trueA[[o]]  #add the values of the elements which must be fused to fusemean   
		denom = denom+fusions[k,o,,]     
	}	
	# now redefine A[[k]]: unchanged from trueA if there's no fusion, and the mean of the fused elements when there is fusion:
	A[[k]] = fusemean/denom
}

#newc: list of matrices.  newc[[k]][i,j] gives the (2*ordermats[[k]]-K-1) adjustment for lam2 at that element.  Basically, how many rank higher vs how many rank lower?
newc = list()
for(k in 1:K)
{
	others = setdiff(1:K,k)
	others.smaller.k = 1:(k-1)
	newc[[k]] = A[[1]]*0   
	for(o in others) {newc[[k]] = newc[[k]] + (A[[o]]-A[[k]]< -1e-4) - (A[[o]]-A[[k]]>1e-4)}  
}

} #end loop here

# final version of betas.g:
for(k in 1:K)
{
	betas.g[[k]] = A[[k]] - lam2/L*newc[[k]]
}
# now soft-threshold the solution matices:
for(k in 1:K)
{
	X[[k]] = soft(betas.g[[k]],lam=lam1/L,penalize.diagonal=penalize.diagonal)
}
return(X)
}




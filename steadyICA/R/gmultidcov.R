#------------------
# Nicholas James
# July 18, 2013


# Give a membership vector group this method returns the group version of the symmetric 
# or asymmetric multivariate distance covariance statistic.

gmultidcov = function(S,group=1:ncol(S),alpha=1,symmetric=TRUE){
	# Perform checks on arguments
	if(!is.matrix(S))
		stop("The argument 'S' must be a matrix.")
	if(!is.vector(group))
		stop("The argument 'group' a vector.")
	if(!(symmetric==TRUE || symmetric==FALSE))
		stop("The argument 'symmetric' must be able to be evaluated as a boolean.")
	if(length(group) != ncol(S))
		stop("The length of 'group' and the number of columns in 'S' must be equal.")
	# Rename group labels so as to be in {1,2,...,C}
	u = unique(group)
	C = length(u)
	ans = 0
	for(i in 1:C)
		group[group==u[i]] = i
	# Calculate symmetric or asymmetric statistic
	if(symmetric)
		for(i in 1:C)
			ans = ans + dcovustat(S[,group==i],S[,group!=i],alpha)
	else
		for(i in 1:(C-1))
			ans = ans + dcovustat(S[,group==i],S[,group>i],alpha)
	return(ans)
}


#------------------
# Nicholas James
# July 15, 2013

# Calculate a complete empirical measure of mutual multivariate independence.
# Makes use of the utils::combn and combinat::permn functions, which are written in R. 
# So a speed up may be possible if it were to be rewritten in C/C++. This also has the 
# potential to be very memory intensive.

compInd = function(S,group=1:ncol(S),alpha=1){
	#Perform checks on arguments
	if(!is.matrix(S))
		stop("The argument 'S' must be a matrix.")
	if(nrow(S) < 2)
		stop("The matrix 'S' must have more than 2 rows.")
	if(!is.vector(group))
		stop("The argument 'group' must a vector.")
	if(length(group) != ncol(S))
		stop("The length of 'group' and the number of columns in 'S' must be equal.")
	if(alpha <= 0 || alpha>2)
		stop("The argument 'alpha' must be in the interval (0,2].")
	n = unique(group)
	m = length(n)
	for(i in 1:m)
		group[group==n[i]] = i
	n = nrow(S)
	if(n < 2*m)
		stop("The number of observations must be at least twice the number of groups.")
	d = ncol(S)
	x1 = x2 = numeric(d)

	J11 = J12 = J22 = 0
	# Calculate within sample distances
	J11 = sum(dist(S)^alpha)
	J11 = 2*J11/(n*(n-1))
	
	U1 = utils::combn(n,2*m)
	for(i in 1:ncol(U1)){
		I = sort(U1[,i]) # sort the times
		J = I[(m+1):(2*m)] # the last m go to the j indices
		I = I[1:m] # the first m go to the i indices
		I2 = matrix(unlist(combinat::permn(I)),byrow=F,nrow=m)
		J2 = matrix(unlist(combinat::permn(J)),byrow=F,nrow=m)
		for(j in 1:ncol(I2)){
			for(k in 1:ncol(J2)){
				for(x in 1:d){
					x1[x] = S[I2[group[x],j],x]
					x2[x] = S[J2[group[x],k],x]
				}
				J22 = J22 + as.numeric(dist(rbind(x1,x2))^alpha)
			}
		}
	}
	J22 = J22 / ( choose(n,2*m)*factorial(m)*factorial(m) ) 

	# Calculate between distances
	for(i in 1:n){
		U2 = utils::combn((1:n)[-i],m)
		x1 = S[i,]
		for(j in 1:ncol(U2)){
			I = matrix(unlist(combinat::permn(U2[,j])),byrow=F,nrow=m)
			for(k in 1:ncol(I)){
				for(x in 1:d)
					x2[x] = S[I[group[x],k],x]
				J12 = J12 + as.numeric(dist(rbind(x1,x2))^alpha)
			}
		}
	}
	J12 = J12/(n * factorial(m) * choose(n-1,m))
	return(2*J12 - J11 - J22)
}

#apxCompInd = function(S,group=1:ncol(S),alpha=1){
#
#}



#------------------
# Nicholas James
# July 31, 2013

# Permutation test used to test for independence between components (or grouped components).
# Test is based on the function gmultidcov.

permTest = function(S, group=1:ncol(S), R=199,FUN = c('gmultidcov','compInd'),...){
	# Perform checks on arguments
	if(!is.matrix(S))
		stop("The argument 'S' must be a matrix.")
	if(!is.vector(group))
		stop("The argument 'group' a vector.")
#	if(!(symmetric==TRUE || symmetric==FALSE))
#		stop("The argument 'symmetric' must be able to be evaluated as a boolean.")
	if(length(group) != ncol(S))
		stop("The length of 'group' and the number of columns in 'S' must be equal.")
	R = floor(R)
	if(R < 0)
		stop("The argument 'R' must be nonnegative.")
	FUN = match.arg(FUN,c('gmultidcov','compInd'))
	FUN = match.fun(FUN)

	# Rename group labels so as to be in {1,2,...,C}
	u = unique(group)
	C = length(u)
	for(i in 1:C)
		group[group==u[i]] = i
	Sgroup = sort(group)

	# Obtain observed statistic
	obs = FUN(S,group,...)

	over = 0
	n = nrow(S)
	perm = matrix(1:n,nrow=n,ncol=C,byrow=FALSE)
	for(r in 1:R){
		S1 = NULL
		perm = apply(perm,2,sample) # Permute time for each group
		for(i in 1:C) # Make permuted observation matrix
			S1 = cbind(S1,S[perm[,i],group==i])
		stat = FUN(S1,Sgroup,...)
		if( stat >= obs)
			over = over+1
	}
	return( (over+1)/(R+1) )
}

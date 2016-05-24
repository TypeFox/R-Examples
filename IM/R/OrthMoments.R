################# Compute Discrete/Continuous Orthogonal moments  ##############################
# 
# Author: Tan Dang, Allison Irvine
###############################################################################



####### can choose one type or a combination of 2 ###################################
# KR: Krawchouk normalized (discrete)
# TC: Tchebichef normalized (discrete)
# H: Hahn normalized (discrete)

# Param:list of length 2, one numeric array for x and one numeric array for y
#first number in params[[i]] is order

# 1 constant between (0,1) for KR
# integer a,b, for Hahn

#case:list of length 2 - type of polynomials for x, for y
#DISCRETE:
#"cheby"		Chebyshev
#"krawt"		Krawtchouk
#"hahn"			Hahn
#CONTINUOUS:
#"gegen"		Gegenbauer
#"chebycont"	ChebyshevCont	(continuous chebyshev)
#"legend"		Legendre
##########################################


# Compute Moment
OrthMoments <- function(I,polynomials){

	# Computing moments
	M=t(polynomials[[2]])%*%I%*%(polynomials[[1]])
	return(M)
}

# Reconstruct image using moments up to a specified order for x and y
OrthReconstruct<- function(M,polynomials,d,order){
	order = order+1
	# get xMat, yMat - polynomial values at x-y direction in various order
	polyX= polynomials[[1]][,1:order[1]];
	polyY= polynomials[[2]][,1:order[2]];

	# create blank image
	IR= mat.or.vec(d[1],d[2])
	#M= zeroMoment(M[1:order[2],1:order[1]]);
	
	# loop through all possible order combination
	IR=polyY%*%(M[1:order[2],1:order[1]])%*%t(polyX)
	
	#IR= (IR-min(IR))/(max(IR)-min(IR))
	
	IR
}




# precompute polynomial values for x and y 
#pass centroid if using a continuous polynomial
OrthPolynomials <- function(case,d,order,param,Invariant=0,center=c(0,0)){
	
	
	# get x matrix
	polyX=switch(case[1], 
				krawt=Krawtchouk(order[1],d[2],param[[1]],center[1],Invariant),
				cheby=Chebyshev(order[1],d[2],center[1],Invariant),
				hahn=DualHahn(order[1],d[2],param[[1]],center[1],Invariant),
				gegen = Gegenbauer(order[1],d[2],param[[1]],center[1],Invariant),
				chebycont = ChebyshevCont(order[1],d[2],center[1],Invariant),
				legend = Legendre(order[1],d[2],center[1],Invariant)
				)
	
	# get y matrix
	polyY=switch(case[2], 
				krawt=Krawtchouk(order[2],d[1],param[[2]],center[2],Invariant),
				cheby=Chebyshev(order[2],d[1],center[2],Invariant),
				hahn=DualHahn(order[2],d[1],param[[2]],center[2],Invariant),
				gegen = Gegenbauer(order[2],d[1],param[[2]],center[2],Invariant),
				chebycont = ChebyshevCont(order[2],d[1],center[2],Invariant),
				legend = Legendre(order[2],d[1],center[2],Invariant)
				)
	
	if (Invariant){
		polyX=polyX/(sqrt(Invariant));
		polyY=polyY/(sqrt(Invariant));
	}
	list(polyX,polyY)
}


# Weighted Krawtchouk
Krawtchouk <- function(order,N,param,center=0,Invariant=0){
	center= round(center)
	p=param[1]
	
	#check parameter constraints
	if ((p<0) || (p>1)){
		return(NULL)
	}
		
	
	n = 0:(order);
	x = 0:(N-1);
	# Scaling invariant
	if (Invariant) x=x/sqrt(Invariant);
	N2=N-1;
	
	w = KrawWeights(N,p,x);
	
	#if (Invariant){
		# Translate to centroid
	#	w = c(w[(N-center+1):N], w[1:(N-center)]); 
	#	x = c(x[(N-center+1):N], x[1:(N-center)]); 
	#}
	#rows are n, columns are x
	#depends on n

	K = mat.or.vec(N,order+1);
	K[,1] = sqrt(w);
	K[,2] = (1-(x/(p*N2))) * sqrt((p*N2)/(1-p)) * sqrt(w);
	
	for (i in 2:(order-1)) {
		#i is equal to n
		#using n1 to be equal to n-1
		n1 = i-1;
		
		c1 = sqrt( (p*(N2-n1))/((1-p)*(i)) ) * (N2*p - 2*n1*p + n1 - x);
		c2 = sqrt( ((p^2)*(N2-n1)*(N-n1))/(((1-p)^2)*i*n1) ) * n1 * (1-p);
		c1 = c1/(p*(N2-n1));
		c2 = c2/(p*(N2-n1));
		
		K[,i+1] = c1*K[,i] - c2*K[,i-1];
		
	}

	K
}

# Weight Function of Krawtchouk
KrawWeights <- function(N,p,x) {
	N2 = N-1;

	w = ((N2-x+1)/x)*(p/(1-p));
	w[1] = (1-p)^N2;
	for (i in 1:N2) {
		w[i+1] = w[i+1] * w[i];
	}
	
	w
	
}

# Weighted Tchebichef
Chebyshev <- function(order,N,center=0,Invariant=0){
	center= round(center);
	x= 0:(N-1);

	TC=as.matrix(mat.or.vec(order+1,N))
	
	TC[1,1]= 1/sqrt(N);
	for (n in 1:(order)){
		TC[n+1,1]= -TC[n,1]*sqrt((N-n)/(N+n))*sqrt((2*n+1)/(2*n-1));
	}
	
	n= 0:(order);
	# Scaling invariant
	if (Invariant) TC[,1]= TC[,1]/((sqrt(Invariant))^(0:order));
	
	TC[,2]= TC[,1]*(1+n*(n+1)/(1-N));
	
	for (X in 2:(N-1)){
		if (X<= ceiling(N/2)){
			TC[,X+1] = 1/(X*(N-X))*(((-n)*(n+1)-(2*X-1)*(X-N-1)-X)*TC[,X]+(X-1)*(X-N-1)*TC[,X-1]);
		} else {
			TC[,X+1]= (-1)^n*TC[,N-X];
		}
	}
	
	#if (Invariant){
	#	# Translate to centroid
	#	TC = cbind(TC[,(N-center+1):N], TC[,1:(N-center)]);
	#}
	TC=t(TC)

	TC
}

# Scale moment of Tchebichef discrete
scaleTC<- function(TC,N,n){
	TC= TC*sqrt((2*n+1)/N)
	
	if (n>0){
		for (i in 1:n){
			TC= TC/sqrt(N^2-i^2)
		}
	}
	
	TC
}


# Dual-Hahn polynomial
#param: a,c are integers, a > -1/2, a > (abs(c)-1)
DualHahn<- function(order,N,param,center=0,Invariant=0){
	center= round(center)
	# Get parameter
	a=param[1]
	b=a+N
	c=param[2]
	
	# Check nmax<= N-1
	if (order > N-1) order = N-1
	
	# Validate constrain
	if (a<= (-1/2)) {
		return(NULL)	
	}
	if (a<= (abs(c)-1)){
		return(NULL)
	}
	
	# Compute s
	S=seq(a,b-1,1)
	
	# Scaling invariant
	if (Invariant) S=S/sqrt(Invariant); 
	
	# Create H matrix to store poly value, x is column, order is row
	H=as.matrix(mat.or.vec(N,order+1))
	
	# Compute rho(s)/d^2 
	rhoD= array(0,c(N,2))
	
	# Compute rho(a)/d(0)^2
	rhoD[1,1]= 1/(2*a+1)
	for (i in 1:(N-1)){
		rhoD[1,1]=rhoD[1,1]*((a+i-c)/(2*a+i+1))	
	}
	# Compute rho(s)/d(0)^2
	for (Is in 1:(N-1)){
		s=S[Is]
		rhoD[Is+1,1]=rhoD[Is,1]*(a+s+1)*(c+s+1)*(b-s-1)/(s-a+1)/(b+s+1)/(s-c+1)	
	}
	
	#if (Invariant){
		# Translate to centroid
	#	rhoD[,1]= c(rhoD[(N-center+1):N,1], rhoD[1:(N-center),1]);
	#	S = c(S[(N-center+1):N], S[1:(N-center)]);
	#}

	# Compute rho(s)/d(1)^2
	rhoD[,2]=rhoD[,1]/(a+c+1)/(N-1)/(b-c-1)
	
	
	# Compute H
	H[,1]=sqrt(rhoD[,1]*(2*S+1))
	H[,2]=(-(a+S+1)*(c+S+1)*(b-S-1)+(S-a)*(S+b)*(S-c))*sqrt(rhoD[,2]/(2*S+1))
	
	for (n in 1:(order-1)){
		A= 1/(n+1)*(S*(S+1)-a*b+a*c-b*c-(N-c-1)*(2*n+1)+2*n^2)
		B= 1/(n+1)*(a+c+n)*(N-n)*(b-c-n)
		T1=(n+1)/(a+c+n+1)/(N-n-1)/(b-c-n-1)
		T2=T1*n/(a+c+n)/(N-n)/(b-c-n)
		H[,n+2]=A*sqrt(T1)*H[,n+1]-B*sqrt(T2)*H[,n]
	}
	
	H
}


# Gegenbauer Polynomial
Gegenbauer <- function(order,N,param,centerX,Invariant){
	#map coordinates around centroid of image
	x = mapCenter(N,centerX,Invariant);

	# get lambda
	lambda=param
	
	# create matrix 
	G= as.matrix(mat.or.vec(N,order+1))
	
	# initial value
	G[,1]= 1
	G[,2]= 2*lambda*x
	
	# recursion
	if (order>=2){
		for (n in 1:(order-1)){
			G[,n+2]=(2*n+lambda)/(n+lambda)*x*G[,n+1]-(n+2*lambda-1)/(n+1)*G[,n]
		}
	}
	
	for (n in 0:order){
		if (ceiling(2*lambda)==(2*lambda)){
			#cnst= (n+1)*...*(n+2*lambda-1)
			cnst=n+1;
			if (2*lambda>=2){
				for (i in 3:ceiling(2*lambda)){
					cnst= cnst*(n+i-1);
				}
			}
		} else {
			cnst= gamma(n+2*lambda)/factorial(n);
		}
		if (lambda==1) cnst= (n+1);
		if (n==0) cnst=gamma(2*lambda);
		G[,n+1]=G[,n+1]*sqrt((1-x^2)^(lambda-1/2)*(gamma(lambda))^2*(n+lambda)
				/cnst/pi/2^(1-2*lambda))
	}
	G
}

# Chebyshev (continuous) polynomial
ChebyshevCont <- function(order,N,centerX,Invariant){
	#map coordinates artound centroid of image
	x = mapCenter(N,centerX,Invariant)

	# create matrix 
	CH= mat.or.vec(N,order+1)
	
	# initial values
	CH[,1]= 1
	CH[,2]= x
	
	# recursion
	if (order>=2){
		for (n in 1:(order-1)){
			CH[,n+2]=2*x*CH[,n+1]-CH[,n]
			#CH[,n+1]= CH[,n-1]-2*x*CH[,n]
		}	
	}
	
	for (n in 1:order){
		CH[,n+1]=CH[,n+1]*sqrt(2/pi)*(1-x^2)^(-1/4);
	}
	
	CH[,1]= CH[,1]*sqrt(1/pi)*(1-x^2)^(-1/4);
	
	CH
}


# Legendre polynomial
Legendre <- function(order,N,centerX,Invariant){
	#map coordinates artound centroid of image
	x = mapCenter(N,centerX,Invariant)

	L= as.matrix(mat.or.vec(N,order+1))
	
	L[,1]=1
	L[,2]=x
	
	for (n in 1:(order-1)){
		L[,n+2]=1/(n+1)*((2*n+1)*x*L[,n+1]-n*L[,n])
	}
	
	for (n in 0:order){
		L[,n+1]=L[,n+1]*sqrt((2*n+1)/2)
	}
	
	L
}

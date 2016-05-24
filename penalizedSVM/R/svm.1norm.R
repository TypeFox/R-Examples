lpsvm<-function (A, d, k=5, nu=0,output=1, delta=10^-3, epsi=10^-4,  seed=123, maxIter=700 ){
#function [w,gamma,trainCorr, testCorr,  nu]=lpsvm(A,d,k,nu,output,delta)
## this implements L1 SVM classification, actually a fast Newton algorithm NLPSVM from Fung and Mangasarian(2004)
# lpsvm = nlp svm = L1 svm = 1-norm svm
# 
#  
# version 1.3
# last revision: 07/07/03
#===========================================================================
# Usage: lpsvm(A,d,k,nu,output,delta);
#
# A and d are both required, everything else has a default
# An example: [w gamma train test time nu] = lpsvm(A,d,10);
#
#  	# separating hyperplane: A *w = gamma !
#
# Input parameters:
#    A: Data points = inputmatrix in R (m,n ); A row= pat, col=feature
#    d: 1's or -1's length= m 
#    k: k-fold for cv 
#     way to divide the data set into test and training set
#       if k = 0: simply run the algorithm without any correctness
#         calculation, this is the default
#       if k = 1: run the algorithm and calculate correctness on
#         the whole data set
#       if k = any value less than the number of rows in the data set:
#         divide up the data set into test and training
#         using k-fold method
#       if k = number of rows in the data set: use the 'leave 1' method
#		epsi: tuning parameter	
#
#    output= verbose: 0 - no output, 1 - produce output, default is 0
#    nu:             weighted parameter
#                    1 - easy estimation
#                    0  - hard estimation 
#                    any other value - used as nu by the algorithm
#                    default - 0
#    delta:  default is 10^-3
#===================================================================
# Output parameters:
#    for separating hyperplane: Ax + b = 0 ! 
#
#       w:              the normal vector of the classifier
#       b= (-1)*  gamma:          the threshold
#       trainCorr:      training set correctness
#       testCorr:       test set correctness
#       nu:             estimated value (or specified value) of nu
#==========================================================================
	
	require(statmod) # for matrixmultplications 
	require(corpcor)
	
	# be shure that d is a vector of 1 and -1
	tmp<-matrix(as.numeric(as.character(d)))
	rownames(tmp)<-names(d)
	d<-tmp
	#str(d)
	# d is a 1-row matrix: 
	
	if (nu==0){nu = .EstNuLong(A,d)}       # default is hard estimation
	if (nu==1){nu = .EstNuShort(A,d)} 
	# else the use defined nu
	
	set.seed(seed)
	r<- sample(1:length(d))
	d<-d[r]
	A <- A[r, ]
	
	trainCorr=0
	testCorr=0
	
	# flag: non empty model for cv
	non.empty<-rep(NA, k)
	
	
	# in any case calculate the classifier 
	List<-  .core(A=A,d=d,nu=nu,delta=delta, epsi=epsi,maxiter=maxIter)#
	if (!is.null(List)){
			w<-  List$w
			gamma <- List$gamma
			iter <- List$iter
			xind<-  List$xind
			
			if (k==1){#if k==1 only training set correctness is calculated
				trainCorr = round(.correctness(AA=A[,xind, drop=FALSE],dd=d,w=w,gamma=gamma), 1);
				if (output==1) print(paste("Training set correctness: ",trainCorr, "%"));
			}  
			if (output==1) 	print(paste("Number of Iterations:" ,iter))
	}else {
			if (output==1) print("empty model")
	} # 
		 
	
	if  (k >1){
	# do k-fold CV	
			
		sm<-nrow(A)
		sn<- ncol(A)
		# accumulative iterations
		accuIter = 0;
		
		num.w<- 0 # average number of genes in the final model
		
		indx = c(0:k);
		indx = floor(sm*indx/k);    #last row numbers for all 'segments'
		# split trainining set from test set
		for  (i in  1:k){
						
			# ist ok to take the first nrow/k rows , because we have already permuted them.
			Ctest = A[(indx[i]+1):indx[i+1], , drop=FALSE];
			dtest = d[(indx[i]+1):indx[i+1]];
			
			# everything ecxept Ctrain 
			Ctrain = A[- ((indx[i]+1):indx[i+1]), , drop= FALSE	];
			dtrain = d[- ((indx[i]+1):indx[i+1])];
			
			# in k fold cv don's store coefs for each iteration !!!!
			List.cv = .core(A=Ctrain,d=dtrain,nu,delta, epsi,maxiter=maxIter);
			
			# if the list is not empty
			if (!is.null(List.cv)){
				non.empty[i]<- TRUE
				w.cv<-  List.cv$w
				gamma.cv <- List.cv$gamma
				iter.cv <- List$.cviter
				xind.cv<-  List.cv$xind
				
				tmp.num.w<-length(w.cv)
				tmpTrainCorr = round(.correctness(Ctrain[,xind.cv, drop=FALSE],dtrain,w=w.cv,gamma=gamma.cv),1)
				tmpTestCorr = round(.correctness(Ctest[,xind.cv, drop=FALSE],dtest,w=w.cv,gamma=gamma.cv),1);
				
				if (output == 1){
					print("________________________________________________");
					print(paste("Fold",i));
					print(paste("Training set correctness:",tmpTrainCorr));
					print(paste("Testing set correctness:",tmpTestCorr));
					print(paste("Number of iterations:",iter.cv));
					print("selected genes:")
					print(w.cv)
				
				}
				
				trainCorr = trainCorr + tmpTrainCorr
				testCorr = testCorr + tmpTestCorr
				accuIter = accuIter + iter.cv # accumulative iterations
				num.w<- num.w + tmp.num.w
			}else{
				non.empty[i]<- FALSE
				if (output == 1){
					print("________________________________________________");
					print(paste("Fold",i));
					print("empty model");
				}
			}
		} # end of for (looping through test sets)
			
			# take in accout the number of non-empty Lists (models)!
			trainCorr = trainCorr/sum(non.empty)
			testCorr = testCorr/sum(non.empty)
			
			if (output == 1){
				print("==============================================");
				print(paste("Training set correctness:",trainCorr))
				print(paste("Testing set correctness:" ,testCorr))
				print(paste("Average number of iterations: ",accuIter/sum(non.empty)))
				print(paste("Average number of genes: ",num.w/sum(non.empty)))
			}
	}
		
	ret<- list(w=w, b=(-1)*gamma, xind=xind, lambda1=epsi, iter =iter, trainCorr=trainCorr,testCorr=testCorr, nu=nu, maxIter=maxIter )
	class(ret) <- "1norm"
	return(ret)
}


`.core` <-
function(A,d,nu,delta, epsi, tol=10^(-3), alpha=10^3, maxiter=50){
	
	num.samp<-dim(A)[1]
	num.clones<-dim(A)[2]
	
	m<-dim(A)[1]
	n<-dim(A)[2]
	en<- matrix(rep(1, n))
	em<- matrix(rep(1, m))
	
	#initial u
	u=matrix(rep(1, m));
	iter=0;
	epsi=epsi*em;
	nu=nu*em;
	diff=1;
	#DA=spdiags(d,0,m,m)*A;  * matrix multiplication 
	# A = spdiags(B,k,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by k.
	# d a matrix, set the columns' values into the main diagonal (k=0), (k <0 below  and k>0 above the main diaganal)
	#  spdiags(d,0,m,m) <=> create m*m matrix with elements of d in the main diagonal -> diag(as.vector(d)) !!!!!
	DA<- vecmat(as.vector(d), A )  
	
	while ( (diff>tol) & (iter<maxiter) ){
		# Stop if ||u -  u_old || <= tol or i = imax. Else, set i = i +1, a = 2a and  run i+1 iteration 
		uold <-u;
		iter <- iter+1;
		if (iter >1) alpha<- 2*alpha
		
		du=d*u;
		Adu=t(A) %*% du;
		# replase negative values by 0
		pp=Adu-en; pp[pp<0]<-0 # positive
		np=-Adu-en; np[np<0]<- 0 # negative 
		# sum of matrix = apply(matrix, 2, sum) !!!! 
		#dd=sum(du)*d;unu=max(u-nu,0);uu=max(-u,0);
		dd=sum(du)*d;
		unu=u-nu; unu[unu <0]<- 0;
		uu=-u; uu[uu <0]<- 0 
		#Gradient
		g= -epsi+(d*(A %*% pp))-(d*(A%*%np))+dd+unu-alpha*uu;
		#Hessian
		E_vec<- sqrt(sign(np)+sign(pp))
		H=cbind(matvec(DA,E_vec), d);
		
		f=as.vector(delta+sign(unu)+alpha*sign(uu) )
		
		#if the Q matrix is not large:  num pat < num clones = n.thr invert it directly
		# else apply  SMW FORMula
				# see Feature SelectionSforSVMC_05204tt_1normSVM.ppt, slide 20 
				# use SMW FORMula : invert a quadratic matrix Q: = U * D * U' + A  
		
		if (num.samp>= num.clones){ 
		#with_smw(A,d,nu,delta, epsi=epsi)
			inv_Q<- .find.inverse(U=H, D_vec=rep(1,n+1), A_vec=f, n.thr=1000)
		}else {
			# calulate Q_1 directly
			F<- diag(f)
			Q = H %*% t(H)+F
			# TODO why?
			# 1 time: error Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'
			# versuche es abzufangen!
			# sometimes inv_Q = only NA'S !!!! or gives as output empty model
			inv_Q<- pseudoinverse(Q) 
		}
		
		if(  all(is.na(inv_Q)) ) {
			# exit the loop
			diff<- tol; flag.failed<- TRUE
		}else{
			flag.failed<- FALSE
			#di=-((H*H'+F)\g);
			di <- (-1)* (inv_Q %*% g )
			u=u+di;
			#diff<- sqrt( sum( (g)^2 ) )
			diff<- sqrt( sum( (u-uold)^2 ) )
		}
	}
	
	du=d*u;
	Adu=t(A)%*%du;
	pp=Adu-en; pp[pp<0]<-0 # positive
	np=-Adu-en; np[np<0]<- 0 # negative 
	
	w=1/epsi[1]*(pp-np);
	gamma=-(1/epsi[1])*sum(du);
	# adapt w: skip all 0!  
  xind<- which(!(w==0))
  w<-w[xind ,]
  # if only one w !=0,  give the name manually
  if (length(xind) == 1 ) names(w)<- colnames(A)[xind]
	
	if (!flag.failed) {	List<-list(w=w,gamma=gamma,xind=xind,epsi=epsi[1],delta=delta, tol=tol, iter=iter)
	} else 	List<- NULL

	return(List)
}

`.correctness` <-
function(AA,dd,w,gamma){
	# separating hyperplane: AA *w = gamma ! 
	p=sign(AA %*%w-gamma)
	corr=sum(p==dd)/nrow(AA)*100;
	return(corr)
}



`.EstNuLong` <-
function( C,d) {
	
	m<-dim(C)[1]
	n<-dim(C)[2]
	e<- rep(1, m)
	
	# H = [] erzeugen einer neuen Matrix
	# H=[C -e]; 
	H <- cbind(C, (-1)*e) 
	if ( m<201){
	H2=H;	d2=as.matrix(as.numeric(as.character(d)));
	}else{
		# rand(n,m) erzeugt eine  nxm - Matrix mit gleichförmig verteilten Zufallszahlen zwischen 0 und 1
		#r=rand(m,1);
		r<- matrix(runif(m),m, 1)
		# [s1,s2]=sort(r); ### s1, sorted vector, s2 - order(r); 
		s1<-sort(r)
		s2<-order(r)
		# take the first 200 clones from random order 
		H2=H[s2[1:200],];
		d2=as.matrix(d[s2[1:200]]);
	}

	lamda<-1;
	#[vu,u]=eig(H2 %*% t(H2));
	# [V,D] = eig(A) produces matrices of eigenvalues (D) and eigenvectors (V) of matrix A, so that A*V = V*D. 
	tt<- eigen(H2 %*% t(H2))
	vu <- tt$vectors
	u<- tt$values
	p=length(u)
	
	yt=t(d2)%*% vu;
	lamdaO=lamda+1;
	
	cnt<-0
	while ( (abs(lamdaO-lamda)>10e-4) & (cnt<100)) {
		cnt=cnt+1;
		nu1<-0; pr<-0; ee<-0; waw<-0;
		lamdaO<-lamda;
		for (i in 1:p){
			nu1= nu1 + lamda/(u[i]+lamda);
			pr= pr + u[i]/(u[i]+lamda)^2;
			ee= ee + u[i]*yt[i]^2/(u[i]+lamda)^3;
			waw= waw + lamda^2*yt[i]^2/(u[i]+lamda)^2;
		} 
		lamda=nu1*ee/(pr*waw);
	}
	value <- lamda
	
	if (cnt==100) value<-1 
	return(value)
}

`.EstNuShort` <-
function(C,d){
		# easy way to estimate nu if not specified by the user
	#	value = 1/(sum(sum(C.^2))/size(C,2));
	# size = dim ; 
	# C^2 = C %*% C  (matrix multiplication), 
	# .^ Potenzierung C.^2 = C * C = C^2 (elementerweise !)
		value = 1/(sum(sum(C^2))/ncol(C));    
		return( value)
}

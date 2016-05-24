splithalfCP <-
function(X,n,m,p,r,centopt,normopt,scaleopt,addanal,conv,maxit,ort1,ort2,ort3,laba,labb,labc){

cat("This procedure performs a SPLIT-HALF analysis on X",fill=TRUE)
cat("NOTE:",fill=TRUE)
cat("In SPLIT-HALF analysis, the A-mode is taken as 'replication mode'",fill=TRUE)
cat("(which means that A-mode entities are considered a random sample",fill=TRUE)
cat("if this does not make sense, you should rearrange your data so that the A-mode is a replication mode)",fill=TRUE)
cat("The splitting into two halves can be done randomly (default), or into odd vs. even sequence numbers",fill=TRUE)

narg=nargs()
if (narg<15){
	laba=paste("a",1:n,sep="")
	labb=paste("b",1:m,sep="")
	labc=paste("c",1:p,sep="")
}

X=as.matrix(X)
cat("If you prefer odd/even split, specify '1':",fill=TRUE)
nsplit=scan("",n=1)
if (length(nsplit)==0){
	nsplit=0
}
if (nsplit==1){		#wi=[1:2:n 2:2:n]
	check=0
	wi=c(seq(1,n,2), seq(2,n,2))
	cat("Splitting has been done into odd vs. even sequence numbers",fill=TRUE)	
} else{
	# create random splits
	w=matrix(runif(n*1,0,1),n,1)
	wi=ord(w)$a
	cat("Splitting has been done randomly",fill=TRUE)	
	}
n1=ceiling(n/2)
n2=n-n1
X1=X[wi[1:n1],]
X2=X[wi[(n1+1):n],]

cat("You will now enter the split half procedure with the options specified so far",fill=TRUE)
cat("However, you may want to modify certain choices here (rather than rerunning the full Candecomp/Parafac)",fill=TRUE)
cat("Do you want to modify certain choices? If so, specify '1':",fill=TRUE)
ccc=scan("",n=1)
if (length(ccc)==0){
	ccc=0
}
if (ccc==1){
	cat("Centering:",fill=TRUE)
	cat(" 0 = none (default)", fill=TRUE)
	cat(" 1= across A-mode",fill=TRUE)
	cat(" 2= across B-mode",fill=TRUE)
	cat(" 3= across C-mode",fill=TRUE)
	cat(" 12= across A-mode and across B-mode",fill=TRUE)
	cat(" 13= across A-mode and across C-mode",fill=TRUE)
	cat(" 23= across B-mode and across C-mode",fill=TRUE)
	cat(paste("Your centering choice was",centopt),fill=TRUE)	
	cat("If you want to change it, specify '1':",fill=TRUE)
	cc=scan("",n=1)
	if (length(cc)==0){
		cc=0
	}
	if (cc==1){
		check=0
		while(check==0){
			cat("Specify centering option:",fill=TRUE)
			centopt=scan("",n=1)
			if (length(centopt)==0){
				centopt=0
			}
			if((centopt!=0) & ((centopt!=1) & (centopt!=2) & (centopt!=3) & (centopt!=12) & (centopt!=13) & (centopt!=23))){ 
				cat(" ",fill=TRUE)
				cat("Error! Make a proper choice for centering the data!")
				cat(" ",fill=TRUE)
			} else{
				check=1
			}
		}
	}

	cat("Normalizing:",fill=TRUE)
	cat(" 0= none (default)", fill=TRUE)
	cat(" 1= within A-mode", fill=TRUE)
	cat(" 2= within B-mode",fill=TRUE)
	cat(" 3= within C-mode",fill=TRUE)
	cat(paste("Your normalizing choice was",normopt),fill=TRUE)
	cat("If you want to change it, specify '1':",fill=TRUE)
	cc=scan("",n=1)
	if (length(cc)==0){
		cc=0
	}
	if (cc==1){
		check=0
		while(check==0){
			cat("Specify normalizing option",fill=TRUE)
			normopt=scan("",n=1)
			if (length(normopt)==0){
				normopt=0
			}
			if ((normopt!=0) & ((normopt!=1) & (normopt!=2) & (normopt!=3))){
				cat(" ",fill=TRUE)
				cat("Error! Make a proper choice for normalizing the data!")
				cat(" ",fill=TRUE)
			} else{
				check=1
			}
		}
	}

  	cat(paste("Number of components for A, B and C was:",r),fill=TRUE)
	cat("If you want to change it, specify '1':",fill=TRUE)
	cc=scan("",n=1)
	if (length(cc)==0){
		cc=0
	}
	if (cc==1){
		cat("How many components do you want to use?",fill=TRUE)
		check=0
		while (check==0){
			r=scan("",n=1)
			if (length(r)==0){
				r=0
			}
			if ((r==0) | ((floor(r)-r)!=0)){
				cat(" ",fill=TRUE)
				cat("Error! How many components do you want to use?",fill=TRUE)
				cat(" ",fill=TRUE)
			} else{
				check=1
			}
		}
	}
 
 	#constraints
 	cat("Constraints:",fill=TRUE)
	cat("1 = No constraints",fill=TRUE)
	cat("2 = Orthogonality constraints",fill=TRUE)
	cat("3 = Zero correlations constraints",fill=TRUE)
 	cat(paste("Constraints on A, B and C were: ",ort1,", ",ort2,", ",ort3,sep=""),fill=TRUE)
	cat("If you want to change them, specify '1':",fill=TRUE)
	cc=scan("",n=1)
	if (length(cc)==0){
		cc=0
	}
	if (cc==1){
		cat("Do you want to use constraints? If so, enter '1':",fill=TRUE)
		c=scan("",n=1)
		if (length(c)==0){			
			c=0
		}
		if (c!=1){
			ort1=1
			ort2=1
			ort3=1
			cat("No constraints imposed",fill=TRUE)
		}
		if (c==1){
			check=0
			while (check==0){
				check=1
				cat("Digit:",fill=TRUE)
				cat("1 = No constraints (default)",fill=TRUE)
				cat("2 = Orthogonality constraints",fill=TRUE)
				cat("3 = Zero correlations constraints",fill=TRUE)
				cat("Specify the A-mode constraints:",fill=TRUE)
				ort1=scan("",n=1)
				cat("Specify the B-mode constraints:",fill=TRUE)
				ort2=scan("",n=1)
				cat("Specify the C-mode constraints:",fill=TRUE)
				ort3=scan("",n=1)

				if (length(ort1)==0){
					ort1=1
				}
				if (length(ort2)==0){
					ort2=1
				}
				if (length(ort3)==0){
					ort3=1
				}
				if ((ort1!=1) & (ort1!=2) & (ort1!=3)){
					ort1=1
					cat(" ",fill=TRUE)
					cat("Error! Incorrect constraint! A will not be constrained!", fill=TRUE)
					cat(" ",fill=TRUE)
				}
				if ((ort2!=1) & (ort2!=2) & (ort2!=3)){
					ort2=1
					cat(" ",fill=TRUE)
					cat("Error! Incorrect constraint! B will not be constrained!", fill=TRUE)
					cat(" ",fill=TRUE)
				}
				if ((ort3!=1) & (ort3!=2) & (ort3!=3)){
					ort3=1
					cat(" ",fill=TRUE)
					cat("Error! Incorrect constraint! C will not be constrained!", fill=TRUE)
					cat(" ",fill=TRUE)
				}
				if ((ort1==2) & (r>n)){
					cat(" ",fill=TRUE)
					cat(paste("Warning: due to the constraint on A, you cannot use more than", n, "components",sep=" "),fill=TRUE)
					cat(paste(n, "components will be used",sep=" "),fill=TRUE)
					cat(" ",fill=TRUE)
					r=n
					check=0
				}
				if ((ort2==2) & (r>m)){
					cat(" ",fill=TRUE)
					cat(paste("Warning: due to the constraint on A, you cannot use more than", m, "components",sep=" "),fill=TRUE)
					cat(paste(m, "components will be used",sep=" "),fill=TRUE)
					cat(" ",fill=TRUE)
					r=m
					check=0
				}
				if ((ort3==2) & (r>p)){
					cat(" ",fill=TRUE)
					cat(paste("Warning: due to the constraint on A, you cannot use more than", p, "components",sep=" "),fill=TRUE)
					cat(paste(p, "components will be used",sep=" "),fill=TRUE)
					cat(" ",fill=TRUE)
					r=p
					check=0
				}
				if ((ort1==3) & (r>n-1)){
					cat(" ",fill=TRUE)
					cat(paste("Warning: due to the constraint on A, you cannot use more than", n-1, "components",sep=" "),fill=TRUE)
					cat(paste(n-1, "components will be used",sep=" "),fill=TRUE)
					cat(" ",fill=TRUE)
					r=n-1
					check=0
				}
				if ((ort2==3) & (r>m-1)){
					cat(" ",fill=TRUE)
					cat(paste("Warning: due to the constraint on A, you cannot use more than", m-1, "components",sep=" "),fill=TRUE)
					cat(paste(m-1, "components will be used",sep=" "),fill=TRUE)
					cat(" ",fill=TRUE)
					r=m-1
					check=0
				}
				if ((ort3==3) & (r>p-1)){
					cat(" ",fill=TRUE)
					cat(paste("Warning: due to the constraint on A, you cannot use more than", p-1, "components",sep=" "),fill=TRUE)
					cat(paste(p-1, "components will be used",sep=" "),fill=TRUE)
					cat(" ",fill=TRUE)
					r=p-1
					check=0
				}
			}
		}
	}

	cat("Scaling:",fill=TRUE)
	cat("0 = None ",fill=TRUE)
	cat("1 = B- and C-modes (scaling compensated in A) ",fill=TRUE)
	cat("2 = A- and C-modes (scaling compensated in B) ",fill=TRUE)
	cat("3 = A- and B-modes (scaling compensated in C) ",fill=TRUE)
	cat(paste("Your choice was",scaleopt),fill=TRUE)
	cat("If you want to change it, specify '1':",fill=TRUE)
	cc=scan("",n=1)
	if (length(cc)==0){
		cc=0
	}
	if (cc==1){
		check1=0
		while (check1==0){
			cat("What modes do you want to scale?",fill=TRUE)
			cat("0 = None (default)",fill=TRUE)
			cat("1 = B- and C-modes (scaling compensated in A) ",fill=TRUE)
			cat("2 = A- and C-modes (scaling compensated in B) ",fill=TRUE)
			cat("3 = A- and B-modes (scaling compensated in C) ",fill=TRUE)
			scaleopt=scan("",n=1)
			if (length(scaleopt)==0){
				scaleopt=0
			}
			if ((scaleopt==0) | (scaleopt==1) | (scaleopt==2) | (scaleopt==3)){
				check1=1
			} else{
				cat(" ",fill=TRUE)
				cat("Error! Select a proper number!")
				cat(" ",fill=TRUE)
			}
		}
	}
	
	cat(paste("Your choice was to use a convergence criterion equal to",conv),fill=TRUE)
	cat("If you want to change it, specify '1':",fill=TRUE)
	cc=scan("",n=1)
	if (length(cc)==0){
		cc=0
	}
	if (cc==1){
		cat("Specify convergence criterion (default=1e-6)", fill=TRUE)
		conv=scan("",n=1)
		if (length(conv)==0){  
		conv=1e-6
		}
	}

	cat(paste("Your choice was to use",addanal,"additional random starts in the analysis"),fill=TRUE)
	cat("If you want to change it, specify '1':",fill=TRUE)
	cc=scan("",n=1)
	if (length(cc)==0){
		cc=0
	}
	if (cc==1){
		check=0
		cat("How many additional runs do you want to use?",fill=TRUE)
		while (check==0){
			addanal=scan("",n=1)
			if (length(addanal)==0){
				addanal=0
				check=1
			}
			if ((floor(addanal)-addanal)!=0){
				cat(" ",fill=TRUE)
				cat("Error! How many additional runs do you want to use?",fill=TRUE)
				cat(" ",fill=TRUE)
			} else{
				check=1
			}
		}
	}
	
	cat(paste("Your choice was to allow a maximum number of iterations equal to",maxit),fill=TRUE)
	cat("If you want to change it, specify '1':",fill=TRUE)
	cc=scan("",n=1)
	if (length(cc)==0){
		cc=0
	}
	if (cc==1){
		check=0
		cat("Specify the maximum number of iterations you allow (default=10000).",fill=TRUE)
		while (check==0){
			maxit=scan("",n=1)
			if (length(maxit)==0){
				maxit=10000
				check=1
			}
			if ((floor(maxit)-maxit)!=0){
				cat(" ",fill=TRUE)
				cat("Error! Specify the maximum number of iterations you allow (default=10000)",fill=TRUE)
				cat(" ",fill=TRUE)
			} else{
				check=1
			}
		}
	}
}

cat("Analysis of FULL DATA",fill=TRUE)
# Full data analysis (reference)
# Full data Preprocessing 
cc=centopt
if ((cc==1) | (cc==12) | (cc==13)){
	X=cent3(X,n,m,p,1)
	cat("X has been centered across A-mode",fill=TRUE)
}
if ((cc==2) | (cc==12) | (cc==23)){
	X=cent3(X,n,m,p,2)
	cat("X has been centered across B-mode",fill=TRUE)
}
if ((cc==3) | (cc==13) | (cc==23)){
	X=cent3(X,n,m,p,3)
	cat("X has been centered across C-mode",fill=TRUE)
} else{
	cat("X has not been centered",fill=TRUE)
}
cc=normopt
if (cc==1){
	X=norm3(X,n,m,p,1)
	cat("X has been normalized within A-mode",fill=TRUE)
}
if (cc==2){
	X=norm3(X,n,m,p,2)
	cat("X has been normalized within B-mode",fill=TRUE)
}
if (cc==3){
	X=norm3(X,n,m,p,3)
	cat("X has been normalized within C-mode",fill=TRUE)
}
if ((cc!=1) & (cc!=2) & (cc!=3)){
	cat("X has not been normalized",fill=TRUE)
}
# Full data analysis 

start=0
if (n>=r){
	A=orth(matrix(runif(n*r,0,1),nrow=n)-0.5)
} else{
	A=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
	A=A[1:n,]
}
if (m>=r){
	B=orth(matrix(runif(m*r,0,1),nrow=m)-0.5)
} else{
	B=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
	B=B[1:m,]
}
if (p>=r){
	C=orth(matrix(runif(p*r,0,1),nrow=p)-0.5)
} else{
	C=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
	C=C[1:p,]
}
Para=CPfunc(X,n,m,p,r,ort1,ort2,ort3,start,conv,maxit,A,B,C)		
A=Para$A
B=Para$B
C=Para$C
f=Para$f
fp=Para$fp
func=vector("numeric",length=1+addanal)
names(func)=paste("Start n.",1:(1+addanal),sep="")
func=Para$fp
print(func)
for (run in 1:addanal){
    cat(paste("Run no.",run+1,sep=" "),fill=TRUE)
	start=1
	Parab=CPfunc(X,n,m,p,r,ort1,ort2,ort3,start,conv,maxit,A,B,C)   
	func[run+1]=Parab$fp
	if (Parab$fp>1.0001*Para$fp){		# if fit more than .01% better is found, replace solution
		A=Parab$A
		B=Parab$B
		C=Parab$C
		f=Parab$f
		fp=Parab$fp
    }
}
if (addanal>=1){
	cat("Fit % values from all runs:",fill=TRUE)
	print(round(t(func), digits=2))
}

if ((scaleopt==1) | (scaleopt==2) | (scaleopt==3)){
	RNsol=renormsolCP(A,B,C,scaleopt)
	A=RNsol$A
	B=RNsol$B
	C=RNsol$C
}

Afull=A
Bfull=B
Cfull=C

# Split 1 data analysis (reference)
# Split 1 Preprocessing 
cat("Analysis of SPLIT 1",fill=TRUE)
n_orig=n
n=n1
X_orig=X
X=X1
cc=centopt
if ((cc==1) | (cc==12) | (cc==13)){
	X=cent3(X,n,m,p,1)
	cat("X has been centered across A-mode",fill=TRUE)
}
if ((cc==2) | (cc==12) | (cc==23)){
	X=cent3(X,n,m,p,2)
	cat("X has been centered across B-mode",fill=TRUE)
}
if ((cc==3) | (cc==13) | (cc==23)){
	X=cent3(X,n,m,p,3)
	cat("X has been centered across C-mode",fill=TRUE)
}
if ((cc!=1) & (cc!=2) & (cc!=3) & (cc!=12) & (cc!=13) & (cc!=23)){
	cat("X has not been centered",fill=TRUE)
}
cc=normopt
if (cc==1 ){
	X=norm3(X,n,m,p,1)
	cat("X has been normalized within A-mode",fill=TRUE)
}
if (cc==2){
	X=norm3(X,n,m,p,2)
	cat("X has been normalized within B-mode",fill=TRUE)
}
if (cc==3){
	X=norm3(X,n,m,p,3)
	cat("X has been normalized within C-mode",fill=TRUE)
}
if ((cc!=1) & (cc!=2) & (cc!=3)){
	cat("X has not been normalized",fill=TRUE)
}
Xs1=X
# Split 1 analysis 
start=0
if (n>=r){
	A=orth(matrix(runif(n*r,0,1),nrow=n)-0.5)
} else{
	A=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
	A=A[1:n,]
}
if (m>=r){
	B=orth(matrix(runif(m*r,0,1),nrow=m)-0.5)
} else{
	B=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
	B=B[1:m,]
}
if (p>=r){
	C=orth(matrix(runif(p*r,0,1),nrow=p)-0.5)
} else{
	C=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
	C=C[1:p,]
}
Para=CPfunc(X,n,m,p,r,ort1,ort2,ort3,start,conv,maxit,A,B,C)		
A=Para$A
B=Para$B
C=Para$C
f=Para$f
fp=Para$fp
func=vector("numeric",length=1+addanal)
names(func)=paste("Start n.",1:(1+addanal),sep="")
func=Para$fp
for (run in 1:addanal){
    cat(paste("Run no.",run+1,sep=" "),fill=TRUE)
	start=1
	Parab=CPfunc(X,n,m,p,r,ort1,ort2,ort3,start,conv,maxit,A,B,C)   
	func[run+1]=Parab$fp
	if (Parab$fp>1.0001*Para$fp){		# if fit more than .01% better is found, replace solution
		A=Parab$A
		B=Parab$B
		C=Parab$C
		f=Parab$f
		fp=Parab$fp
	}
}
if (addanal>=1){
	cat("Fit % values from all runs:",fill=TRUE)
	print(round(t(func), digits=2))
}
# normalize solution
if ((scaleopt==1) | (scaleopt==2) | (scaleopt==3)){
	RNsol=renormsolCP(A,B,C,scaleopt)
	A=RNsol$A
	B=RNsol$B
	C=RNsol$C
}

As1=A
Bs1=B
Cs1=C

# Split 2 data analysis (reference)
# Split 2 Preprocessing 
cat("Analysis of SPLIT 2",fill=TRUE)
n=n2
X=X2
cc=centopt
if ((cc==1) | (cc==12) | (cc==13)){
	X=cent3(X,n,m,p,1)
	cat("X has been centered across A-mode",fill=TRUE)
}
if ((cc==2) | (cc==12) | (cc==23)){
	X=cent3(X,n,m,p,2)
	cat("X has been centered across B-mode",fill=TRUE)
}
if ((cc==3) | (cc==13) | (cc==23)){
	X=cent3(X,n,m,p,3)
	cat("X has been centered across C-mode",fill=TRUE)
}
if ((cc!=1) & (cc!=2) & (cc!=3) & (cc!=12) & (cc!=13) & (cc!=23)){
	cat("X has not been centered",fill=TRUE)
}
cc=normopt
if (cc==1){
	X=norm3(X,n,m,p,1)
	cat("X has been normalized within A-mode",fill=TRUE)
}
if (cc==2){
	X=norm3(X,n,m,p,2)
	cat("X has been normalized within B-mode",fill=TRUE)
}
if (cc==3){
	X=norm3(X,n,m,p,3)
	cat("X has been normalized within C-mode",fill=TRUE)
}
if ((cc!=1) & (cc!=2) & (cc!=3)){
	cat("X has not been normalized",fill=TRUE)
}
Xs2=X
# Split 2 analysis 
start=0
if (n>=r){
	A=orth(matrix(runif(n*r,0,1),nrow=n)-0.5)
} else{
	A=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
	A=A[1:n,]
}
if (m>=r){
	B=orth(matrix(runif(m*r,0,1),nrow=m)-0.5)
} else{
	B=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
	B=B[1:m,]
}
if (p>=r){
	C=orth(matrix(runif(p*r,0,1),nrow=p)-0.5)
} else{
	C=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
	C=C[1:p,]
}
Para=CPfunc(X,n,m,p,r,ort1,ort2,ort3,start,conv,maxit,A,B,C)		
A=Para$A
B=Para$B
C=Para$C
f=Para$f
fp=Para$fp
func=vector("numeric",length=1+addanal)
names(func)=paste("Start n.",1:(1+addanal),sep="")
func=Para$fp
for (run in 1:addanal){
    cat(paste("Run no.",run+1,sep=" "),fill=TRUE)
	start=1
	Parab=CPfunc(X,n,m,p,r,ort1,ort2,ort3,start,conv,maxit,A,B,C)   
	func[run+1]=Parab$fp
	if (Parab$fp>1.0001*Para$fp){		# if fit more than .01% better is found, replace solution
		A=Parab$A
		B=Parab$B
		C=Parab$C
		f=Parab$f
		fp=Parab$fp
	}
}
if (addanal>=1){
	cat("Fit % values from all runs:",fill=TRUE)
	print(round(t(func), digits=2))
}
# normalize solution
if ((scaleopt==1) | (scaleopt==2) | (scaleopt==3)){
	RNsol=renormsolCP(A,B,C,scaleopt)
	A=RNsol$A
	B=RNsol$B
	C=RNsol$C
}

As2=A
Bs2=B
Cs2=C

n=n_orig
X=X_orig

# Compute stability indices:
# Split 1
# permutation
if (r>1){
	coeff=rep(0,factorial(r))
	permutation=perms(r)
	for (j in 1:(nrow(permutation))){
		Bbperm=Bs1[,permutation[j,]]
		Cbperm=Cs1[,permutation[j,]]
		for (r1 in 1:r){
			coeff[j]=coeff[j]+abs(phi(Bfull[,r1],Bbperm[,r1]))*abs(phi(Cfull[,r1],Cbperm[,r1]))
		}
	}
	As1=As1[,permutation[which(coeff==max(coeff)),]]
	Bs1=Bs1[,permutation[which(coeff==max(coeff)),]]
	Cs1=Cs1[,permutation[which(coeff==max(coeff)),]]
}
# reflection
for (r1 in 1:r){
	if (phi(Bs1[,r1],Bfull[,r1])<0){
		Bs1[,r1]=-Bs1[,r1]
		Cs1[,r1]=-Cs1[,r1]
	}
	if (phi(Cs1[,r1],Cfull[,r1])<0){
		Cs1[,r1]=-Cs1[,r1]
		As1[,r1]=-As1[,r1]
	}
}
# Split 2
# permutation
if (r>1){
	coeff=rep(0,factorial(r))
	permutation=perms(r)
	for (j in 1:(nrow(permutation))){
		Bbperm=Bs2[,permutation[j,]]
		Cbperm=Cs2[,permutation[j,]]
		for (r1 in 1:r){
			coeff[j]=coeff[j]+abs(phi(Bfull[,r1],Bbperm[,r1]))*abs(phi(Cfull[,r1],Cbperm[,r1]))
		}
	}
	As2=As2[,permutation[which(coeff==max(coeff)),]]
	Bs2=Bs2[,permutation[which(coeff==max(coeff)),]]
	Cs2=Cs2[,permutation[which(coeff==max(coeff)),]]
}
# reflection
for (r1 in 1:r){
	if (phi(Bs2[,r1],Bfull[,r1])<0){
		Bs2[,r1]=-Bs2[,r1]
		Cs2[,r1]=-Cs2[,r1]
		}
	if (phi(Cs2[,r1],Cfull[,r1])<0){
		Cs2[,r1]=-Cs2[,r1]
		As2[,r1]=-As2[,r1]
	}
}

labComp=paste("Comp.",1:r,sep="")

cat("RESULTS stability analysis",fill=TRUE)
cat("Both splits are analyzed and permuted/scaled towards solution for full data",fill=TRUE)
cat("One can compare complete outcomes for different splits",fill=TRUE)
cat("or inspect correlation/congruence coefficients between corresponding columns of component matrices",fill=TRUE)
cat("Split-half solutions for B",fill=TRUE)
cat("B's next to each other, separated by column of 0's",fill=TRUE)
X=round(cbind(Bs1,0,Bs2),digits=2)
rownames(X)=labb
colnames(X)=c(labComp,"-",labComp)
print(X)
cat("Split-half solutions for C",fill=TRUE)
cat("C's next to each other, separated by column of 0's",fill=TRUE)
X=round(cbind(Cs1,0,Cs2),digits=2)
rownames(X)=labc
colnames(X)=c(labComp,"-",labComp)
print(X)
cat("Congruences for A in splits and in appropriate part of Afull",fill=TRUE)
X=round(cbind(diag(phi(Afull[wi[1:n1],],As1)),diag(phi(Afull[wi[(n1+1):n],],As2))),digits=2)
rownames(X)=labComp
colnames(X)=c("SPL1","SPL2")
print(X)
cat("Correlations for A in splits and in appropriate part of Afull",fill=TRUE)
X=round(cbind(diag(phi(Cc(Afull[wi[1:n1],]),Cc(As1))),diag(phi(Cc(Afull[wi[(n1+1):n],]),Cc(As2)))),digits=2)
rownames(X)=labComp
colnames(X)=c("SPL1","SPL2")
print(X)
cat("Congruence values for B-mode component matrix",fill=TRUE)
X=round(diag(phi(Bs1,Bs2)),digits=2)
names(X)=c(labComp)
print(X)
cat("Congruence values for C-mode component matrix",fill=TRUE)
X=round(diag(phi(Cs1,Cs2)),digits=2)
names(X)=c(labComp)
print(X)

colnames(As1)=labComp
colnames(As2)=labComp
colnames(Afull)=labComp
colnames(Bs1)=labComp
colnames(Bs2)=labComp
colnames(Bfull)=labComp
colnames(Cs1)=labComp
colnames(Cs2)=labComp
colnames(Cfull)=labComp
rownames(As1)=wi[1:n1]
rownames(As2)=wi[(n1+1):n]
rownames(Afull)=laba
rownames(Bs1)=labb
rownames(Bs2)=labb
rownames(Bfull)=labb
rownames(Cs1)=labc
rownames(Cs2)=labc
rownames(Cfull)=labc

out=list()
out$Afull=Afull
out$As1=As1
out$As2=As2
out$Bfull=Bfull
out$Bs1=Bs1
out$Bs2=Bs2
out$Cfull=Cfull
out$Cs1=Cs1
out$Cs2=Cs2	
}

CP <-
function(data,laba,labb,labc){

cat(" ",fill=TRUE)
cat("WELCOME to the interactive CANDECOMP/PARAFAC analysis program",fill=TRUE)
cat("Warning: If you insert an object of mode CHARACTER when not requested, an error occurs and the program stops!")
cat(" ",fill=TRUE)

# check if the dataset is an array or a dataframe or a matrix
check=0
if (is.matrix(data)==FALSE){
	if (is.array(data)==TRUE){
		n=dim(data)[1]
		m=dim(data)[2]
		p=dim(data)[3]
		X=matrix(data,n,m*p)
	}
	if(is.data.frame(data)==TRUE){
		X=as.matrix(data)
		check=1
	} 
} else if (is.matrix(data)==TRUE){
	X=data
	check=1
}
while(check==1){
	cat(" ",fill=TRUE)
	cat("Specify the number of A-mode entities ",fill=TRUE)
	n=scan("",n=1) 
	cat("Specify the number of B-mode entities ",fill=TRUE)
	m=scan("",n=1)
	cat("Specify the number of C-mode entities ",fill=TRUE)
	p=scan("",n=1)
	if ((length(n)!=0) & (length(m)!=0) & (length(p)!=0)) {
		if ((dim(X)[1]==n) & (dim(X)[2]==m*p)){
			check=0
		} else{
			cat(" ",fill=TRUE)
			cat("Error! The array size does not agree with input data. Please specify ALL size!",fill=TRUE)
			cat(" ",fill=TRUE)
		}
	} else if ((length(n)==0) || (length(m)==0) || (length(p)==0)) {
		cat(" ",fill=TRUE)
		cat("Error! Please specify ALL sizes!",fill=TRUE)
		cat(" ",fill=TRUE)
	}
}

# check if labels are available
cat(" ",fill=TRUE)
narg=nargs()
if (narg>1){
	if (length(laba)!=n){
		cat(" ",fill=TRUE)
		cat("Error! The number of labels for the A-mode entities does not agree with input data. Please insert new labels!",fill=TRUE)
		cat(" ",fill=TRUE)
		narg=1
	}
}
if (narg==1){
# create the vectors of label for the A-, B- and C-modes
    cat("Provide the labels for A-mode entities (e.g., a1 a2 a3 ...), otherwise press 'enter' and default labels will be given. ",fill=TRUE)
	laba=scan(what = "character",n=n)
    if ((length(laba)!=n) & (length(laba)>0)){
		cat(" ",fill=TRUE)
		cat("Error! The number of labels for A-mode entities does not agree with input data. Default labels will be given!",fill=TRUE)
        cat(" ",fill=TRUE)
		laba=paste("a",1:n,sep="")
    }
	if (length(laba)==0){
		laba=paste("a",1:n,sep="")
	}
}	
narg=nargs()
if (narg>1){
	if (length(labb)!=m){
		cat(" ",fill=TRUE)
		cat("Error! The number of labels for B-mode entities does not agree with input data. Please insert new labels!",fill=TRUE)
		cat(" ",fill=TRUE)
		narg=1
	}
}
if (narg==1){	
    cat("Provide the labels for B-mode entities (e.g., b1 b2 b3 ...), otherwise press 'enter' and default labels will be given. ",fill=TRUE)
    labb=scan(what = "character",n=m)
    if ((length(labb)!=m) & (length(labb)>0)){
		cat(" ",fill=TRUE)
		cat("Error! The number of labels for B-mode entities does not agree with input data. Default labels will be given!",fill=TRUE)
        cat(" ",fill=TRUE)
		labb=paste("b",1:m,sep="")
    }
	if (length(labb)==0){
		labb=paste("b",1:m,sep="")
	}
}	
narg=nargs()
if (narg>1){
	if (length(labc)!=p){
		cat(" ",fill=TRUE)
		cat("Error! The number of labels for C-mode entities does not agree with input data. Please insert new labels!",fill=TRUE)
		cat(" ",fill=TRUE)
		narg=1
	}
}
if (narg==1){	
    cat("Provide the labels for C-mode entities (e.g., c1 c2 c3 ...), otherwise press 'enter' and authomatic labels will be given. ",fill=TRUE)
    labc=scan(what = "character",n=p)
    if ((length(labc)!=p) & (length(labc)>0)){
		cat(" ",fill=TRUE)
		cat("Error! The number of labels for C-mode entities does not agree with input data. Default labels will be given!",fill=TRUE)
		cat(" ",fill=TRUE)
		labc=paste("c",1:p,sep="")
    }
	if (length(labc)==0){
		labc=paste("c",1:p,sep="")
	}
}

Xprep=X

# Anova
cat(" ",fill=TRUE)
cat("To see ANOVA results, specify 1:",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
if (c==1){
	threewayanova(Xprep,n,m,p)
}

# Centering
cat(" ",fill=TRUE)
check=1
while (check==1){
	cat("How do you want to center your array?",fill=TRUE)
	cat("0 = none (default)", fill=TRUE)
	cat("1 = across A-mode",fill=TRUE)
	cat("2 = across B-mode",fill=TRUE)
	cat("3 = across C-mode",fill=TRUE)
	cat("12 = across A-mode and across B-mode",fill=TRUE)
	cat("13 = across A-mode and across C-mode",fill=TRUE)
	cat("23 = across B-mode and across C-mode",fill=TRUE)
	prep=scan("",n=1)
	if (length(prep)==0){
		prep=0
	}
	if (prep==0){
		Xprep=Xprep
		cat("Data have not been centered",fill=TRUE)
		check=0
	}
	centopt=prep	
	if ((prep==1) || (prep==12) || (prep==13)){
		Xprep=cent3(Xprep,n,m,p,1)
		cat("Data have been centered across A-mode",fill=TRUE)
		check=0
	}
	if ((prep==2) || (prep==12) || (prep==23)){
		Xprep=cent3(Xprep,n,m,p,2)
		cat("Data have been centered across B-mode",fill=TRUE)
		check=0
	}
	if ((prep==3) || (prep==13) || (prep==23)){
		Xprep=cent3(Xprep,n,m,p,3)
		cat("Data have been centered across C-mode",fill=TRUE)
		check=0
	}
	if ((prep!=0) & (prep!=1) & (prep!=2) & (prep!=3) & (prep!=12) & (prep!=13) & (prep!=23)){
		cat(" ",fill=TRUE)
		cat("Error! Make a proper choice for centering the data",fill=TRUE)
		cat(" ",fill=TRUE)
	}
}

# Normalizing
cat(" ",fill=TRUE)
while (check==0){
    cat("How do you want to normalize your array?",fill=TRUE)
	cat("0 = none (default)", fill=TRUE)
	cat("1 = within A-mode", fill=TRUE)
	cat("2 = within B-mode", fill=TRUE)
	cat("3 = within C-mode", fill=TRUE)
	prep=scan("",n=1)
	if (length(prep)==0){	
		prep=0
	}
	if (prep==0){
		Xprep=Xprep
		cat("Data have not been normalized",fill=TRUE)
		check=1
	}
	normopt=prep
	if (prep==1){ 
		Xprep=norm3(Xprep,n,m,p,1)
		cat("Data have been normalized within A-mode",fill=TRUE)
		check=1
	}
	if (prep==2){
		Xprep=norm3(Xprep,n,m,p,2)
		cat("Data have been normalized within B-mode",fill=TRUE)
		check=1
	}
	if (prep==3){
		Xprep=norm3(Xprep,n,m,p,3)
		cat("Data have been normalized within C-mode",fill=TRUE)
		check=1
	}
	if ((prep!=0) & (prep!=1) & (prep!=2) & (prep!=3) & (prep!=12) & (prep!=13) & (prep!=23)){
		cat(" ",fill=TRUE)
		cat("Error! Make a proper choice for normalizing the data",fill=TRUE)
		cat(" ",fill=TRUE)
	}
}

cat(" ",fill=TRUE)
cat("Note: The preprocessed data are now available in Xprep.",fill=TRUE)
cat("Xprep can be used for analyses outside the Candecomp/Parafac program",fill=TRUE)

cat(" ",fill=TRUE)
cat("To see PCA of MEAN, specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){ 
	c=0
}
if (c==1){
	pcam=1
	while (c==1){
		PCAmean=pcamean(Xprep,n,m,p) 
		cat(" ",fill=TRUE)
		cat("If you want another PCA of mean analysis, specify '1':",fill=TRUE)
		c=scan("",n=1)
		if (length(c)==0){
			c=0
		}
	}
} else{
pcam=0
}
if (pcam==1){
	cat("These matrices can be used for analyses outside the Candecomp/Parafac program and are given in $A1, $A2, $B1, $B2, $C1, $C2",fill=TRUE)
} else{
	PCAmean=list()
	PCAmean$A1=NULL
	PCAmean$A2=NULL
	PCAmean$B1=NULL
	PCAmean$B2=NULL
	PCAmean$C1=NULL
	PCAmean$C2=NULL
}

cat(" ",fill=TRUE)
cat("You are now about to do a Candecomp/Parafac analysis.",fill=TRUE)

cat(" ",fill=TRUE)
cat("You have to specify the dimensionality to use.",fill=TRUE)
cat("To search this, you are advised to run a series of (unconstrained)",fill=TRUE)
cat("Candecomp/Parafac analyses using different dimensionalities",fill=TRUE)
cat("The results can next be used to find useful dimensionalities by the scree test.",fill=TRUE)
cat("To reduce the computation time, note that the maximum allowed number of iterations is 5000,",fill=TRUE)
cat("the convergence criterion is 1e-6 and the rational start via eigendecompositions is considered.",fill=TRUE)
cat("If you want to do this for choosing your dimensionality, specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
if (c==1){
	cat("For the scree test it is needed to indicate the maximum number of dimension you want to study",fill=TRUE)
	cat("Up to how many components do you want to use?",fill=TRUE)
	check0=0
	while (check0==0){
		maxC=scan("",n=1)
		if (length(maxC)==0){
			maxC=0
		}
		if ((maxC==0) | ((floor(maxC)-maxC)!=0)){
			cat(" ",fill=TRUE)
			cat("Error! Up to how many components do you want to use?",fill=TRUE)
			cat(" ",fill=TRUE)
		} else{
			check0=1
		}
	}	
	cat(" ",fill=TRUE)
	out=CPrunsFit(Xprep,n,m,p,maxC)
	cat("       Number of      Fit",fill=TRUE)
	outCP=out[,c(1,4)]
	rownames(outCP)=rep("CP",length=dim(outCP)[1])
	colnames(outCP)=c("    components","    (%) ")
	print(round(outCP,digits=2))
	if (nrow(out)>2){
		cat(" ",fill=TRUE)
		cat("Suggestion based on Convex Hull procedure",fill=TRUE)
		DS=DimSelector(out,n,m,p,1)
	}
	cat("Figure 1 now gives fit values of all solutions plotted against S",fill=TRUE)
	cat("Figure 2 gives fit values of all solutions plotted against number of free parameters",fill=TRUE)
	outplot=out[,c(1,4)]
	CPdimensionalityplot(outplot,n,m,p)
}

cat(" ",fill=TRUE)
cat("You can now do a CANDECOMP/PARAFAC ANALYSIS",fill=TRUE)

cat(" ",fill=TRUE)
cat("How many components do you want to use?",fill=TRUE)
check=0
while (check==0){
	r=scan("",n=1)
	if (length(r)==0){
		r=0
	}	
	if ((r==0) | ((floor(r)-r)!=0)){
		cat(" ",fill=TRUE)
		cat("Error! Please specify the number of components!",fill=TRUE)
		cat(" ",fill=TRUE)
	} else{
		check=1
	}
}

#constraints
cat(" ",fill=TRUE)
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

cat(" ",fill=TRUE)
cat("Specify convergence criterion (default=1e-6)", fill=TRUE)
conv=scan("",n=1)
if (length(conv)==0){  
	conv=1e-6
}

cat(" ",fill=TRUE)
cat("By default, only a rationally started analysis run will be carried out.", fill=TRUE)
cat("To decrease the chance of missing the optimal solution, you may use additional, randomly started runs.", fill=TRUE)
cat("If you want additional runs, specify how many (e.g., 4):", fill=TRUE)
check=0
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

cat(" ",fill=TRUE)
cat("Specify the maximum number of iterations you allow (default=10000).",fill=TRUE)
check=0
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

A=matrix(rnorm(n*r),ncol=r)
B=matrix(rnorm(m*r),ncol=r)
C=matrix(rnorm(p*r),ncol=r)	
# algorithm
cat(paste("Run no.",1,sep=" "),fill=TRUE)
PAR=CPfunc(Xprep,n,m,p,r,ort1,ort2,ort3,0,conv,maxit,A,B,C)		
A=PAR$A
B=PAR$B
C=PAR$C
f=PAR$f
iter=PAR$iter
fp=PAR$fp
tripcos=PAR$tripcos
iter=PAR$iter
cput=PAR$cputime
fp=PAR$fp
f=PAR$f
funcF=vector("numeric",length=1+addanal)
funcV=vector("numeric",length=1+addanal)
cputime=vector("numeric",length=1+addanal)
niter=vector("numeric",length=1+addanal)
funcF[1]=fp
funcV[1]=f
cputime[1]=cput
niter[1]=iter
namev=paste("Start n.",1:(1+addanal),sep="")
names(funcF)=namev
names(funcV)=namev
names(cputime)=namev
names(niter)=namev
if (addanal>=1){
	for (run in 1:addanal){
		cat(paste("Run no.",run+1,sep=" "),fill=TRUE)
		start=1
		PAR2=CPfunc(Xprep,n,m,p,r,ort1,ort2,ort3,start,conv,maxit,A,B,C)   
		funcF[run+1]=PAR2$fp
		funcV[run+1]=PAR2$f
		cputime[run+1]=PAR2$cputime
		niter[run+1]=PAR2$iter
		if (PAR2$fp>1.0001*PAR$fp){		# if fit more than .01% better is found, replace solution
			A=PAR2$A
			B=PAR2$B
			C=PAR2$C
			f=PAR2$f
			fp=PAR2$fp
			tripcos=PAR2$tripcos
		}
	}
	cat(" ",fill=TRUE)
	cat("Fit (%) values from all runs:",fill=TRUE)
	print(round(funcF,digits=2))
}
names(fp)="Fit (%)"

cat(" ",fill=TRUE)
cat(paste("Candecomp/Parafac analysis with ",r," components, gave a fit of ",round(fp,digits=2), "%"),fill=TRUE) 
if (r>1){
	cat("Simple check on degeneracy: inspect matrix of triple congruences",fill=TRUE)
	TC=phi(A,A)*phi(B,B)*phi(C,C)
	rownames(TC)=paste("Comp.",1:r,sep="")
	colnames(TC)=rownames(TC)
	print(round(TC,digits=4))
}
H=matrix(0,r,r^2)       # superidentity 3-way array
for (ii in 1:r){
	H[ii,(ii-1)*r+ii]=1
}

cat(" ",fill=TRUE)
cat("It is sometimes useful to SCALE solution, e.g., scale two matrices so that",fill=TRUE)
cat("they have unit sum of squares compensating this scale in remaining matrix.",fill=TRUE)
cat("If you want to scale components, specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
scalemode=0
if (c==1){
	check1=0
	while (check1==0){
		cat("What modes do you want to scale?",fill=TRUE)
		cat("1 = B- and C-modes (scaling compensated in A) ",fill=TRUE)
		cat("2 = A- and C-modes (scaling compensated in B) ",fill=TRUE)
		cat("3 = A- and B-modes (scaling compensated in C) ",fill=TRUE)
		scalemode=scan("",n=1)
		if (length(scalemode)==0){
			scalemode=0
		}
		if ((scalemode==1) | (scalemode==2) | (scalemode==3)){
			check1=1
			RNsol=renormsolCP(A,B,C,scalemode)
			A=RNsol$A
			B=RNsol$B
			C=RNsol$C
		} else{
			cat(" ",fill=TRUE)
			cat("Error! Select a proper number!")
			cat(" ",fill=TRUE)
		}
	}
}

labComp=paste("Comp.",1:r,sep="")
rownames(A)=laba
rownames(B)=labb
rownames(C)=labc
colnames(A)=labComp
colnames(B)=labComp
colnames(C)=labComp

cat(" ",fill=TRUE)
cat("Solution for A, B and C in detail.",fill=TRUE)
cat("A",fill=TRUE)
print(round(A, digits=2))
cat("B",fill=TRUE)
print(round(B, digits=2))
cat("C",fill=TRUE)
print(round(C, digits=2))

# manual permutation and reflection:
cat(" ",fill=TRUE)
cat("You can now manually PERMUTE and REFLECT columns of solution",fill=TRUE)
cat("If you want to reflect/permute columns, specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
disp=1
while (c==1){
	if (disp==1){
		cat("Permuted and/or reflected solution for A, B and C",fill=TRUE)
		cat("A",fill=TRUE)
		print(round(A, digits=2))
		cat("B",fill=TRUE)
		print(round(B, digits=2))
		cat("C",fill=TRUE)
		print(round(C, digits=2))
	}	
	if (r>1){
		cat(" ",fill=TRUE)
		cat("Give a vector with new order of columns of A, B and C  (e.g., 3 1 4 2 ..)",fill=TRUE)
		pp=scan("",n=r)
		if (length(pp)==0){
			pp=c(1:r)
			cat(" ",fill=TRUE)
			cat("The columns of A, B and C will not be permuted", fill=TRUE)
			cat(" ",fill=TRUE)
		} else{
			if ((length(pp)!=r) | (sum((sort(pp)-1:r)^2)>0)){	
				cat(" ",fill=TRUE)
				cat("Error! Incorrect permutation! The columns of A, B and C will not be permuted!",fill=TRUE)
				cat(" ",fill=TRUE)
				pp=c(1:r)
			}
		}	
	} else{
		pp=1
	}
	cat("Give a vector for reflection of columns of A (e.g., 1 -1 -1 1 ..)", fill=TRUE)
	tau=scan("",n=r)
	if (length(tau)==0){
		tau=array(1,r)
		cat(" ",fill=TRUE)
		cat("Warning: the columns of A will not be reflected",fill=TRUE)
		cat(" ",fill=TRUE)
	} else{
		if ((sum(tau^2)!=r) | (sum((tau^2-1)^2)!=0)){
			cat(" ",fill=TRUE)
			cat("Error! Incorrect reflection! The columns of A will not be reflected!",fill=TRUE)
			cat(" ",fill=TRUE)
			tau=array(1,r)
		}
	}
	A=A%*%diag(tau,nrow=r)
	B=B%*%diag(tau,nrow=r)		
	A=as.matrix(A[,pp])
	B=as.matrix(B[,pp])
	C=as.matrix(C[,pp])
	colnames(A)=labComp
	colnames(B)=labComp
	colnames(C)=labComp

	cat("B",fill=TRUE)
	print(round(B, digits=2))
	cat("C",fill=TRUE)
	print(round(C, digits=2))
	cat("Give a vector for reflection of columns of B (e.g., 1 -1 -1 1 ..)",fill=TRUE)
	tau=scan("",n=r)
	if (length(tau)==0){
		tau=array(1,r)
		cat(" ",fill=TRUE)
		cat("Warning: the columns of B will not be reflected",fill=TRUE)
		cat(" ",fill=TRUE)	
	} else{
		if ((sum(tau^2)!=r) | (sum((tau^2-1)^2)!=0)){
			cat(" ",fill=TRUE)
			cat("Error! Incorrect reflection! The columns of B will not be reflected!",fill=TRUE)
			cat(" ",fill=TRUE)
			tau=array(1,r)
		}
	}	
	B=B%*%diag(tau,nrow=r)
	C=C%*%diag(tau,nrow=r)
	colnames(B)=labComp
	colnames(C)=labComp
	cat("C",fill=TRUE)
	print(round(C, digits=2))
	cat("Give a vector for reflection of columns of C (e.g., 1 -1 -1 1 ..)",fill=TRUE)
	cat("Note that this reflection will be compensated in A",fill=TRUE)
	tau=scan("",n=r)
	if (length(tau)==0){
		tau=array(1,r)
		cat(" ",fill=TRUE)
		cat("Warning: the columns of C will not be reflected",fill=TRUE)
		cat(" ",fill=TRUE)
	} else{
		if ((sum(tau^2)!=r) | (sum((tau^2-1)^2)!=0)){
			cat(" ",fill=TRUE)
			cat("Error! Incorrect reflection! The columns of C will not be reflected!",fill=TRUE)
			cat(" ",fill=TRUE)
			tau=array(1,r)
		}
	}
	C=C%*%diag(tau,nrow=r)
	A=A%*%diag(tau,nrow=r)
	colnames(A)=labComp
	colnames(C)=labComp
	cat("C",fill=TRUE)
	print(round(C, digits=2))
	cat(" ",fill=TRUE)
	cat("If you want to further reflect/permute columns/rows, specify '1':",fill=TRUE)
	c=scan("",n=1)
	if (length(c)==0){
		c=0
	}
	if (c!=1){
		cat(" ",fill=TRUE)
		cat("Solution for A, B and C after permutation and reflection.",fill=TRUE)
		cat("A",fill=TRUE)
		print(round(A, digits=2))
		cat("B",fill=TRUE)
		print(round(B, digits=2))
		cat("C",fill=TRUE)
		print(round(C, digits=2))
	}
}	

c=1
while (c==1){
	cat(" ",fill=TRUE)
	cat("If you want to carry out a STABILITY CHECK on current or different solution, specify '1':",fill=TRUE)
	c=scan("",n=1)
	if (length(c)==0){
		c=0
	}
	if (c==1){
		split=1
		SPLIT=splithalfCP(Xprep,n,m,p,r,centopt,normopt,scalemode,addanal,conv,maxit,ort1,ort2,ort3,laba,labb,labc)
	} else{
		split=0
	}	
}
if (split==1){
	cat("Note: for advanced analysis one can further study the results of the splits and of the full data set",fill=TRUE)
	cat("Results for A, B, and C of full data are given in $Afull, $Bfull and $Cfull",fill=TRUE)
	cat("Results for splits are given in $As1, $Bs1, $Cs1 and $As2, $Bs2, $Cs2",fill=TRUE)
} else{
	SPLIT=list()
	SPLIT$Afull=NULL
	SPLIT$As1=NULL
	SPLIT$As2=NULL
	SPLIT$Bfull=NULL
	SPLIT$Bs1=NULL
	SPLIT$Bs2=NULL
	SPLIT$Cfull=NULL
	SPLIT$Cs1=NULL
	SPLIT$Cs2=NULL
}

cat(" ",fill=TRUE)
cat("If you want to carry out a BOOTSTRAP procedure for",fill=TRUE)
cat("computing confidence intervals for the current solution, specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
if (c==1){
	cat("It is assumed that the A-mode is the mode from which you want to resampling",fill=TRUE)
	cat("If this is not the case, first reorder your data",fill=TRUE)
	BOOT=bootstrapCP(Xprep,A,B,C,n,m,p,r,ort1,ort2,ort3,conv,centopt,normopt,scalemode,maxit,laba,labb,labc)	
	# make labels for intervals of columns of B and C
	cat(" ",fill=TRUE)
	cat("Bootstrap confidence intervals for fit percentage",fill=TRUE)
	print(round(BOOT$fpint,digits=2))
	cat("Bootstrap confidence intervals for B, per component next to each other",fill=TRUE)
	print(round(BOOT$Bint,digits=2))
	cat("Bootstrap confidence intervals for C, per component next to each other",fill=TRUE)
	print(round(BOOT$Cint,digits=2))
	cat("Results for the bootstrap confidence intervals are given in $Bint, $Cint, $FITint",fill=TRUE)
} else{
	BOOT=list()
	BOOT$Bint=NULL
	BOOT$Cint=NULL
	BOOT$FITint=NULL
}

# analysis with FITPARTITIONING
# computation of residuals (on Res) and partitioning of fit
cat(" ",fill=TRUE)
cat("If you want to carry out a FITPARTITIONING on current solution, specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
if (c==1){
	cat("Contribution to fit (in %) for all combinations of components",fill=TRUE)
	cat("These contributions can simply be added to get aggregate contributions of in case all components are orthogonal",fill=TRUE)
	cat("If this is not the case, aggregate fit contributions are given separately, below",fill=TRUE)
	CPfit=CPfitpartitioning(Xprep,n,m,p,A,B,C,laba,labb,labc)		
	cat("Relative fit of A-mode entities (in %)", fill=TRUE)
	print(round(CPfit$fitA,digits=2))
	cat("Relative fit of B-mode entities (in %)",fill=TRUE)
	print(round(CPfit$fitB,digits=2))
	cat(" ",fill=TRUE)
	cat("Relative fit of C-mode entities (in %)",fill=TRUE)
	print(round(CPfit$fitC,digits=2))
	cat("Results for the fitpartitioning are given in $fitA, $fitB, $fitC",fill=TRUE)
} else{
	CPfit=list()
	CPfit$fitA=NULL
	CPfit$fitB=NULL
	CPfit$fitC=NULL
}


cat(" ",fill=TRUE)
cat("Press 'return' to conclude the analysis",fill=TRUE)
c=scan("",n=1)
cat("The component matrices are normalized such that A and B have unit sums of squares",fill=TRUE)
cat("To see component matrices and fit and other results, digits: $Xprep, $A, $B, $C, $fit, $fitA, $fitB, $fitC, ...", fill=TRUE)

XprepOut=rarray(Xprep,n,m,p)
dimnames(XprepOut)=list(laba,labb,labc)

out=list()
out$A=A
out$B=B
out$C=C
out$fit=fp
out$tripcos=tripcos
out$fitValues=funcV
out$funcValues=funcF
out$cputime=cputime
out$iter=niter
out$fitA=CPfit$fitA
out$fitB=CPfit$fitB
out$fitC=CPfit$fitC
out$Bint=BOOT$Bint
out$Cint=BOOT$Cint
out$FITint=BOOT$fpint
out$Afull=SPLIT$Afull
out$As1=SPLIT$As1
out$As2=SPLIT$As2
out$Bfull=SPLIT$Bfull
out$Bs1=SPLIT$Bs1
out$Bs2=SPLIT$Bs2
out$Cfull=SPLIT$Cfull
out$Cs1=SPLIT$Cs1
out$Cs2=SPLIT$Cs2	
out$A1=PCAmean$A1
out$B1=PCAmean$B1
out$C1=PCAmean$C1
out$A2=PCAmean$A2
out$B2=PCAmean$B2
out$C2=PCAmean$C2
out$laba=laba
out$labb=labb
out$labc=labc
out$Xprep=XprepOut

return(out)
}

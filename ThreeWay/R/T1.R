T1 <-
function(dati,laba,labb,labc){		

cat(" ",fill=TRUE)
cat("WELCOME to the interactive TUCKER1 analysis program",fill=TRUE)
cat("Warning: If you insert an object of mode CHARACTER when not requested, an error occurs and the program stops!")
cat(" ",fill=TRUE)

# check if the dataset is an array or a dataframe or a matrix
check=0
if (is.matrix(dati)==FALSE){
	if (is.array(dati)==TRUE){
		n=dim(dati)[1]
		m=dim(dati)[2]
		p=dim(dati)[3]
		X=matrix(dati,n,m*p)
	}
	if(is.data.frame(dati)==TRUE){
		X=as.matrix(dati)
		check=1
	}
} else if (is.matrix(dati)==TRUE){
	X=dati
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
	if ((prep!=0) & (prep!=1) & (prep!=2) & (prep!=3)){
		cat(" ",fill=TRUE)
		cat("Error! Make a proper choice for normalizing the data",fill=TRUE)
		cat(" ",fill=TRUE)
	}
}

cat(" ",fill=TRUE)
cat("Note: The preprocessed data are now available in Xprep.",fill=TRUE)
cat("Xprep can be used for analyses outside the Tucker1 program",fill=TRUE)

cat(" ",fill=TRUE)
cat("In Tucker1 analysis, one of the three modes is reduced,",fill=TRUE)
cat("implying that three different types of model can be considered.",fill=TRUE)
check=0
while (check==0){
    cat("What model do you want to run?",fill=TRUE)
	cat("1 = type A (reduction of A-mode)", fill=TRUE)
	cat("2 = type B (reduction of B-mode)", fill=TRUE)
	cat("3 = type C (reduction of C-mode)", fill=TRUE)
	T1=scan("",n=1)
	if (length(T1)==0){	
		T1=0
	}
	if (T1==1){ 
		cat("A-mode will be reduced by Tucker1 (T1A)",fill=TRUE)
		check=1
	}
	if (T1==2){
		cat("B-mode will be reduced by Tucker1 (T1B)",fill=TRUE)
		check=1
	}
	if (T1==3){
		cat("C-mode will be reduced by Tucker1 (T1C)",fill=TRUE)
		check=1
	}
	if ((T1!=1) & (T1!=2) & (T1!=3)){
		cat(" ",fill=TRUE)
		cat("Error! Make a proper choice",fill=TRUE)
		cat(" ",fill=TRUE)
	}
}
	
cat(" ",fill=TRUE)
cat("You are now about to do a Tucker1 analysis.",fill=TRUE)

cat(" ",fill=TRUE)
cat("You have to specify the dimensionality to use.",fill=TRUE)
cat("To search this, you are advised to first run T1 analyses.",fill=TRUE)
cat("Note that T1 is equivalent to PCASUP.",fill=TRUE)
cat("(PCASUP analyses are PCA's of supermatrices with slices of the 3way array next to each other,",fill=TRUE)
cat("If you want to do this for choosing your dimensionality, specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
if (c==1){
	cat("It is needed to indicate the maximum number of dimension for the mode you want to study",fill=TRUE)
	if (T1==1){
		cat("Up to how many A-mode components do you want to use?",fill=TRUE)
		check0=0
		maxa=scan("",n=1)
		while (check0==0){
			if (length(maxa)==0){
				maxa==0
			}
			if ((maxa==0) | ((floor(maxa)-maxa)!=0) | (maxa>n)){
				cat(" ",fill=TRUE)
				cat("Error! Up to how many A-mode components do you want to use?",fill=TRUE)
				cat(" ",fill=TRUE)
			} else{
				check0=1
			}
		}	
	}	
	if (T1==2){
		cat("Up to how many B-mode components do you want to use?",fill=TRUE)
		check0=0
		while (check0==0){
			maxb=scan("",n=1)
			if (length(maxb)==0){
				maxb=0
			}
			if ((maxb==0) | ((floor(maxb)-maxb)!=0) | (maxb>m)){
				cat(" ",fill=TRUE)
				cat("Error! Up to how many B-mode components do you want to use?",fill=TRUE)
				cat(" ",fill=TRUE)
			} else{
				check0=1
			}
		}		
	}	
	if (T1==3){
		cat("Up to how many C-mode components do you want to use?",fill=TRUE)
		check0=0
		while (check0==0){
			maxc=scan("",n=1)
			if (length(maxc)==0){
				maxc=0
			}
			if ((maxc==0) | ((floor(maxc)-maxc)!=0) | (maxc>p)){
				cat(" ",fill=TRUE)
				cat("Error! Up to how many C-mode components do you want to use?",fill=TRUE)
				cat(" ",fill=TRUE)
			} else{
				check0=1
			}
		}
	}
	if (T1==1){
		maxb=m
		maxc=p
	}
	if (T1==2){
		maxa=n
		maxc=p
	}
	if (T1==3){
		maxa=n
		maxb=m
	}
	cat(" ",fill=TRUE)
	out=T1runsFit(Xprep,n,m,p,maxa,maxb,maxc,T1)
	cat(" ",fill=TRUE)
	cat("Note: T1 model written in terms of T3 model with numbers of components equal to numbers of entities for not reduced modes",fill=TRUE)
	cat("       Numbers of components   Fit   Total number of",fill=TRUE)
	rownames(out)=rep("T3",length=dim(out)[1])
	colnames(out)=c("      A","      B","      C","    (%) "," components")
	print(round(out,digits=2))
	if (nrow(out)>2){
		cat(" ",fill=TRUE)
		cat("Suggestion based on Convex Hull procedure",fill=TRUE)
		DS=DimSelector(out,n,m,p,4)
	}
	cat("Figure 1 now gives fit values of all solutions plotted against tnc=P+Q+R",fill=TRUE)
	cat("Figure 2 gives fit values (%) of all solutions plotted against number of free parameters",fill=TRUE)
	T3dimensionalityplot(out,n,m,p)
}

cat(" ",fill=TRUE)
cat("You can now do a TUCKER1 ANALYSIS",fill=TRUE)

cat(" ",fill=TRUE)
if (T1==1){
	cat("How many A-mode components do you want to use?",fill=TRUE)
	check=0
	while (check==0){
		r1=scan("",n=1)
		if (length(r1)==0){
			r1=0
		}
		if ((r1==0) | ((floor(r1)-r1)!=0) | (r1>n)){
			cat(" ",fill=TRUE)
			cat("Error! How many A-mode components do you want to use?",fill=TRUE)
			cat(" ",fill=TRUE)
		} else{
			check=1
		}
	}
}
if (T1==2){
	cat("How many B-mode components do you want to use?",fill=TRUE)
	check=0
	while (check==0){
		r2=scan("",n=1)
		if (length(r2)==0){
			r2=0
		}
		if ((r2==0) | ((floor(r2)-r2)!=0) | (r2>m)){
			cat(" ",fill=TRUE)
			cat("Error! How many B-mode components do you want to use?",fill=TRUE)
			cat(" ",fill=TRUE)
		} else{
			check=1
		}
	}
}
if (T1==3){
	cat("How many C-mode components do you want to use?",fill=TRUE)
	check=0
	while (check==0){
		r3=scan("",n=1)
		if (length(r3)==0){
			r3=0
		}
		if ((r3==0) | ((floor(r3)-r3)!=0) | (r3>p)){
			cat(" ",fill=TRUE)
			cat("Error! How many C-mode components do you want to use?",fill=TRUE)
			cat(" ",fill=TRUE)
		} else{
			check=1
		}
	}
}
if (T1==1){
	r2=m
	r3=p
	Tuck1=pcasup1(Xprep,n,m,p,T1)
	A=as.matrix(Tuck1$A[,1:r1])
	B=diag(1,nrow=r2)
	C=diag(1,nrow=r3)
	fp=sum(Tuck1$la[1:r1])/sum(Tuck1$la)*100
}
if (T1==2){
	r1=n
	r3=p
	A=diag(1,nrow=r1)
	Tuck1=pcasup1(Xprep,n,m,p,T1)
	B=as.matrix(Tuck1$B[,1:r2])
	C=diag(1,nrow=r3)
	fp=sum(Tuck1$lb[1:r2])/sum(Tuck1$lb)*100
}
if (T1==3){
	r1=n
	r2=m
	A=diag(1,nrow=r1)
	B=diag(1,nrow=r2)
	Tuck1=pcasup1(Xprep,n,m,p,T1)
	C=as.matrix(Tuck1$C[,1:r3])
	fp=sum(Tuck1$lc[1:r3])/sum(Tuck1$lc)*100
}

H=permnew(t(A)%*%Xprep,r1,m,p)
H=permnew(t(B)%*%H,r2,p,r1)
H=permnew(t(C)%*%H,r3,r1,r2)
H=as.matrix(H)
HH=H^2/sum(Xprep^2)*100

names(fp)="Fit (%)"

cat(" ",fill=TRUE)
if (T1==1){
	cat(paste("Tucker1A analysis with ",r1," components, gave a fit of ",round(fp,digits=2), "%"),fill=TRUE) 
}
if (T1==2){
	cat(paste("Tucker1B analysis with ",r2," components, gave a fit of ",round(fp,digits=2), "%"),fill=TRUE) 
}
if (T1==3){
	cat(paste("Tucker1C analysis with ",r3," components, gave a fit of ",round(fp,digits=2), "%"),fill=TRUE) 
}

# make labels for columns of core
str=noquote(vector(mode="character",length=r3*r2))
i=1
for (k in 1:r3){
	for (j in 1:r2){
		str[i]=noquote(paste("  B",as.character(j),"xC",as.character(k),sep=""))
		i=i+1
	}
}
corelabelstr=str
labCompA=paste("A",1:r1,sep="")
labCompB=paste("B",1:r2,sep="")
labCompC=paste("C",1:r3,sep="")
rownames(A)=laba
rownames(B)=labb
rownames(C)=labc
rownames(H)=labCompA
colnames(A)=labCompA
colnames(B)=labCompB
colnames(C)=labCompC
colnames(H)=corelabelstr

cat(" ",fill=TRUE)
cat("Solution for A, B and C, core in H in detail.",fill=TRUE)
cat("Only the matrix for reduced mode and the core will be given.",fill=TRUE)
if (T1==1){
	cat("A",fill=TRUE)
	print(round(A,digits=2))
}
if (T1==2){
	cat("B",fill=TRUE)
	print(round(B,digits=2))
}
if (T1==3){
	cat("C",fill=TRUE)
	print(round(C,digits=2))
}
cat("Core",fill=TRUE)
print(round(H,digits=2))

cat(" ",fill=TRUE)
cat("Find useful simple structure rotations of core and components",fill=TRUE) 
# varimax
cat("You can now carry out SIMPLE STRUCTURE rotations by varimax.",fill=TRUE)
cat("If you want to varimax rotate the solution, specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
if (c==1){
	vvc=1
	while(vvc==1){
		if (T1==1){
			if (r1>1){
				check=1
				while (check==1){
					cat("Which matrix do you want to rotate?",fill=TRUE)
					cat("1 = matrix A",fill=TRUE)
					cat("2 = core (row-wise varimax)",fill=TRUE)
					vc=scan("",n=1)
					if (length(vc)==0){
						vc=0
					}
					if ((vc!=1) & (vc!=2)){
						cat(" ",fill=TRUE)
						cat("Error! Make a proper choice",fill=TRUE)
						cat(" ",fill=TRUE)
					} else{
						check=0
					}
				}	
				if (vc==1){
					VAR=varim(A)
					AS=VAR$B
					BT=B
					CU=C
					K=t(VAR$T)%*%H
				} else{
					VAR=varim(t(H))
					AS=A%*%VAR$T
					BT=B
					CU=C
					K=t(VAR$B)
				}
			} else{
				cat(" ",fill=TRUE)
				cat("Warning: as the number of A-mode components is 1, no simple structure will be sought",fill=TRUE)
				cat(" ",fill=TRUE)
				AS=A
				BT=B
				CU=C
				K=H
			}
		}	
		if (T1==2){
			if (r2>1){
				H=permnew(H,r1,r2,r3)
				check=1
				while (check==1){
					cat("Which matrix do you want to rotate?",fill=TRUE)
					cat("1 = matrix B",fill=TRUE)
					cat("2 = core (row-wise varimax)",fill=TRUE)
					vc=scan("",n=1)
					if (length(vc)==0){
						vc=0
					}
					if ((vc!=1) & (vc!=2)){
						cat(" ",fill=TRUE)
						cat("Error! Make a proper choice",fill=TRUE)
						cat(" ",fill=TRUE)
					} else{
						check=0
					}
				}	
				if (vc==1){
					VAR=varim(B)
					AS=A
					BT=VAR$B
					CU=C
					K=t(VAR$T)%*%H
				} else{
					VAR=varim(t(H))
					AS=A
					BT=B%*%VAR$T
					CU=C
					K=t(VAR$B)
				}
				K=permnew(K,r2,r3,r1)
				K=permnew(K,r3,r1,r2)
			} else{
				cat(" ",fill=TRUE)
				cat("Warning: as the number of B-mode components is 1, no simple structure will be sought",fill=TRUE)
				cat(" ",fill=TRUE)
				AS=A
				BT=B
				CU=C
				K=H
			}
		}	
		if (T1==3){
			if (r3>1){
				H=permnew(H,r1,r2,r3)
				H=permnew(H,r2,r3,r1)
				check=1
				while (check==1){
					cat("Which matrix do you want to rotate?",fill=TRUE)
					cat("1 = matrix C",fill=TRUE)
					cat("2 = core (row-wise varimax)",fill=TRUE)
					vc=scan("",n=1)
					if (length(vc)==0){
						vc=0
					}
					if ((vc!=1) & (vc!=2)){
						cat(" ",fill=TRUE)
						cat("Error! Make a proper choice",fill=TRUE)
						cat(" ",fill=TRUE)
					} else{
						check=0
					}
				}	
				if (vc==1){
					VAR=varim(C)
					AS=A
					BT=B
					CU=VAR$B
					K=t(VAR$T)%*%H
				} else{
					VAR=varim(t(H))
					AS=A
					BT=B
					CU=C%*%VAR$T
					K=t(VAR$B)
				}
				K=permnew(K,r3,r1,r2)
			} else{
				cat(" ",fill=TRUE)
				cat("Warning: as the number of C-mode components is 1, no simple structure will be sought",fill=TRUE)
				cat(" ",fill=TRUE)
				AS=A
				BT=B
				CU=C
				K=H
			}
		}
		rownames(AS)=laba
		rownames(BT)=labb
		rownames(CU)=labc
		rownames(K)=labCompA
		colnames(AS)=labCompA
		colnames(BT)=labCompB
		colnames(CU)=labCompC
		colnames(K)=corelabelstr

		cat(" ",fill=TRUE)
		cat("Rotated solution for A, B and C is in AS, BT and CU, rotated core in K.",fill=TRUE)
		cat("Only the matrix for reduced mode and the core will be given.",fill=TRUE)
		cat("Note: Solution is in AS, BT, CU and K even if no simple structure rotation was carried out.",fill=TRUE)
		if (T1==1){
			cat("Rotated A (AS)",fill=TRUE)
			print(round(A,digits=2))
		}
		if (T1==2){
			cat("Rotated B (BT)",fill=TRUE)
			print(round(B,digits=2))
		}
		if (T1==3){
			cat("Rotated C (CU)",fill=TRUE)
			print(round(CU,digits=2))
		}
		cat("Rotated core",fill=TRUE)
		print(round(K,digits=2))
		cat(" ",fill=TRUE)
		cat("The last varimax rotated solution will be used.",fill=TRUE);
		cat("If you want to study one more solution, specify '1':",fill=TRUE)
		vvc=scan("",n=1)
		if (length(vvc)==0){
			vvc=0
		}
	}
} else{
	AS=A
	BT=B
	CU=C
	K=H
	rownames(AS)=laba
	rownames(BT)=labb
	rownames(CU)=labc
	rownames(K)=labCompA
	colnames(AS)=labCompA
	colnames(BT)=labCompB
	colnames(CU)=labCompC
	colnames(K)=corelabelstr
}

# manual permutation and reflection:
cat(" ",fill=TRUE)
cat("You can now manually PERMUTE and REFLECT columns/rows of solution",fill=TRUE)
cat("If you want to reflect/permute columns/rows, specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
while (c==1){
	if (T1==1){
		cat("Give a vector for reflection of columns of A (e.g., 1 -1 -1 1 ..)",fill=TRUE)
		tau=scan("",n=r1)
		if (length(tau)==0){
			tau=array(1,r1)
			cat(" ",fill=TRUE)
			cat("Warning: the columns of A will not be reflected", fill=TRUE)
			cat(" ",fill=TRUE)
		} else{
			if ((sum(tau^2)!=r1) | (sum((tau^2-1)^2)!=0)){
				cat(" ",fill=TRUE)
				cat("Error! Incorrect reflection! The columns of A will not be reflected!", fill=TRUE)
				cat(" ",fill=TRUE)
				tau=array(1,r1)
			}
		}
		if (r1>1){
			cat("Give a vector with new order of columns of A (e.g., 3 1 4 2 ..)",fill=TRUE)
			pp=scan("",n=r1)
			if (length(pp)==0){
				pp=c(1:r1)
				cat(" ",fill=TRUE)
				cat("Warning: the columns of A will not be permuted", fill=TRUE)
				cat(" ",fill=TRUE)
			} else{
				if ((length(pp)!=r1) | (sum((sort(pp)-1:r1)^2)>0)){	
					cat(" ",fill=TRUE)
					cat("Error! Incorrect permutation! The columns of A will not be permuted", fill=TRUE)
					cat(" ",fill=TRUE)
					pp=c(1:r1)
				}
			}
		} else{
			pp=1
		}
		AS=AS%*%diag(tau,nrow=r1)		
		K=diag(tau,nrow=r1)%*%K			
		AS=as.matrix(AS[,pp])
		K=as.matrix(K[pp,])
		if (r1==1){
			K=t(K)
		}
		rownames(K)=labCompA
		colnames(AS)=labCompA
		colnames(K)=corelabelstr
		cat("Rotated A (AS)",fill=TRUE)
		print(round(AS,digits=2))
		cat("Rotated K (K)",fill=TRUE)
		print(round(K,digits=2))
	}
	if (T1==2){	
		cat("Give a vector for reflection of columns of B (e.g., 1 -1 -1 1 ..)",fill=TRUE)
		tau=scan("",n=r2)
		if (length(tau)==0){
			tau=array(1,r2)
			cat(" ",fill=TRUE)
			cat("Warning: the columns of B will not be reflected", fill=TRUE)
			cat(" ",fill=TRUE)
		}	
		else{
			if ((sum(tau^2)!=r2) | (sum((tau^2-1)^2)!=0)){
				cat(" ",fill=TRUE)
				cat("Error! Incorrect reflection! The columns of B will not be reflected", fill=TRUE)
				cat(" ",fill=TRUE)
				tau=array(1,r2)
			}
		}	
		if (r2>1){
			cat("Give vector with new order of columns of B (e.g., 3 1 4 2 ..)",fill=TRUE)
			pp=scan("",n=r2)
			if (length(pp)==0){
				pp=c(1:r2)
				cat(" ",fill=TRUE)
				cat("Warning: the columns of B will not be permuted", fill=TRUE)
				cat(" ",fill=TRUE)
			} else{
				if ((length(pp)!=r2) | (sum((sort(pp)-1:r2)^2)>0)){	
					cat(" ",fill=TRUE)
					cat("Error! Incorrect permutation! The columns of B will not be permuted", fill=TRUE)
					cat(" ",fill=TRUE)
					pp=c(1:r2)
				}
			}	
		} else{
			pp=1
		}
		rownames(K)=labCompA
		colnames(BT)=labCompB
		colnames(K)=corelabelstr
		K=permnew(K,r1,r2,r3)
		BT=BT%*%diag(tau,nrow=r2)				
		K=diag(tau,nrow=r2)%*%K
		BT=as.matrix(BT[,pp])
		K=as.matrix(K[pp,])
		if (r2==1){
			K=t(K)
		}
		K=permnew(K,r2,r3,r1)
		K=permnew(K,r3,r1,r2)
		cat("Rotated B (BT)",fill=TRUE)
		print(round(BT,digits=2))
		cat("Rotated core",fill=TRUE)
		print(round(K,digits=2))
	}
	if (T1==3){	
		cat("Give a vector for reflection of columns of C (e.g., 1 -1 -1 1 ..)",fill=TRUE)
		tau=scan("",n=r3)
		if (length(tau)==0){
			tau=array(1,r3)
			cat(" ",fill=TRUE)
			cat("Warning: the columns of C will not be reflected", fill=TRUE)
			cat(" ",fill=TRUE)
		} else{
			if ((sum(tau^2)!=r3) | (sum((tau^2-1)^2)!=0)){
				cat(" ",fill=TRUE)
				cat("Error! Incorrect reflection! The columns of C will not be reflected", fill=TRUE)
				cat(" ",fill=TRUE)
				tau=array(1,r3)
			}
		}
		if (r3>1){
			cat("Give a vector with new order of columns of C (e.g., 3 1 4 2 ..)",fill=TRUE)
			pp=scan("",n=r3)
			if (length(pp)==0){
				pp=c(1:r3)
				cat(" ",fill=TRUE)
				cat("Warning: the columns of C will not be permuted", fill=TRUE)
				cat(" ",fill=TRUE)
			} else{
				if ((length(pp)!=r3) | (sum((sort(pp)-1:r3)^2)>0)){	
					cat(" ",fill=TRUE)
					cat("Error! Incorrect permutation! The columns of C will not be permuted!", fill=TRUE)
					cat(" ",fill=TRUE)
					pp=c(1:r3)
				}
			}	
		} else{
			pp=1
		}
		K=permnew(K,r1,r2,r3)
		K=permnew(K,r2,r3,r1)
		CU=CU%*%diag(tau,nrow=r3)		
		K=diag(tau,nrow=r3)%*%K
		CU=as.matrix(CU[,pp])
		K=as.matrix(K[pp,])
		K=permnew(K,r3,r1,r2)
		if (r3==1){
			K=t(K)
		}
		rownames(K)=labCompA
		colnames(CU)=labCompC
		colnames(K)=corelabelstr
		cat("Rotated C (CU)",fill=TRUE)
		print(round(CU,digits=2))
		cat("Rotated core",fill=TRUE)
		print(round(K,digits=2))
	}
	cat(" ",fill=TRUE)
	cat("If you want to further reflect/permute columns/rows, specify '1'",fill=TRUE)
	c=scan("",n=1)
	if (length(c)==0){
		c=0
	}
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
	#colnames(K)=corelabelstr
	#print(round((K^2)*100/sum(Xprep^2),digits=2),rownames=labCompA)
	print(round((K^2)*100/sum(Xprep^2),digits=2))
	renormmode=0
	T1fit=T3fitpartitioning(Xprep,n,m,p,AS,BT,CU,K,renormmode,laba,labb,labc)		
	cat("Relative fit of A-mode entities (in %)", fill=TRUE)
	print(round(T1fit$fitA,digits=2))
	cat("Relative fit of B-mode entities (in %)",fill=TRUE)
	print(round(T1fit$fitB,digits=2))
	cat("Relative fit of C-mode entities (in %)",fill=TRUE)
	print(round(T1fit$fitC,digits=2))
	cat("Results for the fitpartitioning are given in $fitA, $fitB, $fitC",fill=TRUE)
} else{
	T1fit=list()
	T1fit$fitA=NULL
	T1fit$fitB=NULL
	T1fit$fitC=NULL
}

cat(" ",fill=TRUE)
cat("Press 'return' to conclude the analysis",fill=TRUE)
c=scan("",n=1)
cat("To see (rotated) component matrices, core matrix, fit and other results: $Xprep, $A, $B, $C, $core, $fit, $fitA, $fitB, $fitC, ...", fill=TRUE)

XprepOut=rarray(Xprep,n,m,p)
dimnames(XprepOut)=list(laba,labb,labc)

out=list()
out$A=AS
out$B=BT
out$C=CU
out$core=K
out$fit=fp
out$fitA=T1fit$fitA
out$fitB=T1fit$fitB
out$fitC=T1fit$fitC
out$laba=laba
out$labb=labb
out$labc=labc
out$Xprep=XprepOut
return(out)
}

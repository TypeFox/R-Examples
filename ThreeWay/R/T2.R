T2 <-
function(dati,laba,labb,labc){

cat(" ",fill=TRUE)
cat("WELCOME to the interactive TUCKER2 analysis program",fill=TRUE)
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
cat("Xprep can be used for analyses outside the Tucker2 program",fill=TRUE)

cat(" ",fill=TRUE)
cat("In Tucker2 analysis, two of the three modes are reduced,",fill=TRUE)
cat("implying that three different types of model can be considered.",fill=TRUE)
check=0
while (check==0){
    cat("What model do you want to run?",fill=TRUE)
	cat("1 = type AB (reduction of A-mode and B-mode)", fill=TRUE)
	cat("2 = type AC (reduction of A-mode and C-mode)", fill=TRUE)
	cat("3 = type BC (reduction of B-mode and C-mode)", fill=TRUE)
	T2=scan("",n=1)
	if (length(T2)==0){	
		T2=0
	}
	if (T2==1){ 
		cat("A-mode and B-mode will be reduced by Tucker2 (T2AB)",fill=TRUE)
		check=1
	}
	if (T2==2){
		cat("A-mode and C-mode will be reduced by Tucker2 (T2AC)",fill=TRUE)
		check=1
	}
	if (T2==3){
		cat("B-mode and C-mode will be reduced by Tucker2 (T2BC)",fill=TRUE)
		check=1
	}
	if ((T2!=1) & (T2!=2) & (T2!=3)){
		cat(" ",fill=TRUE)
		cat("Error! Make a proper choice",fill=TRUE)
		cat(" ",fill=TRUE)
	}
}
	
cat(" ",fill=TRUE)
cat("You are now about to do a Tucker2 analysis.",fill=TRUE)

cat(" ",fill=TRUE)
cat("You have to specify the dimensionalities to use.",fill=TRUE)
cat("To search these, you are advised to first run PCASUP analyses.",fill=TRUE)
cat("(PCASUP analyses are PCA's of supermatrices with slices of the 3way array next to each other,",fill=TRUE)
cat("thus 2 supermatrices are analyed by PCA (considering the two reduced modes)",fill=TRUE) 
cat("and for each we get a component matrix, and eigenvalues.",fill=TRUE)
cat("The eigenvalues may give an indication as to the required number of components for each mode.)",fill=TRUE)
cat("The results can next be used to find useful dimensionalities by the generalized scree test",fill=TRUE)
cat("(see Timmerman & Kiers, 2000, Kiers & der Kinderen, 2003)",fill=TRUE)
cat("If you want to do the PCASUPs for choosing your dimensionality , specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
if (c==1){
	cat("For the generalized scree test it is needed to indicate the maximum number of dimensions for each mode you want to study",fill=TRUE)
	if (T2!=3){
		cat("Up to how many A-mode components do you want to use?",fill=TRUE)
		check0=0
		while (check0==0){
			maxa=scan("",n=1)
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
	if (T2!=2){
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
		if (T2!=1){
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
	if (T2==1){
		maxc=p
	}
	if (T2==2){
		maxb=m
	}
	if (T2==3){
		maxa=n
	}
	cat(" ",fill=TRUE)
	out=T2runsApproxFit(Xprep,n,m,p,maxa,maxb,maxc,T2)
	cat(" ",fill=TRUE)
	cat("Note: T2 model written in terms of T3 model with number of components equal to number of entities for not reduced mode",fill=TRUE)
	cat("       Numbers of components   Fit   Total number of",fill=TRUE)
	rownames(out)=rep("T2",length=dim(out)[1])
	colnames(out)=c("      A","      B","      C","    (%) "," components")
	print(round(out,digits=2))
	if (nrow(out)>2){
		cat(" ",fill=TRUE)
		cat("Suggestion based on Convex Hull procedure",fill=TRUE)
		DS=DimSelector(out,n,m,p,3)
	}
	cat("Figure 1 now gives fit values of all solutions plotted against tnc=P+Q+R",fill=TRUE)
	cat("Figure 2 gives fit values (%) of all solutions plotted against number of free parameters",fill=TRUE)
	T3dimensionalityplot(out,n,m,p)
}

cat(" ",fill=TRUE)
cat("You can now do a TUCKER2 ANALYSIS",fill=TRUE)

cat(" ",fill=TRUE)
if (T2!=3){
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
if (T2!=2){
	cat("How many B-mode components do you want to use?",fill=TRUE)
	check=0
	while (check==0){r2=scan("",n=1)
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
if (T2!=1){
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
if (T2==1){
	r3=p
}
if (T2==2){
	r2=m
}
if (T2==3){
	r1=n
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

cat(paste("Run no.",1,sep=" "),fill=TRUE)
Tuck2=T2func(Xprep,n,m,p,r1,r2,r3,0,conv,T2)		
A=Tuck2$A
B=Tuck2$B
C=Tuck2$C
H=Tuck2$H
La=Tuck2$La
Lb=Tuck2$Lb
Lc=Tuck2$Lc
iter=Tuck2$iter
cput=Tuck2$cputime
fp=Tuck2$fp
f=Tuck2$f
funcF=vector("numeric",length=1+addanal)
funcV=vector("numeric",length=1+addanal)
cputime=vector("numeric",length=1+addanal)
niter=vector("numeric",length=1+addanal)
namev=paste("Start n.",1:(1+addanal),sep="")
names(funcF)=namev
names(funcV)=namev
names(cputime)=namev
names(niter)=namev
funcF[1]=fp
funcV[1]=f
cputime[1]=cput
niter[1]=iter
if (addanal>=1){
	for (run in 1:addanal){
		cat(paste("Run no.",run+1,sep=" "),fill=TRUE)
		start=1
		Tuck2b=T2func(Xprep,n,m,p,r1,r2,r3,start,conv,T2)			
		funcF[run+1]=Tuck2b$fp
		funcV[run+1]=Tuck2b$f
		cputime[run+1]=Tuck2b$cputime
		niter[run+1]=Tuck2b$iter
		if (Tuck2b$fp>1.0001*Tuck2$fp){		# if fit more than .01% better is found, replace solution
			A=Tuck2b$A
			B=Tuck2b$B
			C=Tuck2b$C
			H=Tuck2b$H
			fp=Tuck2b$fp
			La=Tuck2b$La
			Lb=Tuck2b$Lb
			Lc=Tuck2b$Lc
			fp=Tuck2b$fp
			f=Tuck2b$f
		}
	}
	cat(" ",fill=TRUE)
	cat("Fit (%) values from all runs:",fill=TRUE)
	print(round(funcF,digits=2))
}
names(fp)="Fit (%)"
cat(" ",fill=TRUE)
if (T2==1){
	cat(paste("Tucker2AB analysis with ",r1,"x",r2," components, gave a fit of ",round(fp,digits=2), "%"),fill=TRUE) 
}
if (T2==2){
	cat(paste("Tucker2AC analysis with ",r1,"x",r3," components, gave a fit of ",round(fp,digits=2), "%"),fill=TRUE) 
}
if (T2==3){
	cat(paste("Tucker2BC analysis with ",r2,"x",r3," components, gave a fit of ",round(fp,digits=2), "%"),fill=TRUE) 
}

# Prepare plotting matrices
cat(" ",fill=TRUE)
cat("PLOTTING",fill=TRUE)
cat("No special plotting routines are incorporated, but you can use any scatter-plotting routine you like",fill=TRUE)
cat("If you want to make plots of entities of one mode, you can use $Aplot, $Bplot, or $Cplot, and use labels $laba, $labb, $labc.",fill=TRUE)
Aplot=A%*%diag(diag(La)^0.5,nrow=r1)
Bplot=B%*%diag(diag(Lb)^0.5,nrow=r2)
Cplot=C%*%diag(diag(Lc)^0.5,nrow=r3)
cat("If you want to make plots of entities of two modes, with the third projected in it as axes, you can use as coordinates: ",fill=TRUE)
cat("($BAplot,$C),($ACplot,$B), or ($CBplot,$A)",fill=TRUE)
Ha=H
Hb=permnew(Ha,r1,r2,r3)
Hc=permnew(Hb,r2,r3,r1)
CBplot=kronecker(C,B)%*%t(Ha)
ACplot=kronecker(A,C)%*%t(Hb)
BAplot=kronecker(B,A)%*%t(Hc)

out=list()

cat(" ",fill=TRUE)
cat("It is sometimes useful to RENORMALIZE components, e.g., such that",fill=TRUE)
cat(" the influence of the different A-mode components on the core is equalized",fill=TRUE)
cat("This will usually lead to correlated components in that mode, and",fill=TRUE)
cat(" to a more simple structure for the components in the OTHER modes, ",fill=TRUE)
cat("If you want to renormalize components, specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
renormmode=0
if (c==1){
	check1=0
	while (check1==0){
		cat(" ",fill=TRUE)
		if (T2==1){
			cat("Which mode do you want to renormalize? Enter 1 (=A) or 2 (=B):",fill=TRUE)
			renormmode=scan("",n=1)
			if (length(renormmode)==0){
				renormmode=0
			}
			if ((renormmode==1) | (renormmode==2)){
				check1=1
				RNsol=renormsolT3(A,B,C,H,renormmode)
			} else{
				cat(" ",fill=TRUE)
				cat("Error! Select a proper number!")
				cat(" ",fill=TRUE)
			}
		}		
		if (T2==2){
			cat("Which mode do you want to renormalize? Enter 1 (=A) or 3 (=C):",fill=TRUE)
			renormmode=scan("",n=1)
			if (length(renormmode)==0){
				renormmode=0
			}
			if ((renormmode==1) | (renormmode==3)){
				check1=1
				RNsol=renormsolT3(A,B,C,H,renormmode)
			} else{
				cat(" ",fill=TRUE)
				cat("Error! Select a proper number!")
				cat(" ",fill=TRUE)
			}
		}				
		if (T2==3){
			cat("Which mode do you want to renormalize? Enter 2 (=B) or 3 (=C):",fill=TRUE)
			renormmode=scan("",n=1)
			if (length(renormmode)==0){
				renormmode=0
			}
			if ((renormmode==2) | (renormmode==3)){
				check1=1
				RNsol=renormsolT3(A,B,C,H,renormmode)
			} else{
				cat(" ",fill=TRUE)
				cat("Error! Select a proper number!")
				cat(" ",fill=TRUE)
			}
		}
	}
	A=RNsol$A
	B=RNsol$B
	C=RNsol$C
	H=RNsol$H
}

cat(" ",fill=TRUE)
cat("Find useful simple structure rotation of core and components",fill=TRUE) 
c=1
while (c==1){
	# varimcocoruns
	# program to run varimcoco repeatedly with different relative weights for A,B,C
	cat("You can now carry out SIMPLE STRUCTURE rotations with varying weights.",fill=TRUE)
	cat("If desired, you can specify a range of different weights for each mode.",fill=TRUE)
#	cat("To specify a range of values you type, e.g., 1:4 (then the values 1,2,3 and 4 are chosen),",fill=TRUE)
#	cat(" or seq(5,25,5) (then the values 5,10,15,20 and 25 are chosen),",fill=TRUE)
#	cat(" or c(1,2,5,10,100) (then the values 1,2,5,10 and 100 are chosen).",fill=TRUE)
	cat("You can also specify a single weight (e.g., just 1).",fill=TRUE)
	cat("Analyses will be carried out with all combinations of relative weights.",fill=TRUE)
	if (T2!=3){
		if (r1>1){
			cat("Specify (range of) relative weight(s) for A (default=0):",fill=TRUE)
			aa=scan()
			if (length(aa)==0){
				aa=0
			}
		} else{	
			cat(" ",fill=TRUE)
			cat("Warning: as the number of A-mode components is 1, no simple structure for A will be sought (relative weight=0)",fill=TRUE)
			cat(" ",fill=TRUE)
			aa=0
		}	
	} else{
		cat(" ",fill=TRUE)
		cat("Warning: as the A-mode was not reduced, no simple structure for A will be sought (relative weight=0)",fill=TRUE)
		cat(" ",fill=TRUE)
		aa=0
	}
	if (T2!=2){
		if (r2>1){
			cat("Specify (range of) relative weight(s) for B (default=0):",fill=TRUE)
			bb=scan()
			if (length(bb)==0){
				bb=0
			}
		} else{	
			cat(" ",fill=TRUE)
			cat("Warning: as the number of B-mode components is 1, no simple structure for B will be sought (relative weight=0)",fill=TRUE)
			cat(" ",fill=TRUE)
			bb=0
		}	
	} else{
		cat(" ",fill=TRUE)
		cat("Warning: as the B-mode was not reduced, no simple structure for B will be sought (relative weight=0)",fill=TRUE)
		cat(" ",fill=TRUE)
		bb=0
	}
	if (T2!=1){
		if (r3>1){
			cat("Specify (range of) relative weight(s) for C (default=0):",fill=TRUE)
			cc=scan()
			if (length(cc)==0){
				cc=0
			}
		} else{	
			cat(" ",fill=TRUE)
			cat("Warning: as the number of C-mode components is 1, no simple structure for C will be sought (relative weight=0)",fill=TRUE)
			cat(" ",fill=TRUE)
			cc=0
		}
	} else{
		cat(" ",fill=TRUE)
		cat("Warning: as the C-mode was not reduced, no simple structure for C will be sought (relative weight=0)",fill=TRUE)
		cat(" ",fill=TRUE)
		cc=0
	}
	for (wa_rel in aa){
		for (wb_rel in bb){
			for (wc_rel in cc){
				cat(" ",fill=TRUE)
				cat(paste("Simple structure rotation with relative weights ", wa_rel, wb_rel, wc_rel),fill=TRUE)
				noplus=0
				if (T2==1){
					VAR=varimcoco(A,B,C,H,wa_rel,wb_rel,wc_rel,1,1,0)
					if (max(r1,r2)==1){
						noplus=1
						c=0
					}
				}
				if (T2==2){
					VAR=varimcoco(A,B,C,H,wa_rel,wb_rel,wc_rel,1,0,1)	
					if (max(r1,r3)==1){
						noplus=1
						c=0
					}
				}
				if (T2==3){
					VAR=varimcoco(A,B,C,H,wa_rel,wb_rel,wc_rel,0,1,1)	
					if (max(r2,r3)==1){
						noplus=1
						c=0
					}
				}
				AS=VAR$AS
				BT=VAR$BT
				CU=VAR$CU
				K=VAR$K
				out=c(wa_rel, wb_rel, wc_rel, VAR$f1, VAR$f2a, VAR$f2b, VAR$f2c)
			}
		}
	}
	cat("         relative weights          simplicity values",fill=TRUE)
	cat("        A       B       C       core     A       B       C",fill=TRUE)
    cat(paste("       ",round(out[1],digits=2),"     ",round(out[2],digits=2),"     ",round(out[3],digits=2),"     ",round(out[4],digits=2),"   ",round(out[5],digits=2),"   ",round(out[6],digits=2),"  ",round(out[7],digits=2)),fill=TRUE)
	cat(" ",fill=TRUE)
	if (noplus==0){
		cat("If you want another series of simple structure analyses, specify '1':",fill=TRUE)
		c=scan("",n=1)
		if (length(c)==0){
			c=0
		}
	}
}

cat(" ",fill=TRUE)
cat("SIMPLE STRUCTURE rotated solution in detail",fill=TRUE)
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

c=1
while (c==1){
	cat("The last simple structure rotated solution will be used.",fill=TRUE);
	if (noplus==0){
		cat("To specify a different simple structure rotation, specify '1':",fill=TRUE)
		c=scan("",n=1)
		if (length(c)==0){
			c=0
		}
		if (c==1){
			cat("Specify your choice for the preferred simple structure rotation of core and components",fill=TRUE)
			if (T2!=3){
				if (r1>1){
					cat("Specify relative weight for simplicity of A-mode (default=0):",fill=TRUE)
					wa_rel=scan("",n=1)
					if (length(wa_rel)==0){
						wa_rel=0
					}
				} else{	
					wa_rel=0
				}	
			}
			if (T2!=2){
				if (r2>1){
					cat("Specify relative weight for simplicity of B-mode (default=0):",fill=TRUE)
					wb_rel=scan("",n=1)
					if (length(wb_rel)==0){
						wb_rel=0
					}
				} else{	
					wb_rel=0
				}	
			}
			if (T2!=1){
				if (r3>1){
					cat("Specify relative weight for simplicity of C-mode (default=0):",fill=TRUE)
					wc_rel=scan("",n=1)
					if (length(wc_rel)==0){
						wc_rel=0
					}
				} else{	
					wc_rel=0
				}	
			}
			if (T2==1){
				VAR=varimcoco(A,B,C,H,wa_rel,wb_rel,wc_rel,1,1,0)	
			}
			if (T2==2){
				VAR=varimcoco(A,B,C,H,wa_rel,wb_rel,wc_rel,1,0,1)	
			}
			if (T2==3){
				VAR=varimcoco(A,B,C,H,wa_rel,wb_rel,wc_rel,0,1,1)	
			}
			AS=VAR$AS
			BT=VAR$BT
			CU=VAR$CU
			K=VAR$K
		}	
	}	
	labCompA=paste("A",1:r1,sep="")
	labCompB=paste("B",1:r2,sep="")
	labCompC=paste("C",1:r3,sep="")
	rownames(AS)=laba
	rownames(BT)=labb
	rownames(CU)=labc
	rownames(K)=labCompA
	colnames(AS)=labCompA
	colnames(BT)=labCompB
	colnames(CU)=labCompC
	colnames(K)=corelabelstr
	cat("Backnormalize A,B,C, and H",fill=TRUE)
	if (renormmode==1){
		Ds=diag(SUM(AS)$col^.5,nrow=r1)
		AS=AS%*%solve(Ds)
		K=Ds%*%K
		rownames(K)=labCompA
		colnames(AS)=labCompA
		colnames(K)=corelabelstr
	}
	if (renormmode==2){
		Ds=diag(SUM(BT)$col^.5,nrow=r2)
		BT=BT%*%solve(Ds)
		K=K%*%kronecker(diag(r3),Ds)
		rownames(K)=labCompA
		colnames(BT)=labCompB
		colnames(K)=corelabelstr
	}
	if (renormmode==3){
		Ds=diag(SUM(CU)$col^.5,nrow=r3)
		CU=CU%*%solve(Ds)
		K=K%*%kronecker(Ds,diag(r2))
		rownames(K)=labCompA
		colnames(CU)=labCompC
		colnames(K)=corelabelstr
	}
	cat(" ",fill=TRUE)
	cat("Rotated solution for A, B and C is in AS, BT and CU, rotated core in K.",fill=TRUE)
	cat("Only the matrices for reduced modes and the core will be given.",fill=TRUE)
	cat("Note: Solution is in AS, BT, CU and K even if no simple structure rotation was carried out.",fill=TRUE)
	if (T2!=3){
		cat("Rotated A (AS)",fill=TRUE)
		print(round(AS,digits=2))
	}
	if (T2!=2){
		cat("Rotated B (BT)",fill=TRUE)
		print(round(BT,digits=2))
	}
	if (T2!=1){
		cat("Rotated C (CU)",fill=TRUE)
		print(round(CU,digits=2))
	}
	cat("Rotated core",fill=TRUE)
	print(round(K,digits=2))
	cat(" ",fill=TRUE)

	if (noplus==0){
		cat("You can keep present as final solution, or study a different one.",fill=TRUE)
		cat("If you want to study one more solution, specify '1':",fill=TRUE)
		c=scan("",n=1)
		if (length(c)==0){
			c=0
		}
	} else{
		c=0
	}
}

# manual permutation and reflection:
cat(" ",fill=TRUE)
cat("You can now manually PERMUTE and REFLECT columns/rows of solution",fill=TRUE)
cat("If you want to reflect/permute columns/rows, specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
disp=0
while (c==1){
	if (disp>0){
		cat("Rotated solution for A, B and C is in AS, BT and CU, rotated core in K",fill=TRUE)
		cat("Only the matrices for reduced modes and the core will be given.",fill=TRUE)
		cat("Note: Solution will be in AS, BT, CU and K even if no simple structure rotation was carried out.",fill=TRUE)
		if (T2!=3){
			cat("Rotated A (AS)",fill=TRUE)
			print(round(AS,digits=2))
		}
		if (T2!=2){
			cat("Rotated B (BT)",fill=TRUE)
			print(round(BT,digits=2))
		}
		if (T2!=1){
			cat("Rotated C (CU)",fill=TRUE)
			print(round(CU,digits=2))
		}
		cat("Rotated core",fill=TRUE)
		print(round(K,digits=2))
	}
	disp=1
	if (T2!=3){
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
					cat("Error! Incorrect permutation! The columns of A will not be permuted!", fill=TRUE)
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
		if (T2==1){
			cat(" ",fill=TRUE)
			cat("Rotated B (BT)",fill=TRUE)
			print(round(BT,digits=2))
		} else{
			cat(" ",fill=TRUE)
			cat("Rotated C (CU)",fill=TRUE)
			print(round(CU,digits=2))
		}
		cat("Rotated K (K)",fill=TRUE)
		print(round(K,digits=2))
	}
	if (T2!=2){	
		cat("Give a vector for reflection of columns of B (e.g., 1 -1 -1 1 ..)",fill=TRUE)
		tau=scan("",n=r2)
		if (length(tau)==0){
			tau=array(1,r2)
			cat(" ",fill=TRUE)
			cat("Warning: the columns of B will not be reflected", fill=TRUE)
			cat(" ",fill=TRUE)
		} else{	
			if ((sum(tau^2)!=r2) | (sum((tau^2-1)^2)!=0)){
				cat(" ",fill=TRUE)
				cat("Error! Incorrect reflection! The columns of B will not be reflected!", fill=TRUE)
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
					cat("Error! Incorrect permutation! The columns of B will not be permuted!", fill=TRUE)
					cat(" ",fill=TRUE)
					pp=c(1:r2)
				}
			}
		} else{
			pp=1
		}
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
		rownames(K)=labCompA
		colnames(BT)=labCompB
		colnames(K)=corelabelstr
		if (T2==3){
			cat(" ",fill=TRUE)
			cat("Rotated C (CU)",fill=TRUE)
			print(round(CU,digits=2))
		}
		cat(" ",fill=TRUE)
		cat("Rotated core",fill=TRUE)
		print(round(K,digits=2))
	}
	if (T2!=1){	
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
				cat("Error! Incorrect reflection! The columns of C will not be reflected!", fill=TRUE)
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
		if (r3==1){
			K=t(K)
		}
		K=permnew(K,r3,r1,r2)
		rownames(K)=labCompA
		colnames(CU)=labCompC
		colnames(K)=corelabelstr
		cat(" ",fill=TRUE)
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
cat(" ",fill=TRUE)
cat("Rotated solution for A, B and C is in AS, BT and CU, rotated core in K",fill=TRUE)
cat("Only the matrices for reduced modes and the core will be given.",fill=TRUE)
if (T2!=3){
	cat("Rotated A (AS)",fill=TRUE)
	print(round(AS,digits=2))
}
if (T2!=2){
	cat("Rotated B (BT)",fill=TRUE)
	print(round(BT,digits=2))
}
if (T2!=1){
	cat("Rotated C (CU)",fill=TRUE)
	print(round(CU,digits=2))
}
cat("Rotated core",fill=TRUE)
print(round(K,digits=2))

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
	T2fit=T3fitpartitioning(Xprep,n,m,p,AS,BT,CU,K,renormmode,laba,labb,labc)		
	cat("Relative fit of A-mode entities (in %)", fill=TRUE)
	#print(cbind(round(T2fit$fitA,digits=2)),rownames=laba)
	print(round(T2fit$fitA,digits=2))
	cat("Relative fit of B-mode entities (in %)",fill=TRUE)
	#print(cbind(round(T2fit$fitB,digits=2)),rownames=labb)
	print(round(T2fit$fitB,digits=2))
	cat("Relative fit of C-mode entities (in %)",fill=TRUE)
	#print(cbind(round(T2fit$fitC,digits=2)),rownames=labc)
	print(round(T2fit$fitC,digits=2))
	if (renormmode==1){
		cat("Contributions to total fit (in %) of each B- & C-mode component",fill=TRUE)
		print(round(T2fit$BCcontr,digits=2))
	}
	if (renormmode==2){
		cat("Contributions to total fit (in %) of each A- & C-mode component",fill=TRUE)
		print(round(T2fit$ACcontr,digits=2))
	}
	if (renormmode==3){
		cat("Contributions to total fit (in %) of each A- & B-mode component",fill=TRUE)
		print(round(T2fit$ABcontr,digits=2))
	}
	cat("Results for the fitpartitioning are given in $fitA, $fitB, $fitC, $fitAB, $fitAC, $fitBC",fill=TRUE)
} else{
	T2fit=list()
	T2fit$fitA=NULL
	T2fit$fitB=NULL
	T2fit$fitC=NULL
	T2fit$ABcontr=NULL
	T2fit$BCcontr=NULL
	T2fit$ACcontr=NULL
}

cat(" ",fill=TRUE)
cat("Press 'return' to conclude the analysis",fill=TRUE)
c=scan("",n=1)
cat("To see (rotated) component matrices, core matrix, fit and other results: $Xprep, $A, $B, $C, $K, $fit, $fitA, $fitB, $fitC, ...", fill=TRUE)

XprepOut=rarray(Xprep,n,m,p)
dimnames(XprepOut)=list(laba,labb,labc)
strCB=noquote(vector(mode="character",length=p*m))
i=1
for (k in 1:p){
	for (j in 1:m){
		strCB[i]=noquote(paste(" B",as.character(j),"xC",as.character(k),sep=""))
		i=i+1
	}
}
rownames(CBplot)=strCB
colnames(CBplot)=labCompA
strAC=noquote(vector(mode="character",length=n*p))
i=1
for (k in 1:n){
	for (j in 1:p){
		strAC[i]=noquote(paste(" C",as.character(j),"xA",as.character(k),sep=""))
		i=i+1
	}
}
rownames(ACplot)=strAC
colnames(ACplot)=labCompB
strBA=noquote(vector(mode="character",length=m*n))
i=1
for (k in 1:m){
	for (j in 1:n){
		strBA[i]=noquote(paste(" A",as.character(j),"xB",as.character(k),sep=""))
		i=i+1
	}
}
rownames(BAplot)=strBA
colnames(BAplot)=labCompC

out=list()
out$A=AS
out$B=BT
out$C=CU
out$core=K
out$fit=fp
out$fitValues=funcV
out$funcValues=funcF
out$cputime=cputime
out$iter=niter
out$fitA=T2fit$fitA
out$fitB=T2fit$fitB
out$fitC=T2fit$fitC
out$fitAB=T2fit$fitABcontr
out$fitAC=T2fit$fitACcontr
out$fitBC=T2fit$fitBCcontr
out$laba=laba
out$labb=labb
out$labc=labc
out$Xprep=XprepOut
return(out)
}

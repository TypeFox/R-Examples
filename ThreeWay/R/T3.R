T3 <-
function(data,laba,labb,labc){

cat(" ",fill=TRUE)
cat("WELCOME to the interactive TUCKER3 analysis program",fill=TRUE)
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
	if ((prep!=0) & (prep!=1) & (prep!=2) & (prep!=3)){
		cat(" ",fill=TRUE)
		cat("Error! Make a proper choice for normalizing the data",fill=TRUE)
		cat(" ",fill=TRUE)
	}
}

cat(" ",fill=TRUE)
cat("Note: The preprocessed data are now available in Xprep.",fill=TRUE)
cat("Xprep can be used for analyses outside the Tucker3 program",fill=TRUE)

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
	cat("These matrices can be used for analyses outside the Tucker3 program and are given in $A1, $A2, $B1, $B2, $C1, $C2",fill=TRUE)
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
cat("You are now about to do a Tucker3 analysis.",fill=TRUE)

cat(" ",fill=TRUE)
cat("You have to specify the dimensionalities to use.",fill=TRUE)
cat("To search these, you are advised to first run PCASUP analyses.",fill=TRUE)
cat("(PCASUP analyses are PCA's of supermatrices with slices of the 3way array next to each other,",fill=TRUE)
cat("thus 3 supermatrices are analyed by PCA, and for each we get a component matrix, and eigenvalues.",fill=TRUE)
cat("The eigenvalues may give an indication as to the required number of components for each mode.)",fill=TRUE)
cat("The results can next be used to find useful dimensionalities by the generalized scree test",fill=TRUE)
cat("If you want to do the PCASUPs for choosing your dimensionality, specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
if (c==1){
	cat("For the generalized scree test it is needed to indicate the maximum number of dimensions for each mode you want to study",fill=TRUE)
	cat("Up to how many A-mode components do you want to use?",fill=TRUE)
	check0=0
	while (check0==0){
		maxa=scan("",n=1)
		if (length(maxa)==0){
			maxa=0
		}
		if ((maxa==0) | ((floor(maxa)-maxa)!=0) | (maxa>n)){
			cat(" ",fill=TRUE)
			cat("Error! Up to how many A-mode components do you want to use?",fill=TRUE)
			cat(" ",fill=TRUE)
		} else{
			check0=1
		}
	}	
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
	cat(" ",fill=TRUE)
	out=T3runsApproxFit(Xprep,n,m,p,maxa,maxb,maxc)
	cat(" ",fill=TRUE)
	cat("       Numbers of components   Fit   Total number of",fill=TRUE)
	rownames(out)=rep("T3",length=dim(out)[1])
	colnames(out)=c("      A","      B","      C","    (%) "," components")
	print(round(out,digits=2))
	if (nrow(out)>2){
		cat(" ",fill=TRUE)
		cat("Suggestion based on Convex Hull procedure",fill=TRUE)
		DS=DimSelector(out,n,m,p,2)
	}
	cat("Figure 1 now gives fit values of all solutions plotted against tnc=P+Q+R",fill=TRUE)
	cat("Figure 2 gives fit values (%) of all solutions plotted against number of free parameters",fill=TRUE)
	T3dimensionalityplot(out,n,m,p)
}

cat(" ",fill=TRUE)
cat("You can now do a TUCKER3 ANALYSIS",fill=TRUE)

cat(" ",fill=TRUE)
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
Tuck3=T3func(Xprep,n,m,p,r1,r2,r3,0,conv)		
A=Tuck3$A
B=Tuck3$B
C=Tuck3$C
H=Tuck3$H
La=Tuck3$La
Lb=Tuck3$Lb
Lc=Tuck3$Lc
iter=Tuck3$iter
cput=Tuck3$cputime
fp=Tuck3$fp
f=Tuck3$f
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
		Tuck3b=T3func(Xprep,n,m,p,r1,r2,r3,start,conv)			
		funcF[run+1]=Tuck3b$fp
		funcV[run+1]=Tuck3b$f
		cputime[run+1]=Tuck3b$cputime
		niter[run+1]=Tuck3b$iter
		if (Tuck3b$fp>1.0001*Tuck3$fp){		# if fit more than .01% better is found, replace solution
			A=Tuck3b$A
			B=Tuck3b$B
			C=Tuck3b$C
			H=Tuck3b$H
			f=Tuck3b$f
			fp=Tuck3b$fp
			La=Tuck3b$La
			Lb=Tuck3b$Lb
			Lc=Tuck3b$Lc
		}
	}
	cat(" ",fill=TRUE)
	cat("Fit (%) values from all runs:",fill=TRUE)
	print(round(funcF,digits=2))
}
names(fp)="Fit (%)"
cat(" ",fill=TRUE)
cat(paste("Tucker3 analysis with ",r1,"x",r2,"x",r3," components, gave a fit of ",round(fp,digits=2), "%"),fill=TRUE) 

out=list()

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
		cat("Which mode do you want to renormalize? Enter 1 (=A), 2 (=B) or 3 (=C):",fill=TRUE)
		renormmode=scan("",n=1)
		if (length(renormmode)==0){
			renormmode=0
		}
		if ((renormmode==1) | (renormmode==2) | (renormmode==3)){
			check1=1
			RNsol=renormsolT3(A,B,C,H,renormmode)
			A=RNsol$A
			B=RNsol$B
			C=RNsol$C
			H=RNsol$H
		} else{
			cat(" ",fill=TRUE)
			cat("Error! Select a proper number!")
			cat(" ",fill=TRUE)
		}
	}
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
	if (r1>1){
		cat("Specify (range of) relative weight(s) for A (default=0):",fill=TRUE)
		aa=scan()
	} else{	
		cat(" ",fill=TRUE)
		cat("Warning: as the number of A-mode components is 1, no simple structure for A will be sought (relative weight=0)",fill=TRUE)
		cat(" ",fill=TRUE)
		aa=0
	}	
	if (r2>1){
		cat("Specify (range of) relative weight(s) for B (default=0):",fill=TRUE)
		bb=scan()
	} else{	
		cat(" ",fill=TRUE)
		cat("Warning: as the number of B-mode components is 1, no simple structure for B will be sought (relative weight=0)",fill=TRUE)
		cat(" ",fill=TRUE)
		bb=0
	}	
	if (r3>1){
		cat("Specify (range of) relative weight(s) for C (default=0):",fill=TRUE)
		cc=scan()
	}
	else{
		cat(" ",fill=TRUE)
		cat("Warning: as the number of C-mode components is 1, no simple structure for C will be sought (relative weight=0)",fill=TRUE)
		cat(" ",fill=TRUE)
		cc=0
	}	
	if (length(aa)==0){
		aa=0
	}
	if (length(bb)==0){
		bb=0
	}
	if (length(cc)==0){
		cc=0
	}
	for (wa_rel in aa){
		for (wb_rel in bb){
			for (wc_rel in cc){
				cat(" ",fill=TRUE)
				cat(paste("Simple structure rotation with relative weights ", wa_rel, wb_rel, wc_rel),fill=TRUE)
				VAR=varimcoco(A,B,C,H,wa_rel,wb_rel,wc_rel)		
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
	if (max(r1,r2,r3)==1){
		c=0
	} else{	
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
	if (max(r1,r2,r3)>1){
		cat("To specify a different simple structure rotation, specify '1':",fill=TRUE)
		c=scan("",n=1)
		if (length(c)==0){
			c=0
		}
		if (c==1){
			cat("Specify your choice for the preferred simple structure rotation of core and components",fill=TRUE)
			if (r1>1){
				cat("Specify relative weight for simplicity of A-mode (default=0):",fill=TRUE)
				wa_rel=scan("",n=1)
				if (length(wa_rel)==0){
					wa_rel=0
				}
			} else{	
				wa_rel=0
			}	
			if (r2>1){
				cat("Specify relative weight for simplicity of B-mode (default=0):",fill=TRUE)
				wb_rel=scan("",n=1)
				if (length(wb_rel)==0){
					wb_rel=0
				}
			} else{	
				wb_rel=0
			}	
			if (r3>1){
				cat("Specify relative weight for simplicity of C-mode (default=0):",fill=TRUE)
				wc_rel=scan("",n=1)
				if (length(wc_rel)==0){
					wc_rel=0
				}
			} else{	
				wc_rel=0
			}	
			VAR=varimcoco(A,B,C,H,wa_rel,wb_rel,wc_rel)		
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
		rownames(CU)=labc
		colnames(CU)=labCompC
		colnames(K)=corelabelstr
	}
	cat("Rotated solution for A, B and C is in AS, BT and CU, rotated core in K",fill=TRUE)
	cat("Note: Solution is in AS, BT, CU and K even if no simple structure rotation was carried out.",fill=TRUE)
	cat("Rotated A (AS)",fill=TRUE)
	print(round(AS,digits=2))
	cat("Rotated B (BT)",fill=TRUE)
	print(round(BT,digits=2))
	cat("Rotated C (CU)",fill=TRUE)
	print(round(CU,digits=2))
	cat("Rotated core",fill=TRUE)
	print(round(K,digits=2))

	if (max(r1,r2,r3)>1){
		cat(" ",fill=TRUE)
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
		cat("Note: Solution will be in AS, BT, CU and K even if no simple structure rotation was carried out.",fill=TRUE)
		cat("Rotated A (AS)",fill=TRUE)
		print(round(AS,digits=2))
		cat("Rotated B (BT)",fill=TRUE)
		print(round(BT,digits=2))
		cat("Rotated C (CU)",fill=TRUE)
		print(round(CU,digits=2))
		cat("Rotated core",fill=TRUE)
		print(round(K,digits=2))
	}
	disp=1
	cat("Give a vector for reflection of columns of A (e.g., 1 -1 -1 1 ..)",fill=TRUE)
	tau=scan("",n=r1)
	if (length(tau)==0){
		cat(" ",fill=TRUE)
		cat("Warning: the columns of A will not be reflected", fill=TRUE)
		cat(" ",fill=TRUE)
		tau=array(1,r1)
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
	colnames(AS)=labCompA
	rownames(K)=labCompA
	colnames(K)=corelabelstr
	cat(" ",fill=TRUE)
	cat("Rotated B (BT)",fill=TRUE)
	print(round(BT,digits=2))
	cat("Rotated C (CU)",fill=TRUE)
	print(round(CU,digits=2))
	cat("Rotated K (K)",fill=TRUE)
	print(round(K,digits=2))
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
				cat("Incorrect permutation", fill=TRUE)
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
	colnames(BT)=labCompB
	rownames(K)=labCompA
	colnames(K)=corelabelstr
	cat(" ",fill=TRUE)
    cat("Rotated C (CU)",fill=TRUE)
	print(round(CU,digits=2))
	cat(" ",fill=TRUE)
	cat("Rotated core",fill=TRUE)
	print(round(K,digits=2))
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
	colnames(CU)=labCompC
	rownames(K)=labCompA
	colnames(K)=corelabelstr
	cat(" ",fill=TRUE)
	cat("Rotated C (CU)",fill=TRUE)
	print(round(CU,digits=2))
	cat("Rotated core",fill=TRUE)
	print(round(K,digits=2))
	cat(" ",fill=TRUE)
	cat("If you want to further reflect/permute columns/rows, specify '1'",fill=TRUE)
	c=scan("",n=1)
	if (length(c)==0){
		c=0
	}
}
cat(" ",fill=TRUE)
cat("Rotated solution for A, B and C is in AS, BT and CU, rotated core in K",fill=TRUE)
cat("Rotated A (AS)",fill=TRUE)
rownames(AS)=laba
print(round(AS,digits=2))
cat("Rotated B (BT)",fill=TRUE)
rownames(BT)=labb
print(round(BT,digits=2))
cat("Rotated C (CU)",fill=TRUE)
rownames(CU)=labc
print(round(CU,digits=2))
cat("Rotated core",fill=TRUE)
rownames(K)=labCompA
colnames(K)=corelabelstr
print(round(K,digits=2))

cat(" ",fill=TRUE)
cat("You can produce a JOINT PLOT now",fill=TRUE)
cat("In joint plots, one mode is fixed, while the entities of the other two are plotted",fill=TRUE)
cat("If you want to produce a joint plot for the current solution, specify '1':",fill=TRUE)
c=scan("",n=1)
if (length(c)==0){
	c=0
}
while (c==1){
	cat("Which mode do you want to keep fixed ('1' for A-mode, '2' for B-mode or '3' for C-mode)?",fill=TRUE)
	fixmode=scan("",n=1)
	if ((fixmode==1) | (fixmode==2) | (fixmode==3)){
		cat("For which component do you want to see the joint plot?",fill=TRUE)
		fixunit=scan("",n=1)
			if (((fixmode==1) & (fixunit<=r1) & (floor(fixunit)==fixunit)) | ((fixmode==2) & (fixunit<=r2) & (floor(fixunit)==fixunit)) | ((fixmode==3) & (fixunit<=r3))  & (floor(fixunit)==fixunit)){
				JTplot=jointplotgen(K,AS,BT,CU,fixmode,fixunit,laba,labb,labc)
				cat(paste("The plot displays",round(JTplot,digits=2),"of the information represented by the component of your choice"),fill=TRUE)
			} else{
				cat(" ",fill=TRUE)
				cat("Error! Number of component specifyed wrongly! No joint plot produced!",fill=TRUE)
				cat(" ",fill=TRUE)
			}
		} else{
		cat(" ",fill=TRUE)
		cat("Error: Fixed mode specified wrongly! No joint plot produced!",fill=TRUE)
		cat(" ",fill=TRUE)
	}
	c=0
	cat(" ",fill=TRUE)
	cat("Do you want to produce another joint plot? If so, type '1':",fill=TRUE)
	c=scan("",n=1)
	if (length(c)==0){
		c=0
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
		SPLIT=splithalfT3(Xprep,n,m,p,r1,r2,r3,centopt,normopt,renormmode,wa_rel,wb_rel,wc_rel,addanal,conv,laba,labb,labc)
	} else{
		split=0
	}	
}
if (split==1){
	cat("Note: for advanced analysis one can further study the results of the splits and of the full data set",fill=TRUE)
	cat("Results for A, B, C, and Core of full data are given in $Afull, $Bfull, $Cfull and $Kfull",fill=TRUE)
	cat("Results for splits are given in $As1, $Bs1, $Cs1, $Ks1 and $As2, $Bs2, $Cs2, $Ks2",fill=TRUE)
	cat("Results for core in weak stability check are given in $Hss1 and $Hss2",fill=TRUE)
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
	SPLIT$Kfull=NULL
	SPLIT$Ks1=NULL
	SPLIT$Ks2=NULL
	SPLIT$Kss1=NULL
	SPLIT$Kss2=NULL
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
	A=AS
	B=BT
	C=CU
	G=K
	cat("Procedure uses optimal matching of bootstrap solutions to original solution",fill=TRUE)
	check1=0
	while(check1==0){
		cat("How do you like to arrive towards full solutions?",fill=TRUE)
		cat("0 = matching via orthogonal rotation towards full solutions",fill=TRUE)
		cat("1 = matching via optimal transformation towards full solutions",fill=TRUE)
		optimalmatch=scan("",n=1)
		if (length(optimalmatch)==0){
			optimalmatch=2
		}
		if ((optimalmatch==0) | (optimalmatch==1)){
			check1=1
			BOOT=bootstrapT3(Xprep,A,B,C,G,n,m,p,r1,r2,r3,conv,centopt,normopt,optimalmatch,laba,labb,labc)	
		} else{
			if ((optimalmatch!=0) & (optimalmatch!=1)){
				cat(" ",fill=TRUE)
				cat("Error! Specify a proper choice!",fill=TRUE)
				cat(" ",fill=TRUE)
			}
		}
	}

	# make labels for intervals of columns of B and C
	cat(" ",fill=TRUE)
	cat("Bootstrap confidence intervals for fit percentage",fill=TRUE)
	print(round(BOOT$fpint,digits=2))
	cat("Bootstrap confidence intervals for B, per component next to each other",fill=TRUE)
	print(round(BOOT$Bint,digits=2))
	cat("Bootstrap confidence intervals for C, per component next to each other",fill=TRUE)
	print(round(BOOT$Cint,digits=2))
	cat("Bootstrap confidence intervals for Core",fill=TRUE)
	print(round(BOOT$Gint,digits=2))
	cat("Results for the bootstrap confidence intervals are given in $Bint, $Cint, $Kint, $FITint",fill=TRUE)
} else{
	BOOT=list()
	BOOT$Bint=NULL
	BOOT$Cint=NULL
	BOOT$Kint=NULL
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
	cat("These contributions can simply be added to get aggregate contributions in case all components are orthogonal",fill=TRUE)
	cat("If this is not the case, aggregate fit contributions are given separately, below",fill=TRUE)
	print(round((K^2)*100/sum(Xprep^2),digits=2))
	T3fit=T3fitpartitioning(Xprep,n,m,p,AS,BT,CU,K,renormmode,laba,labb,labc)		
	cat("Relative fit of A-mode entities (in %)", fill=TRUE)
	print(round(T3fit$fitA,digits=2))
	cat("Relative fit of B-mode entities (in %)",fill=TRUE)
	print(round(T3fit$fitB,digits=2))
	cat("Relative fit of C-mode entities (in %)",fill=TRUE)
	print(round(T3fit$fitC,digits=2))
	if (renormmode==1){
		cat("Contributions to total fit (in %) of each B- & C-mode component",fill=TRUE)
		print(round(T3fit$BCcontr,digits=2))
	}
	if (renormmode==2){
		cat("Contributions to total fit (in %) of each A- & C-mode component",fill=TRUE)
		print(round(T3fit$ACcontr,digits=2))
	}
	if (renormmode==3){
		cat("Contributions to total fit (in %) of each A- & B-mode component",fill=TRUE)
		print(round(T3fit$ABcontr,digits=2))
	}
	cat("Results for the fitpartitioning are given in $fitA, $fitB, $fitC, $fitAB, $fitAC, $fitBC",fill=TRUE)
} else{
	T3fit=list()
	T3fit$fitA=NULL
	T3fit$fitB=NULL
	T3fit$fitC=NULL
	T3fit$ABcontr=NULL
	T3fit$BCcontr=NULL
	T3fit$ACcontr=NULL
}



cat(" ",fill=TRUE)
cat("Press 'return' to conclude the analysis",fill=TRUE)
c=scan("",n=1)
cat("To see (rotated) component matrices, core matrix, fit and other results: $Xprep, $A, $B, $C, $core, $fit, $fitA, $fitB, $fitC, ...", fill=TRUE)

XprepOut=rarray(Xprep,n,m,p)
dimnames(XprepOut)=list(laba,labb,labc)
rownames(Aplot)=laba
rownames(Bplot)=labb
rownames(Cplot)=labc
colnames(Aplot)=labCompA
colnames(Bplot)=labCompB
colnames(Cplot)=labCompC
strCB=noquote(vector(mode="character",length=p*m))
i=1
for (j in 1:p){
	for (k in 1:m){
		strCB[i]=noquote(paste(" c",as.character(j),"xb",as.character(k),sep=""))
		i=i+1
	}
}
rownames(CBplot)=strCB
colnames(CBplot)=labCompA

strAC=noquote(vector(mode="character",length=n*p))
i=1
for (j in 1:n){
	for (k in 1:p){
		strAC[i]=noquote(paste(" a",as.character(j),"xc",as.character(k),sep=""))
		i=i+1
	}
}
rownames(ACplot)=strAC
colnames(ACplot)=labCompB

strBA=noquote(vector(mode="character",length=m*n))
i=1
for (j in 1:m){
	for (k in 1:n){
		strBA[i]=noquote(paste(" b",as.character(j),"xa",as.character(k),sep=""))
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
out$fitA=T3fit$fitA
out$fitB=T3fit$fitB
out$fitC=T3fit$fitC
out$fitAB=T3fit$fitABcontr
out$fitAC=T3fit$fitACcontr
out$fitBC=T3fit$fitBCcontr
out$Bint=BOOT$Bint
out$Cint=BOOT$Cint
out$Kint=BOOT$Gint
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
out$Kfull=SPLIT$Kfull
out$Ks1=SPLIT$Ks1
out$Ks2=SPLIT$Ks2	
out$Kss1=SPLIT$Hss1
out$Kss2=SPLIT$Hss2
out$Aplot=Aplot 
out$Bplot=Bplot
out$Cplot=Cplot
out$CBplot=CBplot
out$ACplot=ACplot
out$BAplot=BAplot
out$A1=PCAmean$A1
out$A2=PCAmean$A2
out$B1=PCAmean$B1
out$B2=PCAmean$B2
out$C1=PCAmean$C1
out$C2=PCAmean$C2
out$laba=laba
out$labb=labb
out$labc=labc
out$Xprep=XprepOut
return(out)
}

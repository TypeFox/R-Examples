splithalfT3 <-
function(X,n,m,p,r1,r2,r3,centopt,normopt,renormmode,wa_rel,wb_rel,wc_rel,addanal,conv,laba,labb,labc){

cat("This procedure performs a SPLIT-HALF analysis on X",fill=TRUE)
cat("NOTE:",fill=TRUE)
cat("In SPLIT-HALF analysis, the A-mode is taken as 'replication mode'",fill=TRUE)
cat("(which means that A-mode entities are considered a random sample",fill=TRUE)
cat("if this does not make sense, you should rearrange your data so that the A-mode is a replication mode)",fill=TRUE)
cat("The splitting into two halves can be done randomly (default), or into odd vs. even sequence numbers",fill=TRUE)

narg=nargs()
if (narg<16){
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
cat("However, you may want to modify certain choices here (rather than rerunning the full Tucker3)",fill=TRUE)
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
				cat("Error! Make a proper choice for centering the data")
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
				cat("Error! Make a proper choice for normalizing the data")
				cat(" ",fill=TRUE)
			} else{
				check=1
			}
		}
	}

	if (renormmode==1){
		cat("Your choice was to renormalize A-mode",fill=TRUE)
	}
	if (renormmode==2){
		cat("Your choice was to renormalize B-mode",fill=TRUE)
	}
	if (renormmode==3){
		cat("Your choice was to renormalize C-mode",fill=TRUE)
	} else{
		if ((renormmode!=1) & (renormmode!=2) & (renormmode!=3)){
			cat("Your choice was not to renormalize solution",fill=TRUE)
		}
	}
	cat("If you want to change it, specify '1':",fill=TRUE)
	cc=scan("",n=1)
	if (length(cc)==0){
		cc=0
	}
	if (cc==1){
		check2=0
		while (check2==0){
			cat("Which mode do you want to renormalize? Enter 1 (=A), 2 (=B) or 3 (=C):", fill=TRUE)
			renormmode=scan("",n=1)
			if (length(renormmode)==0){
				renormmode=0
			}
			if ((renormmode!=1) &(renormmode!=2) & (renormmode!=3)){
				cat(" ",fill=TRUE)
				cat("Error! Make a proper choice for normalizing the data")
				cat(" ",fill=TRUE)
			} else{
				check2=1
			}
		}
	}
		
  	cat(paste("Numbers of components for A, B and C were: ",r1,", ",r2,", ",r3,sep=""),fill=TRUE)
	cat("If you want to change them, specify '1':",fill=TRUE)
	cc=scan("",n=1)
	if (length(cc)==0){
		cc=0
	}
	if (cc==1){
		cat("How many A-mode components do you want to use?",fill=TRUE)
		check=0
		while (check==0){
			r1=scan("",n=1)
			if (length(r1)==0){
				r1=0
			}
			if ((r1==0) | ((floor(r1)-r1)!=0)){
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
			if ((r2==0) | ((floor(r2)-r2)!=0)){
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
			if ((r3==0) | ((floor(r3)-r3)!=0)){
				cat(" ",fill=TRUE)
				cat("Error! How many C-mode components do you want to use?",fill=TRUE)
				cat(" ",fill=TRUE)
			} else{
				check=1
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
	
	cat(paste("Relative simplicity rotation weights for A, B and C are: ",wa_rel,", ",wb_rel,", ",wc_rel,sep=""),fill=TRUE)
	cat("If you want to change them, specify '1':",fill=TRUE)
	cc=scan("",n=1)
	if (length(cc)==0){
		cc=0
	}
	if (cc==1){
		cat("Specify relative weight for simplicity of A-mode (default=0):",fill=TRUE)
		wa_rel=scan("",n=1)
		cat("Specify relative weight for simplicity of B-mode (default=0):",fill=TRUE)
		wb_rel=scan("",n=1)
		cat("Specify relative weight for simplicity of C-mode (default=0):",fill=TRUE)
		wc_rel=scan("",n=1)
		if (length(wa_rel)==0){
			wa_rel=0
		}
		if (length(wb_rel)==0){
			wb_rel=0
		}
		if (length(wc_rel)==0){
			wc_rel=0
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
Tuck3=T3func(X,n,m,p,r1,r2,r3,0,1e-6)
A=Tuck3$A
B=Tuck3$B
C=Tuck3$C
H=Tuck3$H
f=Tuck3$f
iter=Tuck3$iter
fp=Tuck3$fp
La=Tuck3$La
Lb=Tuck3$Lb
Lc=Tuck3$Lc
func=vector("numeric",length=1+addanal)
names(func)=paste("Start n.",1:(1+addanal),sep="")
func=Tuck3$fp
print(func)
for (run in 1:addanal){
    cat(paste("Run no.",run+1,sep=" "),fill=TRUE)
	Tuck3b=T3func(X,n,m,p,r1,r2,r3,1,1e-6)
	func[run+1]=Tuck3b$fp
	if (Tuck3b$fp>1.0001*Tuck3$fp){		# if fit more than .01% better is found, replace solution
		A=Tuck3b$A
		B=Tuck3b$B
		C=Tuck3b$C
		H=Tuck3b$H
		f=Tuck3b$f
		iter=Tuck3b$iter
		fp=Tuck3b$fp
		La=Tuck3b$La
		Lb=Tuck3b$Lb
		Lc=Tuck3b$Lc
    }
}
if (addanal>=1){
	cat("Fit % values from all runs:",fill=TRUE)
	print(round(t(func), digits=2))
}
RNsol=renormsolT3(A,B,C,H,renormmode) 	
# Full data Rotation
VARM=varimcoco(RNsol$A,RNsol$B,RNsol$C,RNsol$H,wa_rel,wb_rel,wc_rel)
AS=VARM$AS
BT=VARM$BT
CU=VARM$CU
K=VARM$K
cat("Backnormalize A,B,C, and H",fill=TRUE)
if (renormmode==1){
	Ds=diag(SUM(AS)$col^.5,nrow=r1)
	AS=AS%*%solve(Ds)
	K=Ds%*%K
}
if (renormmode==2){
	Ds=diag(SUM(BT)$col^.5,nrow=r2)
	BT=BT%*%solve(Ds)
	K=K%*%kronecker(diag(r3),Ds)
}
if (renormmode==3){
	Ds=diag(SUM(CU)$col^.5,nrow=r3)
	CU=CU%*%solve(Ds)
	K=K%*%kronecker(Ds,diag(r2))
}
Afull=AS
Bfull=BT
Cfull=CU
Kfull=K

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
Tuck3=T3func(X,n,m,p,r1,r2,r3,0,1e-6)				
A=Tuck3$A
B=Tuck3$B
C=Tuck3$C
H=Tuck3$H
f=Tuck3$f
iter=Tuck3$iter
fp=Tuck3$fp
La=Tuck3$La
Lb=Tuck3$Lb
Lc=Tuck3$Lc
func=vector("numeric",length=1+addanal)
names(func)=paste("Start n.",1:(1+addanal),sep="")
func=Tuck3$fp
for (run in 1:addanal){
    cat(paste("Run no.",run+1,sep=" "),fill=TRUE)
	Tuck3b=T3func(X,n,m,p,r1,r2,r3,1,1e-6)			
	func[run+1]=Tuck3b$fp
	if (Tuck3b$fp>1.0001*Tuck3$fp){		# if fit more than .01% better is found, replace solution
		A=Tuck3b$A
		B=Tuck3b$B
		C=Tuck3b$C
		H=Tuck3b$H
		f=Tuck3b$f
		iter=Tuck3b$iter
		fp=Tuck3b$fp
		La=Tuck3b$La
		Lb=Tuck3b$Lb
		Lc=Tuck3b$Lc
    }
}
if (addanal>=1){
	cat("Fit % values from all runs:",fill=TRUE)
	print(round(t(func), digits=2))
}
cat("No Postprocessing used in splits, taken care off by reference rotation!",fill=TRUE)
As1=A
Bs1=B
Cs1=C
Ks1=H

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
Tuck3=T3func(X,n,m,p,r1,r2,r3,0,1e-6)				
A=Tuck3$A
B=Tuck3$B
C=Tuck3$C
H=Tuck3$H
f=Tuck3$f
iter=Tuck3$iter
fp=Tuck3$fp
La=Tuck3$La
Lb=Tuck3$Lb
Lc=Tuck3$Lc
func=vector("numeric",length=1+addanal)
names(func)=paste("Start n.",1:(1+addanal),sep="")
func=Tuck3$fp
for (run in 1:addanal){
	cat(paste("Run no.",run+1,sep=" "),fill=TRUE)
	Tuck3b=T3func(X,n,m,p,r1,r2,r3,1,1e-6)			
	func[run+1]=Tuck3b$fp
	if (Tuck3b$fp>1.0001*Tuck3$fp){		# if fit more than .01% better is found, replace solution
		A=Tuck3b$A
		B=Tuck3b$B
		C=Tuck3b$C
		H=Tuck3b$H
		f=Tuck3b$f
		iter=Tuck3b$iter
		fp=Tuck3b$fp
		La=Tuck3b$La
		Lb=Tuck3b$Lb
		Lc=Tuck3b$Lc
	}
}
if (addanal>=1){
	cat("Fit % values from all runs:",fill=TRUE)
	print(round(t(func), digits=2))
}

cat("No Postprocessing used in splits, taken care off by reference rotation!",fill=TRUE)
As2=A
Bs2=B
Cs2=C
Ks2=H

n=n_orig
X=X_orig

# Compute stability indices:
Ss1=solve(t(As1)%*%As1)%*%t(As1)%*%Afull[wi[1:n1],]		
Ss2=solve(t(As2)%*%As2)%*%t(As2)%*%Afull[wi[(n1+1):n],]
Ts1=solve(t(Bs1)%*%Bs1)%*%t(Bs1)%*%Bfull
Ts2=solve(t(Bs2)%*%Bs2)%*%t(Bs2)%*%Bfull
Us1=solve(t(Cs1)%*%Cs1)%*%t(Cs1)%*%Cfull
Us2=solve(t(Cs2)%*%Cs2)%*%t(Cs2)%*%Cfull

As1=As1%*%Ss1
As2=As2%*%Ss2
Bs1=Bs1%*%Ts1
Bs2=Bs2%*%Ts2
Cs1=Cs1%*%Us1
Cs2=Cs2%*%Us2
Ks1=solve(Ss1)%*%Ks1%*%solve(kronecker(t(Us1),t(Ts1)))		
Ks2=solve(Ss2)%*%Ks2%*%solve(kronecker(t(Us2),t(Ts2)))

labCompA=paste("A",1:r1,sep="")
labCompB=paste("B",1:r2,sep="")
labCompC=paste("C",1:r3,sep="")
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

cat("RESULTS stability analysis",fill=TRUE)
cat("Both splits are analyzed and rotated towards solution for full data",fill=TRUE)
cat("One can compare complete outcomes for different splits",fill=TRUE)
cat("or inspect correlation/congruence coefficients between corresponding columns of component matrices",fill=TRUE)
cat("Split-half solutions for B",fill=TRUE)
cat("B's next to each other, separated by column of 0's",fill=TRUE)
X=round(cbind(Bs1,0,Bs2),digits=2)
rownames(X)=labb
colnames(X)=c(labCompB,"-",labCompB)
print(X)
cat("Split-half solutions for C",fill=TRUE)
cat("C's next to each other, separated by column of 0's",fill=TRUE)
X=round(cbind(Cs1,0,Cs2),digits=2)
rownames(X)=labc
colnames(X)=c(labCompC,"-",labCompC)
print(X)
cat("Congruences for A in splits and in appropriate part of Afull",fill=TRUE)
X=round(cbind(diag(phi(Afull[wi[1:n1],],As1)),diag(phi(Afull[wi[(n1+1):n],],As2))),digits=2)
rownames(X)=labCompA
colnames(X)=c("SPL1","SPL2")
print(X)
cat("Correlations for A in splits and in appropriate part of Afull",fill=TRUE)
X=round(cbind(diag(phi(Cc(Afull[wi[1:n1],]),Cc(As1))),diag(phi(Cc(Afull[wi[(n1+1):n],]),Cc(As2)))),digits=2)
rownames(X)=labCompA
colnames(X)=c("SPL1","SPL2")
print(X)
cat("Congruence values for B-mode component matrix",fill=TRUE)
X=round(diag(phi(Bs1,Bs2)),digits=2)
names(X)=c(labCompB)
print(X)
cat("Congruence values for C-mode component matrix",fill=TRUE)
X=round(diag(phi(Cs1,Cs2)),digits=2)
names(X)=c(labCompC)
print(X)
cat("Relatively strong stability check of Core:",fill=TRUE)
cat("(simply comparing core matrices for two splits)",fill=TRUE)
cat("Core for split1",fill=TRUE)
Ks1=as.matrix(Ks1,digits=2)
rownames(Ks1)=labCompA
colnames(Ks1)=corelabelstr
print(round(Ks1,digits=2))
cat("Core for split2",fill=TRUE)
Ks2=as.matrix(Ks2,digits=2)
rownames(Ks2)=labCompA
colnames(Ks2)=corelabelstr
print(round(Ks2,digits=2))

cat("Weaker but Sufficient stability check of Core:",fill=TRUE)
cat("(computing cores in two splits, using full data solutions for A,B and C)",fill=TRUE)
A1=Afull[wi[1:n1],]
A2=Afull[wi[(n1+1):n],]
Hss1=solve(t(A1)%*%A1)%*%t(A1)%*%Xs1
Hss1=solve(t(Bfull)%*%Bfull)%*%t(Bfull)%*%permnew(Hss1,r1,m,p)
Hss1=solve(t(Cfull)%*%Cfull)%*%t(Cfull)%*%permnew(Hss1,r2,p,r1)
Hss1=permnew(Hss1,r3,r1,r2)
Hss2=solve(t(A2)%*%A2)%*%t(A2)%*%Xs2
Hss2=solve(t(Bfull)%*%Bfull)%*%t(Bfull)%*%permnew(Hss2,r1,m,p)
Hss2=solve(t(Cfull)%*%Cfull)%*%t(Cfull)%*%permnew(Hss2,r2,p,r1)
Hss2=permnew(Hss2,r3,r1,r2)
cat("Core for split1",fill=TRUE)
Hss1=as.matrix(Hss1)
rownames(Hss1)=labCompA
colnames(Hss1)=corelabelstr
print(round(Hss1,digits=2))
cat("Core for split2",fill=TRUE)
Hss2=as.matrix(Hss2)
rownames(Hss2)=labCompA
colnames(Hss2)=corelabelstr
print(round(Hss2,digits=2))
cat("Core for full data",fill=TRUE)
Kfull=as.matrix(Kfull)
rownames(Kfull)=labCompA
colnames(Kfull)=corelabelstr
print(round(Kfull,digits=2))

colnames(As1)=labCompA
colnames(As2)=labCompA
colnames(Afull)=labCompA
colnames(Bs1)=labCompB
colnames(Bs2)=labCompB
colnames(Bfull)=labCompB
colnames(Cs1)=labCompC
colnames(Cs2)=labCompC
colnames(Cfull)=labCompC
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
out$Kfull=Kfull
out$Ks1=Ks1
out$Ks2=Ks2	
out$Kss1=Hss1
out$Kss2=Hss2
return(out)
}

dsFC <-
function(X,Crit,dim=NA){
	NR=nrow(X) #number of respondents
	N=ncol(X) #number of multipleChoice items
	if(Crit>N){stop('\n dsHELP: Item Criterium bigger than ', N)}
	XX<-X
	X<-matrix(as.matrix(X),nrow=NR,ncol=N)
	K<-5200 
	Crit<-Crit # the criterion
  #
  itemOp<-c(1:N)	
	for(j in 1:N){itemOp[j]<-max((X[,j]),na.rm=TRUE)
		}
	M<-sum(itemOp) # total number of options, M
	#
	if(is.na(dim)) 
		solMC<-(M-N)
	else solMC<-dim # Number of solution for MC
  #
SolFC<-(itemOp[Crit]-1) # the number of proper solutions

#########
# Checking the data
#
	if(sum(unlist(apply(X,2,tabulate))==0)>0){stop('\n dsHELP: codification problems. Consider use dsCHECK ')}
	if(sum(is.na(X))>0){stop('\n dsHELP: You have NA values')}
#
	dimnames(X)[[2]]<-paste("q",1:N,sep="")
	cn<-dimnames(X)[[2]]
	ItOpNa<-9
	for(j in 1:N){
		ION<-paste(cn[j], 1:itemOp[j], sep=":")
		ItOpNa<-c(ItOpNa,ION)	
	}
	ItOpNa<-ItOpNa[-1]
	dimnames(X)[[1]]<-paste("s",1:NR, sep="")
	SubNa<-dimnames(X)[[1]]
#
dimnames(X)[[2]]<-paste("q.",1:N,sep="")
	cn<-dimnames(X)[[2]]
#finding Fk and F
#	
	C<-matrix(0,M,M) # a cero matrix to be filled MxM
	#
	Fk<-array(99, dim=c(NR,M)) # F with k
	for(i in 1:NR){ 
		Fe<-c(99) 
		for(j in 1:N){ 
			F1<-array(0,itemOp[j])
			if(is.na(X[i,j])){F1[X[i,j]]<-0}#if NA, 0		
			ifelse(j==Crit,F1[X[i,j]]<-K,F1[X[i,j]]<-1) #if not, K or 1
			Fe<-c(Fe,F1)}
		Fe<-Fe[-1] 
		C1<-Fe%*%t(Fe)
		C<-C+C1
		Fk[i,]<-Fe}
	# 
	F<-array(99, dim=c(NR,M))  #F without K
		for(i in 1:NR){  
			Fe<-c(99)
			for(j in 1:N){ 
				F1<-array(0,itemOp[j]) 
          		if(is.na(X[i,j])){F1[X[i,j]]<-0}#if NA, 0		
          		F1[X[i,j]]<-1 
				Fe<-c(Fe,F1)}
			Fe<-Fe[-1]
			F[i,]<-Fe}
			#
	S<-C 
	Ci<-S/(N-1+K) #
#
if(solMC>(M-N))
		stop("\n dsHELP: The number of dimensions have to be less than or equal to ",(M-N))
		#
	f_t<-NR*(N-1+K) 
	f_c<-apply(Fk,2,sum) #the column sum
	ffp<-f_c%*%t(f_c) #ff'
	D<-diag(f_c)
	D5<-diag(f_c^.5)
	ones<-array(1,dim=c(M,M))
	CorrTerm<-(D5%*%ones%*%D5)/f_t
#
#getting matrix A
#
	Dm5<-diag(f_c^-.5,M,M)
	Dm.5CDm.5<-Dm5%*%Ci%*%Dm5	
	A<-(Dm.5CDm.5-CorrTerm)
#
#the SVD of A
#
	ola<-svd(A)
	eigenval<-ola$d #eigenvalues of A
nieve<-ola$d#########
	#
	control<-sum(eigenval[1:solMC]<1/10^15)
	if(control>0){
		print(round(nieve,3))
		stop('Your data need less components. Check the eigenvalues above.')}
	#
	singuval<-eigenval^.5
	sve<-ola$u #eigenvectors of A
#
#Options weights
#
	w<-array(99,dim=c(M,(M-N))) 
			for(i in 1:(M-N)){
			w[,i]<-(f_t/t(sve[,i])%*%sve[,i])^.5*sve[,i]
			}
	x=Dm5%*%w #NORMED weights	
sing2<-singuval #ORIGINAL singuval
	xAdj=x%*%diag(singuval[1:(M-N)]) 	#PROJECTED weigths
#
#The correlation matrix. STEP[C]
#
#matrix of item scores+row and column totals JUST FOR PROPER SOLUTIONS
#
	cm1<-array(99, dim=c(NR, N, SolFC)) #step1
		for(j in 1:SolFC){
			for(i in 1:NR){
				p<-F[i,]*xAdj[,j] 
				for(k in 1:N)
				if(is.na(X[i,k])){p[sum(itemOp[1:k])]<-9.9}#this is to avoid NA problems (I)
				p<-p[p!=0]
				for(k in 1:N) 
				if(p[k]==9.9){p[k]<-0}#once it is done, we change 9.9 to 0
	cm1[i,,j]<-t(p)}}
	#
	cmRT<-array(99, dim=c(NR,SolFC))
	for(i in 1:SolFC){
		cmRT[,i]<-apply(cm1[,,i],1,sum)} #row totals
	#
	cm2<-array(99,dim=c(NR,N+1,SolFC)) #step2
	for(i in 1:SolFC){
		cm2[,,i]<-cbind(cm1[,,i],cmRT[,i])} #items+rowTotal	
	#
	cmCT<-array(99, dim=c(N+1,SolFC))
	for(i in 1:SolFC){
		cmCT[,i]<-apply(cm2[,,i],2,sum)} #col totals
		#
	cm3<-array(99,dim=c(NR+1,N+1,SolFC))	#step3 
		for(i in 1:SolFC){
		cm3[,,i]<-rbind(cm2[,,i],cmCT[,i])} #item+rowTot+colTot 
		#
	corM<-array(99,dim=c(N+1,N+1,SolFC))
		for(i in 1:SolFC){
				corM[,,i]<-cor(cm2[,,i])} #corM is the correlation of items+rowT
#
############### redoing it with Fk #############################################
#
cm1k<-array(99, dim=c(NR, N, SolFC)) #step1
		for(j in 1:SolFC){
			for(i in 1:NR){
				p<-Fk[i,]*xAdj[,j] #### HERE!
				for(k in 1:N)
				if(is.na(X[i,k])){p[sum(itemOp[1:k])]<-9.9}#
				p<-p[p!=0]
				for(k in 1:N) 
				if(p[k]==9.9){p[k]<-0}#
	cm1k[i,,j]<-t(p)}}
	cmRTk<-array(99, dim=c(NR,SolFC))
	for(i in 1:SolFC){
		cmRTk[,i]<-apply(cm1k[,,i],1,sum)} #row totals
	cm2k<-array(99,dim=c(NR,N+1,SolFC)) #step2
	for(i in 1:SolFC){
		cm2k[,,i]<-cbind(cm1k[,,i],cmRTk[,i])} #items+rowTotal	
	cmCTk<-array(99, dim=c(N+1,SolFC))
	for(i in 1:SolFC){
		cmCTk[,i]<-apply(cm2k[,,i],2,sum)} #col totals
	cm3k<-array(99,dim=c(NR+1,N+1,SolFC))	#step3 
		for(i in 1:SolFC){
		cm3k[,,i]<-rbind(cm2k[,,i],cmCTk[,i])} #item+rowTot+colTot 
	corMk<-array(99,dim=c(N+1,N+1,SolFC))
		for(i in 1:SolFC){
				corMk[,,i]<-cor(cm2k[,,i])} #corM is the correlation of items+rowT
####################################################################################

	tabE1<-array(99,dim=c(SolFC,N+1))
	for(i in 1:SolFC){
		tabE1[i,]<-(corMk[ ,N+1,i])^2}
	tabE1[,N+1]<-apply(tabE1[,1:N],1,mean)
	tabE3<-array(99,dim=c(SolFC+1,N+1))
	tabE3[,1]<-c(0:SolFC)
	for(i in 1:SolFC){
		tabE3[i+1,]<-(corMk[ ,N+1,i])^2}
	tabE3[,N+1]<-apply(tabE3[,1:N],1,mean)
#		
#Adjusting Eigenvalues
# 		
	NotAdj<-eigenval[1:8]
	AdjEigen<-array(99,itemOp[Crit]-1) #Adjusted Eigenvaluer for Result_A
		for (i in 1:(itemOp[Crit]-1))
			AdjEigen[i]<-(sum(tabE3[i+1,1:N])-1)/(N-1)
	#
	Keep<-AdjEigen
	#
	####PAGE 163				
	labels<-c(paste("Item", 1:N), "Average")
	tabE3<-list(labels,tabE3[2:(SolFC+1),1:(N+1)])		
	tabE2<-array(99,dim=c(N+1,3,SolFC))
		for(i in 1:SolFC){
			tabE2[,1,i]<-apply(cm3[1:NR,,i]^2,2,sum)/eigenval[i]# SSj
			tabE2[,2,i]<-tabE1[i,1:(N+1)] #	 R2
			tabE2[,3,i]<-tabE2[,2,i]^.5  #  R
			}
#
	delta<-AdjEigen/sum(AdjEigen)*100
	cumdelta<-cumsum(delta)/sum(delta)*100
#
	result<-data.frame(Component=c(1:SolFC), Eigenvalue=round(AdjEigen,4), SingValue=sqrt(AdjEigen), Delta=delta, CumDelta=cumdelta)
#
#
#
# results
#
	Table5<-array(99,dim=c(N,4,SolFC))
	for(i in 1:SolFC){
		Tab<-data.frame(Item=c(1:N), SSj=tabE2[1:N,1,i], R2=tabE2[1:N,2,i], R=tabE2[1:N,3,i])
		Tab<-as.matrix(Tab)
		Table5[,,i]<-Tab
		} #table of SSj, R2...	
#############################################
#
complement<-dsMC(X[,-Crit],dim)#the ordinary dsMC 
# of non-criterion items, by ignoring 
#the criterion item (type B)
#
#############################################

########THE COMPLETE SPACE or 
#
cm1c<-array(99, dim=c(NR, N, (M-N))) #step1
		for(j in 1:(M-N)){
			for(i in 1:NR){
				p<-Fk[i,]*xAdj[,j]
				for(k in 1:N)
				if(is.na(X[i,k])){p[sum(itemOp[1:k])]<-9.9}#this is to avoid NA problems (I)
				p<-p[p!=0]
				for(k in 1:N) 
				if(p[k]==9.9){p[k]<-0}#once it is done, we change 9.9 to 0
	cm1c[i,,j]<-t(p)}}
	cmRTc<-array(99, dim=c(NR,(M-N)))
	for(i in 1:(M-N)){
		cmRTc[,i]<-apply(cm1c[,,i],1,sum)} #row totals
	cm2c<-array(99,dim=c(NR,N+1,(M-N))) #step2
	for(i in 1:(M-N)){
		cm2c[,,i]<-cbind(cm1c[,,i],cmRTc[,i])} #items+rowTotal	
	cmCTc<-array(99, dim=c(N+1,(M-N)))
	for(i in 1:(M-N)){
		cmCTc[,i]<-apply(cm2c[,,i],2,sum)} #col totals
	cm3c<-array(99,dim=c(NR+1,N+1,(M-N)))	#step3 
		for(i in 1:(M-N)){
		cm3c[,,i]<-rbind(cm2c[,,i],cmCTc[,i])} #item+rowTot+colTot 
	corMc<-array(99,dim=c(N+1,N+1,(M-N)))
		for(i in 1:(M-N)){
				corMc[,,i]<-cor(cm2c[,,i])} #corM is the correlation of items+rowT
				
###############################################
#
	tabE1<-array(99,dim=c((M-N),N+1))
	for(i in 1:(M-N)){
		tabE1[i,]<-(corMc[ ,N+1,i])^2}
	tabE1[,N+1]<-apply(tabE1[,1:N],1,mean)
	tabE3<-array(99,dim=c((M-N)+1,N+1))
	tabE3[,1]<-c(0:(M-N))
	for(i in 1:(M-N)){
		tabE3[i+1,]<-(corMc[ ,N+1,i])^2}
	tabE3[,N+1]<-apply(tabE3[,1:N],1,mean)
#		
#Adjusting Eigenvalues
# 		
	AdjEigen<-eigenval[1:(M-N)]
	tabE3<-as.data.frame(round(tabE3[2:((M-N)+1),1:(N+1)],3))
	dimnames(tabE3)[[2]]<-c(paste("q.",1:N,sep=""),"Avge")
	#
	pp<-tabE3
	tabE2<-array(99,dim=c(N+1,3,(M-N)))
		for(i in 1:(M-N)){
			tabE2[,1,i]<-apply(cm3c[1:NR,,i]^2,2,sum)/eigenval[i]# SSj
			tabE2[,2,i]<-tabE1[i,1:(N+1)] #	 R2
			tabE2[,3,i]<-tabE2[,2,i]^.5  #  R
			}
#
#	alpha<-(1-(1-AdjEigen)/((N-1)*AdjEigen))
	delta<-AdjEigen/sum(AdjEigen)*100
	cumdelta<-cumsum(delta)/sum(delta)*100
#
	resultc<-data.frame(Component=c(1:(M-N)), Eigenvalue=round(AdjEigen,4), SingValue=sqrt(AdjEigen), Delta=delta, CumDelta=cumdelta)
# results
#
	Table5c<-array(99,dim=c(N,4,(M-N))) #####(para M-N) (8)
	for(i in 1:(M-N)){
		Tab<-data.frame(Item=c(1:N), SSjc=tabE2[1:N,1,i], R2=tabE2[1:N,2,i], R=tabE2[1:N,3,i])
		Tab<-as.matrix(Tab)
		Table5c[,,i]<-Tab}#table of SSj, R2...1,2,3,4,5,6,7,8
		#

for(i in 1:(M-N-SolFC)){#1,2,3,4,5,6
		SSjc<-tabE2[1:N,1,i+SolFC]#3,4,5,6,7,8
		sumSSjc<-sum(SSjc)
		SSj<-complement$ForSolution[,2,i] #1,2,3,4,5,6
		sumSSj<-sum(complement$ForSolution[,2,i])
		Table5c2<-Table5c[,2,i+SolFC]#3,4,5,6,7,8 
		NewSSjc<-(sumSSj/sumSSjc)*Table5c2#[,2,i] #1,2,3,4,5,6
		Table5c[,2,i+SolFC]<-NewSSjc} #table of SSj, R2... 3,4,5,6,7,8
#############################
		
AdjEigen<-Keep 					######### GOING BACK

#Projected Space scores
#
	hasta<-cumsum(itemOp)
	desde<-hasta+1
	desde<-c(1,desde[1:(N-1)])
	FN<-F[,-(desde[Crit]:hasta[Crit])]		
	xN<-x[-(desde[Crit]:hasta[Crit]),] 				#x=standard weights
#
	singuval[1:SolFC]<-AdjEigen[1:SolFC]^.5 #using adjusted eigenvalues & others
	singuvalFC<-singuval[1:SolFC]
	singuval<-c(singuvalFC,complement$singuval)
####
#
	f_r<-apply(F,1,sum) #the row sum			
	NM<-diag(f_r)
	f_tstar<-sum(f_r)
	invNM<-diag(1/f_r)
	y=invNM%*%FN%*%xN%*%diag(1/singuval[1:(M-N)])
yG1=invNM%*%FN%*%xN%*%diag(1/sing2[1:(M-N)]) ### Proyec W for Sub original sing

constant<-diag(SolFC)
	for(j in 1:SolFC){
		constant[j,j]<-diag(f_tstar*matrix.inverse(t(y[,j])%*%NM%*%y[,j]))^.5
	}
	y<-y[,1:SolFC]%*%constant

sing3<-singuval
#standard scores
	yAdj=y%*%diag(singuval[1:SolFC]) #projected scores for Su with adjusted singv
xAdj=x%*%diag(singuval[1:(M-N)]) #projected weight for Op with adjusted singv
#
xAdjC=x%*%diag(sing2[1:(M-N)]) #same but times original singuval
yAdjC=yG1%*%diag(sing2[1:(M-N)]) #same but times original singuval

#Match-Mismatch Tables:
#
	MaMiTf<-array(99,dim=c(itemOp[Crit], itemOp[Crit], SolFC))
	CoPred<-array(99,SolFC)
	limits<-cumsum(itemOp)
	from<-limits[Crit-1]+1
	if(Crit==1)from<-1
	to<-limits[Crit]
	for(k in 1:SolFC){
		Table=cbind(y[,k], X[,Crit])
		TableSort<-Table[order(Table[,1]),]
		lab<-paste("w",1:itemOp[Crit],sep="")
		nl<-tabulate(TableSort[,2])
		index<-x[from:to,k]
		nl<-nl[order(-rank(index))]
		Weight<-c(99)
		for(i in itemOp[Crit]:1){
			We<-array(lab[i],nl[i])
			Weight<-c(Weight,We)
			}
		Weight<-Weight[-1]
		Weight<-factor(Weight)
		Option<-factor(TableSort[,2])
		for(i in 1:NR){
		if(is.na(X[i,Crit])){Option<-Option[-1]}}
			MMTable<-table(Option, Weight)
			MMT<-rbind(MMTable,index)
			MMTo<-MMT[,order(-index)]
			MMTo<-MMTo[1:itemOp[Crit],]
			MaMiTf[,,k]<-MMTo
			CoPred[k]<-sum(diag(MMTo))/NR}
#
Table5final<-array(99,dim=c(N,4,(M-N))) 
for(i in 1:SolFC) {Table5final[,,i]<-Table5[,,i]}
for(i in (SolFC+1):(M-N)){Table5final[,,i]<-Table5c[,,i]}
#
###### Cramer V #######
cramer<-array(0,dim=c(N,N))#Esta es SSj
for(i in 1:N){
  for(j in 1:N){
    cramer[i,j]<-round(assocstats(table(X[,i],X[,j]))$cramer,2)
  }
}
##################
## % Component Contamination
####################
AvgeMC<-complement$Inf_O[,N]
desde<-itemOp[Crit]
AvgeFC<-tabE3[desde:(M-N),N+1]*N/(N-1)
q<-(M-N)-SolFC
if(q==solMC){CompContam<-100*((AvgeMC-AvgeFC)/AvgeFC)
               TotContam<-100*(sum(AvgeMC)-sum(AvgeFC))/sum(AvgeFC)}
else{cat(' dsHELP: Contamination indexes only available for dim =',q,'\n\n')
CompContam<-0
TotContam<-0}
#####################
if(K>0){tipo<-"A"}
#
#OUTPUT
#		
dsFC.output<-list(
#
IniDat=X,
Call=match.call(),
CramerV=cramer,
CritItem=Crit,
Eigenval=AdjEigen,
ItONa=ItOpNa,
Match=MaMiTf,
N.Comp=(M-N),
N.Item=N,
NoOpI=itemOp,
NSs=NR,
Predict=CoPred,
SolFC=SolFC,
SolMC=solMC,
SubNa=SubNa, 
Tot.Op=sum(itemOp),
tipo=tipo,
CompContam=round(CompContam,2),
TotContam=round(TotContam,2),
dim=dim,

Inf_AC=tabE3,
Out_A=round(result[1:SolFC,],3),
ItemStat_A=Table5,
Rij_A=round(corM[1:N,1:N,1:SolFC],3),
Proj.Op_A=round(xAdj[,1:SolFC],3),
Proj.Su_A=round(yAdj[,1:SolFC],3),
Norm.Op_A=round(x[,1:SolFC],3),
Norm.Su_A=round(y[,1:SolFC],3),

Inf_B=complement$Inf_O,
Out_B=complement$Out_O,
ItemStat_B=complement$ItemStat_O,
Rij_B=complement$Rij_O[1:(N-1),1:(N-1),],
Proj.Op_B=complement$Proj.Op_O,
Proj.Su_B=complement$Proj.Su_O,
Norm.Op_B=complement$Norm.Op_O,
Norm.Su_B=complement$Norm.Su_O,

Inf_C=tabE3,# or Inf_AC
Out_C=round(resultc,3),
Norm.Op_C=round(x,3), 
Proj.Op_C=round(xAdjC,3),
Norm.Su_C=round(yG1,3),
Proj.Su_C=round(yAdjC,3)
)
class(dsFC.output)<-"ds"	
return(dsFC.output)}

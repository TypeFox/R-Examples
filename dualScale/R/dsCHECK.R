dsCHECK <-
function(X, mode='rad'){
  NR=nrow(X) #number of respondents
  N=ncol(X) #number of multipleChoice items
#
if(sum(is.na(X))>0){# there are NA values
    cat('\n------------')
    cat('\nYou had NA values in your Initial Dataset. \n Chose rad/act correction mode.')
    cat('\n------------\n')
    if(mode=="rad"){ #radical action
    cat('Beware, radical action has been taken!\n')
    XF<-na.omit(X)
    }
    if(mode=="act"){ #active inputations
      for(j in 1:N){
        if(sum(is.na(X[,j]))>0){
          maxOpt<-max(X[,j], na.rm=T)
          for(i in 1:NR){
            if(sum(is.na(X[i,j]))>0){ 
           X[i,j]<-maxOpt+1
        }}}}
    XF<-X}}
#
#########
#
if(sum(is.na(X))==0){# there aren't NA values
 X<-matrix(as.matrix(X),nrow=NR,ncol=N)
 itemOp<-c(1:N)
	for(j in 1:N)
	{itemOp[j]<-max((X[,j]),na.rm=TRUE)}
 max<-max(itemOp)
#
 XF<-X
#
 for(j in 1:N){
	Xt<-tabulate(XF[,j])
	cero<-Xt==0 
	control<-sum(as.numeric(cero))
	#
	if(control>0){ #if an option is missed
		pos<-c(1:NR)
		Xpos<-cbind(XF,pos)
		Xord<-Xpos[order(Xpos[,j]),]
		Table<-table(Xord[,j],useNA="ifany")
		Up2<-length(Table) #included NA
		TabNA<-table(Xord[,j])
		UpNA<-length(TabNA) #Not included NA
		ifelse((Up2>UpNA),NewOpt<-c(1:UpNA,NA),NewOpt<-c(1:Up2))
		NewCode<-99
	for(k in 1:Up2){
		NewCode<-c(NewCode,array(NewOpt[k],Table[k]))}
		NewCode<-c(NewCode[-1])
		Xord[,j]<-NewCode
		XF<-Xord[order(Xord[,(N+1)]),]
		XF<-XF[,-(N+1)]}}}
#cat('\n------------')
#cat('\n Now, you can use the transformed data in $TData as input for dsMC or dsFC analysis.\n')
#cat('------------\n')
dsCHECK.output<-list(InitialData=X,TData=as.data.frame(XF)) 
return(dsCHECK.output)}

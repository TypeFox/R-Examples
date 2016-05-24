bovar <-
function(y,Xmat,tr,gm,dmatsg,dmatsel,parents,parvec,w,nsplit,crit=crit)  {
   #selects best observed splitting candidate
     #create matrix to keep the highest value of the criterion
   critmax<-matrix(0,nrow=dim(Xmat)[2],ncol=6)
    critmax[,1]<-1:ncol(Xmat)
     	#only admissible rows of the designmatrix for this node
 	 dmats<-dmatsg[dmatsel==1,]
 	 index<-1:nrow(dmatsg)
 	 #split only if there are more than 1 persons in the node
 	  if(sum(gm)>1){
    critmax[,c(2:6)]<-t(sapply(1:dim(Xmat)[2],function(kk,y,Xmat,tr,gm,dmats,parents,parvec,w,nsplit,crit){
     bos(y,Xmat[,kk],tr,gm,dmats,parents,parvec,w,nsplit,crit=crit) },
     Xmat=Xmat,gm=gm,y=y,tr=tr,parents=parents,nsplit=nsplit,w=w,parvec=parvec,dmats=dmats,crit=crit) )
            #which predictor is the best spllting candidate for this split?
    bestrow<-which(critmax[,4]==max(critmax[,4]))[1] }
    else{bestrow<-1}
    #print(critmax[bestrow,])
           #rowstarnew is rownumber in the general designmatrix dmatsg
     if( critmax[bestrow,3]!=0){
     critmax[bestrow,3]<-index[dmatsel==1][critmax[bestrow,3]] }
     #print(critmax[bestrow,])
      return(critmax[bestrow,])}

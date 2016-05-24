bos <-
function(y,x,tr,gm,dmats,parents,parvec,w,nsplit,crit="es",plotc=FALSE){
  #computes best observed split for a splitting variable candidate
  #ouput is result (rowvector of length 5) with: optimal split point, rownumber of dmats, max value of C, Cdif, and Ccard
 	#y = outcome variable
 	#x = predictor to be used for the split
 	#tr = treatment variable(1=treatment A; 2=treatment B)
 	#gm = indicator vector of persons in the particular parent node that is split
 	#parents contains information of rest of the parent nodes (cardinalities t0, cardinalities t1, and sum y|t=0, sum y|t=1)
 	#dmats=designmatrix with admissible assignments 
 	#parvec=parameter a1 and a2; 
 	#w=vector with weights w1,and w2.
 	#nsplit: number of split
 	#if plotc is T, then the value of the criterium is plotted for all splitpoints and for all possible partitions
  z<-unique(sort(x[gm==1]))
 	n<-length(z)

  #n is  the number of distinct values (split points) of predictor x
  crittot<-matrix(0,nrow=dim(dmats)[1],ncol=n)
	crit1<-matrix(0,nrow=dim(dmats)[1],ncol=n)
	crit2<-matrix(0,nrow=dim(dmats)[1],ncol=n)
  #crittot is a matrix with number of rows is number of possible assignments (rows of Ds)
 	#and number of columns is n = number of splitpoints; the cells contain the value of the criterion C
  #Perform for each possible split point vik:
  if(n>1){
	      for  (i in 1:(n-1)){
	      #print(i)
        	#create indicator matrix of child nodes (Gch) after split on z[i]:
     	splitpoint<-(z[i]+z[i+1]) /2
        Gch<- makeGchmat(gm,x,splitpoint)
        #child will be filled with information of the two childnodes (cardinalities t0, cardinalities t1, and sum y|t=0, sum y|t=1)
        child<-cpmat(Gch,y,tr,crit=crit)
        if(nsplit==1){End<-child}
        if(nsplit>1){
        End<-as.matrix(rbind(parents,child) )}
           		dimnames(End)<-NULL
    		#print(End)
        #End is matrix E that contains the information of all the end nodes after a split (cardinalities t0, cardinalities t1, and sum y|t=0, sum y|t=1)
        ##check conditions for each row of dmats and if these are met compute the value of C which will be collected in matrix crittot
     		check1<-ctc(End,parvec)
     		 #cat("treatment cardinality cond", check1, "\n")
     		 if(check1==1){
     		 check2<-cmd(End,dmats) 
         #cat("mean difference cond", check2, "\n")
         if( sum(check2)!=0) {
           dmats2<-dmats[check2==1,]
          if( sum(check2)==1) {    dmats2<-t(dmats2)}
          cdat<-computeC(End,dmats2,w)
          crittot[ check2==1,i]<-cdat$crittot
          crit1[ check2==1,i]<-cdat$critdif
          crit2[ check2==1,i]<-cdat$critcard }
          }
          } 
 #print(crittot)
 #cat("crittot is",crittot,"\n")
  ccmax<-apply(crittot,2,max)
  
 	#ccmax is vector with maximum value of C for each split point
	colstar<-which(ccmax==max(crittot))[1]
  #colstar is the columnnumber referring to the splitpoint on variable Xk that results in the highest value of C
   rcmax<-apply(crittot,1,max)
	#rcmax is vector with maximum value of C for each row of the design matrix Ds
	rowstar<- which(rcmax==max(crittot))[1]
	if(plotc==T){
	vec<-1:dim(dmats)[1]
	par(mfrow=c(1,2))
	plot(z[ccmax!=0],ccmax[ccmax!=0])
  plot(vec,rcmax) }
  #rowstar is the rownumber referring to the row of dmats (Ds) that results in the highest value of C for this particular predictor variable
     splitpoint<- (z[colstar]+z[colstar+1] )/2
   result<-c(splitpoint,rowstar,crittot[rowstar,colstar],crit1[rowstar,colstar],crit2[rowstar,colstar])}
   if (n==1){result<-numeric(5)}
   return(result)
	}

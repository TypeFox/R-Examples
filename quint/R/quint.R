quint <-
function(formula,data,control=NULL)   {
    #Dataformat without use of formula:
    #dat:data; first column in dataframe = the response variable
    #second column in dataframe = the dichotomous treatment vector (coded with treatment A=1 and treatment B=2)
    #rest of the columns in dataframe are the predictors
    
    tcall<-match.call()
    #Specify data:  	
    dat<-as.data.frame(data);
    if(missing(formula)){
      y<-dat[,1];     tr<-dat[,2];      Xmat<-dat[,-c(1,2)] 
      dat<-na.omit(dat); #omit rows with missing values 
       if(length(levels(as.factor(tr)))!=2){
        stop("Quint cannot be performed. The number of treatment conditions does not equal 2.")
      }}
    else{
      F1 <- Formula(formula)
      mf1 <- model.frame(F1, data = dat)
      y<-as.matrix(mf1[,1])
      origtr<-as.factor(mf1[,2]) 
      tr<-as.numeric(origtr)
      if(length(levels(origtr))!=2){
        stop("Quint cannot be performed. The number of treatment conditions does not equal 2.")
      }
       Xmat<-mf1[,3:dim(mf1)[2]] 
       dat<-cbind(y,tr,Xmat) ;
      dat<-na.omit(dat); #omit rows with missing values
     
      cat("Treatment variable (T) equals 1 corresponds to",attr(F1,"rhs")[[1]], "=", levels(origtr)[1],"\n")
      cat("Treatment variable (T) equals 2 corresponds to",attr(F1,"rhs")[[1]], "=", levels(origtr)[2],"\n") 
      names(dat)[1:2]<-names(mf1)[1:2] 
     }
    cat("The sample size in the analysis is",dim(dat)[1],"\n")
    #Check if the background characteristics are not categorical
    classx<-sapply(1:dim(Xmat)[2],function(kk,Xmat){class(Xmat[,kk])},Xmat=Xmat)
    if(sum(classx=="factor")!=0){
      stop("Quint cannot be performed. Some background characteristics are not numeric.")
    }
       N<-length(y)
       if(is.null(control)) { 	
      #Use default control parameters and criterion
    	control <- quint.control()
      	} 
      	#specify criterion , parameters a and b  (parvec)  and  weights:
    crit<-control$crit    ;     parvec<-control$parvec ;    w<-control$w
    #specify Lmax
    maxl<-control$maxl
    #if no control argument was specified ,use default parameter values
     #Default parameters a1 and a2 for treatment cardinality condition:
 
    if(is.null(parvec)){
    a1<- round(sum(tr==1)/10);a2<- round(sum(tr==2)/10)
    parvec<-c(a1,a2)
    control$parvec<-parvec}
    #print(parvec)
    
    if(is.null(w)){
    #edif=expected mean difference between treatment and control; default value for effect size criterion: edif = 3 (=Cohen's d), 
  #and for difference in means criterion: edif= IQR(Y)  
     edif<-ifelse(crit=="es",3,IQR(y))
  w1<- 1/log(1+edif)
  w2<-1/log(length(y)/2) 
  w<-c(w1,w2)
  control$w<-w}
    
               
     
     ##Create matrix for results
     allresults<-matrix(0,nrow=maxl-1,ncol=6); 
    
    ##Start of the tree growing: all persons are in the rootnode.   L=1  ; Criterion value (cmax)=0
    root<-rep(1,length(y))
    cmax<-0
 
 
    #Step 1
    #Generate design matrix D with admissable assignments after first split
    dmat1<-matrix(c(1,2,2,1),nrow=2)
      
    
     #Select the optimal triplet for the first split: the triplet resulting in the maximum value of the criterion (critmax1) 
      #use the rootnode information: cardinality t=1, cardinality t=2, meant1 ,meant0
     rootvec<-c(sum(tr==1),sum(tr==2),mean(y[tr==1])-mean(y[tr==2]))    ;   
     
     critmax1<-bovar(y,Xmat,tr,root,dmat1,rep(1,nrow(dmat1)),rootvec,parvec,w,1,crit=crit)
     #critmax1: splitting variable, splitpoint, row of dmats, total value of C, values of Difference in treatment outcome component and cardinality component (critdifmean,critcard)
     
     #Make the first split
     Gmat<-makeGchmat(root,Xmat[,critmax1[1]],critmax1[2])
    cat("split 1","\n") ;   cat("#leaves is 2","\n")
    
    ##Keep the child node numbers  nnum; #ncol(Gmat) is current number of leaves (=number of candidate parentnodes)=L ;    #ncol(Gmat)+1 is total number of leaves after the split  (Lafter)
    nnum<-c(2,3)  ;    L<-ncol(Gmat)
    
     ##Keep the results (split information, fit information, end node information) after the first split    
     if(critmax1[4]!=0){
     allresults[1,]<-c(1,critmax1[-3])  
     dmatrow<-dmat1[critmax1[3],]
     cmax<-allresults[1,4]
     endinf<-ctmat(Gmat,y,tr)
     }  
     else{ ##if there is no optimal triplet for the first split:
     Gmat<-Gmat*0; 
     dmatrow<-c(0,0)
     endinf<-matrix(0,ncol=8,nrow=2)}
     
     ##Check the qualitative interaction condition:  Cohen's d in the leafs after the first split >=dmin
     qualint<-"Present" 
    if(abs(endinf[1,7])<control$dmin | abs(endinf[2,7])<control$dmin  ) {L<-maxl
    stop("The qualitative interaction condition is not satisfied: One or both of the effect sizes are lower than the absolute value of ",control$dmin,". There is no qualitative interaction present in the data.","\n") 
     
                                                                         }
    
      ##Perform bias-corrected bootstrapping for the first split:             
 
      if(control$Boot==TRUE&cmax!=0){
      #initiate bootstrap with stratification on treatment groups:
      indexboot<-Bootstrap(y,control$B,tr)
      critmax1boot<-matrix(0,ncol=6,nrow=control$B)
  
      #initialize matrices to keep results
      Gmattrain<-array(0,dim=c(N,maxl,control$B))
      Gmattest<-array(0,dim=c(N,maxl,control$B))
       allresultsboot<-array(0,dim=c(maxl-1,9,control$B))
      #find best first split for the k training sets
        for  (b in 1:control$B ){
         cat("Bootstrap sample ",b,"\n")
         ##use the bootstrap data as training set
        critmax1boot[b,]<-bovar(y[indexboot[,b]],Xmat[indexboot[,b],],tr[indexboot[,b]],root,dmat1,rep(1,nrow(dmat1)),rootvec,parvec,w,1,crit=crit)
           Gmattrain[,c(1:2),b]<-makeGchmat(root,Xmat[indexboot[,b],critmax1boot[b,1]],critmax1boot[b,2])
        ##use the original data as testset
        Gmattest[,c(1:2),b]<-makeGchmat(root,Xmat[,critmax1boot[b,1]],critmax1boot[b,2]) 
         End<-cpmat(Gmattest[,c(1:2),b],y,tr,crit=crit)
        #select the right row in the design matrix
        dmatsel<-t(dmat1[critmax1boot[b,3],])
    
       allresultsboot[1,c(1:8),b] <- c(1,critmax1boot[b,c(1:2)],computeCtest(End,dmatsel,w))
        allresultsboot[1,9,b]<-critmax1boot[b,4]-allresultsboot[1,4,b]
        if(critmax1boot[b,4]==0){allresultsboot[1,,b]<-NA}
      } 
    }
  
      
    #Repeat the tree growing procedure
    stopc<-0
    
    while(L<maxl){
        cat("current value of C", cmax,"\n")
                cat("split",L,"\n")
        Lafter<-ncol(Gmat)+1
        cat("#leaves is",Lafter,"\n")
        ##make a designmatrix (dmat) for the admissible assignments of the leaves after the split
        dmat<-makedmat(Lafter)
        dmatsg<-makedmats(dmat)
        #make parentnode information matrix, select best observed parent node (with optimal triplet)
        parent<-cpmat(Gmat,y,tr,crit=crit)
        critmax<-bonode(Gmat,y,Xmat,tr,dmatrow,dmatsg,parent,parvec,w,L,crit=crit)
         ##Perform the best split and keep results
         Gmatch<-makeGchmat(Gmat[,critmax[1]],Xmat[,critmax[2]],critmax[3])
        Gmatnew<-cbind(Gmat[,-critmax[1]],Gmatch) 
        allresults[L,]<-c(nnum[critmax[1]],critmax[2:3],critmax[5:7]) 
        dmatrownew<-dmatsg[critmax[4],]
        
        #check if cmax new is higher than current value
        
          if(allresults[L,4]<=cmax){
        cat("splitting process stopped after number of leaves equals",L,"because new value of C was not higher than current value of C","\n")
        stopc<-1}
        
        ##repeat this procedure for the bootstrap samples
      if(control$Boot==TRUE&stopc!=1){
        critmaxboot<-matrix(0,nrow=control$B,ncol=7)
     for (b in 1:control$B){
          cat("Bootstrap sample ",b,"\n")
            #make parentnode information matrix pmat
        parent<-cpmat(Gmattrain[,c(1:(Lafter-1)),b],y[indexboot[,b]],tr[indexboot[,b]],crit=crit)
   
            critmaxboot[b,]<- bonode(Gmattrain[,c(1:(Lafter-1)),b],y[indexboot[,b]],Xmat[indexboot[,b],],tr[indexboot[,b]],dmatrow,dmatsg,parent,parvec,w,L,crit=crit)
            #best predictor and node of this split for the training samples
             Gmattrainch<-makeGchmat(Gmattrain[,critmaxboot[b,1],b],Xmat[indexboot[,b],critmaxboot[b,2]],critmaxboot[b,3])
    Gmattrain[,c(1:Lafter),b]<-cbind(Gmattrain[,c(1:(Lafter-1))[-critmaxboot[b,1]],b],Gmattrainch)
    Gmattestch<-makeGchmat(Gmattest[,critmaxboot[b,1],b],Xmat[,critmaxboot[b,2]],critmaxboot[b,3])
    Gmattest[,c(1:Lafter),b]<-cbind(Gmattest[,c(1:(Lafter-1))[-critmaxboot[b,1]],b],Gmattestch)
    ##compute criterion value for the test sets
     End<-cpmat(Gmattest[,c(1:Lafter),b],y,tr,crit=crit)
     #select the right row in the design matrix
     if(critmaxboot[b,5]!=0){
     dmatsel<-t(dmatsg[critmaxboot[b,4],]) 
     allresultsboot[L,c(1:8),b] <- c(nnum[critmaxboot[b,1]],critmaxboot[b,2],critmaxboot[b,3],computeCtest(End,dmatsel,w))
        allresultsboot[L,9,b]<-critmaxboot[b,5]-allresultsboot[L,4,b] }
      if(critmaxboot[b,5]==0){ 
      allresultsboot[L,,b] <-NA
      }  
     
    } 
     if(sum(is.na(allresultsboot[L,9,]))/control$B > .10 ){
    warning("After split ",L,", the partitioning criterion cannot be computed in more than 10 percent of the bootstrap samples. The split is unstable." )
     }
    }
  
         
      #update the parameters after the split:
   
     if(stopc==0) {
     Gmat<-Gmatnew
     dmatrow<-dmatrownew
     cmax<-allresults[L,4] 
     L<-ncol(Gmat)
     nnum<-c(nnum[-critmax[1]],nnum[critmax[1]]*2,nnum[critmax[1]]*2+1)
           }  
    else{
      L<-maxl}

        #end of while loop
    }
    
     
     Lfinal<-ncol(Gmat)  #Lfinal=final number of leaves of the tree
   
    #create  endnode information of the tree
    endinf<-matrix(0,nrow=length(nnum),ncol=10)
    if(cmax!=0){
    endinf[,c(2:9)]<-ctmat(Gmat,y,tr)}
    endinf<-data.frame(endinf)
    endinf[,10]<-dmatrow
    endinf[,1]<-nnum 
    index<-leafnum(nnum)
    endinf<-endinf[index,]
    rownames(endinf)<-paste("Leaf ",1:Lfinal,sep="")
    colnames(endinf)<-c("node","#(T=1)","meanY|T=1","SD|T=1","#(T=2)","meanY|T=2","SD|T=2","d","se","class")
    
    if(Lfinal==2){allresults<-c(2,allresults[1,])}
    if(Lfinal>2){
    allresults<-cbind(2:Lfinal,allresults[1:(Lfinal-1),])}
    
    #compute final estimate of optimism and standard error:
    if(control$Boot==TRUE){
       #raw mean and sd:
    opt<-sapply(1:(Lfinal-1),function(kk,allresultsboot){mean(allresultsboot[kk,9,],na.rm=T)},allresultsboot=allresultsboot)
     se_opt<-sapply(1:(Lfinal-1),function(kk,allresultsboot){sd(allresultsboot[kk,9,],na.rm=T)/sqrt(sum(!is.na(allresultsboot[kk,9,])))},allresultsboot=allresultsboot)
    if(Lfinal==2){allresults<-c(allresults[1:5],allresults[5]-opt,opt,se_opt)
     allresults<-data.frame( t(allresults))}
    if(Lfinal>2){
     allresults<-cbind(allresults[,1:5],allresults[,5]-opt,opt,se_opt)
       allresults<-data.frame( allresults)}
       allresults[,3]<-colnames(Xmat)[allresults[,3]]
     splitnr<- 1:(Lfinal-1)
    allresults<-cbind(splitnr,allresults)
     colnames(allresults)<-c("split","#leaves","parentnode","splittingvar","splitpoint","apparent","biascorrected","opt","se")
     }
       
    if(control$Boot==FALSE){
    if(Lfinal>2){
      allresults<-data.frame( allresults) }
      if(Lfinal==2){allresults<-data.frame( t(allresults))}
    allresults[,3]<-colnames(Xmat)[allresults[,3]]
     splitnr<- 1:(Lfinal-1)
    allresults<-cbind(splitnr,allresults)
     colnames(allresults)<-c("split","#leaves","parentnode","splittingvar","splitpoint","apparent","difcomponent","cardcomponent")
     }
     colnames(Gmat)<-nnum
     
     ##split information (si): also include childnode numbers
      si<-allresults[,3:5]
      cn<-paste(si[,1]*2,si[,1]*2+1,sep="," )
      si<-cbind(parentnode=si[,1],childnodes=cn,si[,2:3])
      rownames(si)<-paste("Split ",1:(Lfinal-1),sep="")
    if(control$Boot==FALSE){
      object<-list(call=tcall,crit=crit,control=control,data=dat,si=si,fi=allresults[,c(1:2,6:8)],li=endinf, nind=Gmat[,index])}
     if(control$Boot==TRUE){ 
       nam<-c("parentnode","splittingvar","splitpoint",
                               "C_boot","C_compdif","checkdif","C_compcard","checkcard","opt")
       dimnames(allresultsboot)<-list(NULL,nam,NULL)
       object<-list(call=tcall,crit=crit,control=control,indexboot=indexboot,data=dat,si=si,fi=allresults[,c(1:2,6:9)],li=endinf, nind=Gmat[,index],siboot=allresultsboot)} 
   
   # object$data=dat  ## cor ## data toegevoegen, dit kan beter
    class(object)<-"quint"
    return(object)
    }

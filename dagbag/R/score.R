score<-function(Y, n.boot=0, score.type="BIC", threshold=0, max.step=500,  ini.adj.matrix=NULL, blacklist=NULL, whitelist=NULL, standardize=TRUE,  standardize.boot=TRUE, random.forest=FALSE, random.step.length=NULL, nrestart=0, perturb=0, shuffle=FALSE, print=FALSE, EPS=1e-06){

##when n.boot==0, other results, such as movement, delta.min, and final.step are returned.

 n=nrow(Y)                ## sample size 
 p=ncol(Y)                ## number of variables

 if(score.type=="BIC"){   ## BIC score 
   score.type=1
 }
 if(score.type=="likelihood"){ ## likelihood score 
   score.type=2
 }

 if(is.null(blacklist)){
   blacklist=matrix(0,p,p)
 }

 if(is.null(whitelist)){
   whitelist=matrix(0,p,p)
 }

 if(shuffle==TRUE){        ##shuffle nodes order
   order.s=sample(1:p,p,replace=FALSE)
   Y=Y[,order.s]
   blacklist=blacklist[order.s,order.s]
   whitelist=whitelist[order.s,order.s]
 }

 if(n.boot>0){            ## bootstrap 
   adj.matrix=array(0,c(p,p,n.boot))
 }
 else{
   adj.matrix=matrix(0,p,p)
 }

 b.ini.adj=0
 if(is.null(ini.adj.matrix)){ ## initial graph 
   ini.adj.matrix=matrix(0,p,p)
 }
 else{
   b.ini.adj=1
 }

 if(random.forest==FALSE){##if random.forest is not used, set random.step.length as arbitrary.
   random.step.length=rep(1e+04,2)
 }  ##end if random.froest 

 final.step=0       ## record the index of the final step 
 bic.ok=numeric(max.step)-99     ## recoed whether there is error of score calculation at each step
 movement=matrix(0,max.step,3)   ## record selected operation (move) at each search step
 delta.min=numeric(max.step)-99  ## record the minimum change of the score at each step
 BICscore=numeric(p)        ## record the current BIC score for each neighborhood 



#### 

  #dyn.load("BIC_comp_dposv.so")
  junk<-.C("score_interf", 
           as.integer(n), 
           as.integer( p), 
           as.double(Y), 
           as.integer(score.type),
           as.integer(standardize),
           as.double(threshold), 
           as.double(EPS), 
           as.integer(max.step),
           as.integer(print), 
           as.integer(ini.adj.matrix),
           as.integer(b.ini.adj),
           as.integer(blacklist),
           as.integer(whitelist),
           as.integer(n.boot),
           as.integer(standardize.boot),
           as.integer(random.forest),
           as.double(random.step.length), 
           as.integer(nrestart),
           as.integer(perturb),
           movement=as.integer(movement), 
           adj.matrix=as.integer(adj.matrix), 
           delta.min=as.double(delta.min), 
           bic.ok=as.integer(bic.ok),
           final.step=as.integer(final.step)
         )


#############
####outputing results 
 temp<-new.env()

 if(n.boot>1){  ##with bootstrap
     m.temp=array(junk$adj.matrix,c(p,p,n.boot))

     if(shuffle==TRUE){     ### shuffle back nodes order 
       res.m=array(0,c(p,p,n.boot))
       res.m[order.s,order.s,]=m.temp
       temp$adj.matrix=res.m
       #temp$shuffle.order=order.s
     }else{
       temp$adj.matrix=m.temp      ## for bootstrap, only return the adjacency matrices (n.boot in total)
     }

 }else{          ## no bootstrap
     m.temp=matrix(junk$adj.matrix,p,p)
     if(shuffle==TRUE){ ### shuffle back nodes order 
       res.m=matrix(0,p,p)
       res.m[order.s,order.s]=m.temp
       temp$adj.matrix=res.m
       #temp$shuffle.order=order.s
       shuffle.movement=matrix(junk$movement,max.step,3) 
       shuffle.movement=shuffle.movement[1:(junk$final.step-1),]  
       new.movement=matrix(0,(junk$final.step-1),3)
       new.movement[,3]=shuffle.movement[,3]
       new.movement[,1]=order.s[shuffle.movement[1:(junk$final.step-1),1]]
       new.movement[,2]=order.s[shuffle.movement[1:(junk$final.step-1),2]]
       temp$movement=new.movement
     }else{
       temp$adj.matrix=m.temp
       temp$movement=matrix(junk$movement,max.step,3)
     }

     temp$final.step=junk$final.step            ## return adjacency matrix, movement, etc. 
     temp$bic.ok=junk$bic.ok
     temp$delta.min=junk$delta.min
 }

res=as.list(temp)
return(res)

}


############## auxiliary functions
#####
y.stand<-function(Y){
###standardize data matrix Y (n by p)   to have mean zero and sd 1
  n=nrow(Y)   ## sample size
  p=ncol(Y)   ## number of variables 
  Y.gm=apply(Y, 2, mean)   
  Y.gsd=apply(Y, 2, sd)   
  Y.new=(Y-matrix(Y.gm, n, p, byrow=T))/matrix(Y.gsd, n, p, byrow=T) ##standardization to have mean 0 and sd 1  
  return(Y.n=Y.new)
}

###
result_skeleton<-function(adj.m, true.ske){
## Find the total number of skeleton edges and the number of correct skeleton edges from the adjacency matrix  of an estimated DAG 
## compared with the true.ske. Useful for summarize simulation results. 
##parameters: adj.m: adjacency matrix (a 0-1 matrix), true.ske -- true skeleton ( a symmetric 0-1 matrix)

  diag(adj.m)=0
  tt=adj.m+t(adj.m) 
  correct.c=sum((tt>0)&(true.ske>0))/2    
  total.c=sum(tt>0)/2
  return(c(total.c, correct.c))
}

###
#######
move_adj<-function(p,movement){
###generate the adjacency matrix  of the final learned DAG from the recorded movements of the search aglorithm.
   adj.result=NULL
   adj.matrix=matrix(0,p,p)
   step=nrow(movement)
   for(i in 1:step){
     if(movement[i,3]>0){
       adj.matrix[movement[i,1],movement[i,2]]=(2-movement[i,3])>0   ##updated by jie on 1/10/2012 for efficiency 
       adj.matrix[movement[i,2],movement[i,1]]=(2-movement[i,3])<0
     }
     adj.result[[i]]=adj.matrix
   }
   return(adj.result)
}

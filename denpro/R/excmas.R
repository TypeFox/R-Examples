excmas<-function(lst){
#
parents<-lst$parent
volumes<-lst$volume
levels<-lst$level
#
nodelkm<-length(parents)
excmasses<-matrix(1,nodelkm,1)
#
mlkm<-moodilkm(parents)
modloc<-mlkm$modloc    #pointers to modes
lkm<-mlkm$lkm       
#
added<-matrix(0,nodelkm,1)  #1 if we have visited this node
#
for (i in 1:lkm){
    node<-modloc[i]
    # calculate curexc
    par<-parents[node]
    if (par==0) valpar<-0 else valpar<-levels[par] 
    curexc<-(levels[node]-valpar)*volumes[node]
    #
    excmasses[node]<-curexc
    while (parents[node]>0){
         node<-parents[node]
         if (added[node]==0){   
           # calculate curexc
           par<-parents[node]
           if (par==0) valpar<-0 else valpar<-levels[par] 
           curexc<-curexc+(levels[node]-valpar)*volumes[node] 
           #
           excmasses[node]<-curexc 
           added[node]<-1 
         }
         else{   #add only previous cumulative 
            excmasses[node]<-excmasses[node]+curexc
         }
    }
}
return(t(excmasses))
}





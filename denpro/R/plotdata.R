plotdata<-function(roots,child,sibling,sibord,levels,volumes,vecs)
{
#plots level-set profile

#parents<-c(0,1,1,0,4,2)
#levels<-c(1,2,2,1,2,3)
#volumes<-c(4,2,1,2,1,1)

itemnum<-length(volumes)

#vecs<-matrix(NA,itemnum,4)
#vecs<-alloroot(vecs,roots,sibord,levels,volumes)

rootnum<-length(roots)
left<-child
right<-sibling

for (i in 1:rootnum){
    pino<-matrix(0,itemnum,1)
    pino[1]<-roots[i]  
    pinin<-1
    while (pinin>0){
        cur<-pino[pinin]      #take from stack
        pinin<-pinin-1
        if (left[cur]>0){     #if not leaf (root may be leaf)
           vecs<-allokoi(vecs,cur,child,sibling,sibord,levels,volumes)   
        }
        if (right[cur]>0){    #if right exists, put to stack
            pinin<-pinin+1
            pino[pinin]<-right[cur]
        }
        while (left[cur]>0){    #go to leaf and put right nodes to stack
             cur<-left[cur]
             if (left[cur]>0){  #if not leaf
                vecs<-allokoi(vecs,cur,child,sibling,sibord,levels,volumes)
             }
             if (right[cur]>0){ #if right exists, put to stack
                pinin<-pinin+1
                pino[pinin]<-right[cur]
             }
        }
    }
}       
#
return(vecs)
}








excmas.bylevel<-function(lst,levnum)
{
#source("~/denpro/R/excmas.bylevel.R")
#excmas.bylevel(lst,20)

levexc<-matrix(0,levnum,1)

maxlev<-max(lst$level)
step<-maxlev/levnum
nodelkm<-length(lst$parent)

mlkm<-moodilkm(lst$parent)
modloc<-mlkm$modloc    #pointers to modes
lkm<-mlkm$lkm       

added<-matrix(0,nodelkm,1)  #1 if we have visited this node

i<-1
while (i<=lkm){
    node<-modloc[i]
    # calculate curexc
    par<-lst$parent[node]
    if (par==0) valpar<-0 else valpar<-lst$level[par] 
    curexc<-(lst$level[node]-valpar)*lst$volume[node]
    
    nodelevind<-min(max(round(lst$level[node]/step),1),levnum)    
    levexc[1:nodelevind]<-levexc[1:nodelevind]+curexc

    while (lst$parent[node]>0){
         node<-lst$parent[node]
         if (added[node]==0){   
           # calculate curexc
           par<-lst$parent[node]
           if (par==0) valpar<-0 else valpar<-lst$level[par] 
           curexc<-(lst$level[node]-valpar)*lst$volume[node] 
           
           nodelevind<-min(max(round(lst$level[node]/step),1),levnum)    
           levexc[1:nodelevind]<-levexc[1:nodelevind]+curexc

           added[node]<-1 
         }
    }
    i<-i+1
}

levexc<-levexc/levexc[1]

diffe<-matrix(0,length(levexc),1)
for (i in 1:(length(levexc)-1)) diffe[i]<-(levexc[i+1]-levexc[i])/step
diffe[length(diffe)]<-diffe[length(diffe)-1]

return(list(levexc=levexc,diffe=diffe))
}





taillevel<-function(root,#child,sibling,
parent,volume,proba)
{
mt<-multitree(parent)
child<-mt$child
sibling<-mt$sibling

nodenum<-length(child)
level<-matrix(0,nodenum,1)
pino<-matrix(0,nodenum,1)

pino[1]<-root
pinin<-1
while (pinin>0){
      cur<-pino[pinin]      #take from stack
      pinin<-pinin-1

      chi<-child[cur]
      pare<-parent[cur]

      prochi<-0
      nexchi<-chi
      while (nexchi>0){
           prochi<-prochi+proba[nexchi]
           nexchi<-sibling[nexchi]
      }
      if (pare==0) levelpare<-0 else levelpare<-level[pare]

      level[cur]<-levelpare+(proba[cur]-prochi)/volume[cur]

      if (sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-sibling[cur]
      }
      while (child[cur]>0){    #go to left and put right nodes to stack
            cur<-child[cur]

            chi<-child[cur]
            pare<-parent[cur]

            prochi<-0
            nexchi<-chi
            while (nexchi>0){
                 prochi<-prochi+proba[nexchi]
                 nexchi<-sibling[nexchi]
            }
            if (pare==0) levelpare<-0 else levelpare<-level[pare]

            level[cur]<-levelpare+(proba[cur]-prochi)/volume[cur]

            if (sibling[cur]>0){  #if candi has siblings
                pinin<-pinin+1
                pino[pinin]<-sibling[cur]
            } 
      }
}
return(level)

}






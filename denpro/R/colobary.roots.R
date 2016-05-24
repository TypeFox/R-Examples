colobary.roots<-function(parent,level,colothre=min(level),paletti=NULL)
{
mt<-multitree(parent) #roots<-mt$roots child<-mt$child sibling<-mt$sibling

itemnum<-length(mt$child)
colot<-matrix("",itemnum,1)
rootnum<-length(mt$roots)
if (is.null(paletti))
 paletti<-c("red","blue","green",
 "orange","navy","darkgreen",
 "orchid","aquamarine","turquoise",
 "pink","violet","magenta","chocolate","cyan",
 colors()[50:657],colors()[50:657])
hep<-1

for (i in 1:rootnum){
   curroot<-mt$roots[i]
   colot[curroot]<-paletti[hep]
   hep<-hep+1
   if (mt$child[curroot]>0){
      pino<-matrix(0,itemnum,1)
      pino[1]<-mt$child[curroot]
      pinin<-1
      while (pinin>0){
          cur<-pino[pinin]      #take from stack
          pinin<-pinin-1
          if ((mt$sibling[cur]==0)||(level[parent[cur]]<colothre)) 
              colot[cur]<-colot[parent[cur]]
          else{ 
                colot[cur]<-paletti[hep]
                hep<-hep+1
          } 
          # put to the stack 
          if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
          }
          # go to left and put right nodes to the stack
          while (mt$child[cur]>0){   
             cur<-mt$child[cur]
             if ((mt$sibling[cur]==0)||(level[parent[cur]]<colothre)) 
                       colot[cur]<-colot[parent[cur]]
             else{ 
                    colot[cur]<-paletti[hep]
                    hep<-hep+1
             } 
             if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
             }
           }
       }#while (pinin>0)
   }
}                           
return(colot)
}





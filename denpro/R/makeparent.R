makeparent<-function(left,right)
{
parent<-matrix(0,length(left),1)

pino<-matrix(0,length(left),1)
pinin<-1
pino[1]<-1

while (pinin>0){

    node<-pino[pinin]
    pinin<-pinin-1

    if (left[node]>0){
       parent[left[node]]<-node
       parent[right[node]]<-node

       pinin<-pinin+1
       pino[pinin]<-right[node]
    }

    while (left[node]>0){
       
        node<-left[node]

        if (left[node]>0){
           parent[left[node]]<-node
           parent[right[node]]<-node

           pinin<-pinin+1
           pino[pinin]<-right[node]
        }
    }

}

return(t(parent))
}

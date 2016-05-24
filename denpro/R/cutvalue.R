cutvalue<-function(roots,child,sibling,level,component,
          AtomlistAtom,AtomlistNext,valnum){
#
#from the cutted multitree, form a "newvalue",
#which gives quantized values for the kernel estimate,
#in addition the values are cutted, so that one mode is 
#removed (input is cutted multitree)
#
itemnum<-length(child)
rootnum<-length(roots)
newvalue<-matrix(0,valnum,1)

for (i in 1:rootnum){
    pino<-matrix(0,itemnum,1)
    pino[1]<-roots[i]
    #    
    pinin<-1
    while (pinin>0){
        cur<-pino[pinin]      #take from stack
        pinin<-pinin-1
        #
        node<-cur
        compo<-component[node]
        ato<-compo                          #ato is pointer to "value"
        while (ato>0){
           newvalue[AtomlistAtom[ato]]<-level[node]
           ato<-AtomlistNext[ato]
        }
        #
        if (sibling[cur]>0){
              pinin<-pinin+1
              pino[pinin]<-sibling[cur]
        }
        while (child[cur]>0){    #go to left and put right nodes to stack
              cur<-child[cur]
              #
              node<-cur
              compo<-component[node]
              ato<-compo                    #ato is pointer to "value"
              while (ato>0){
                  newvalue[AtomlistAtom[ato]]<-level[node]
                  ato<-AtomlistNext[ato]
              }
              #
              if (sibling[cur]>0){  #if cur has siblings
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
             }
        }
    }
}
#
return(newvalue)
}



declevgen<-function(tobehandled,curterminalnum,
left,right,val,vec,infopointer,parent,
low,upp,
d)
# source("~/denpro/R/declevgen.R")
# dl<-declevgen(tobehandled,curterminalnum,left,right,val,vec,
# infopointer,parent,low,upp,d)

{
nodelkm<-length(left)

separy<-matrix(0,nodelkm,1)          #=infopointer?

begsSepaNext<-matrix(0,curterminalnum,1)
begsSepaBegs<-matrix(0,curterminalnum,1)    #pointers to begsSepaAtom
begsLeftBoun<-matrix(0,curterminalnum,1)
begsRighBoun<-matrix(0,curterminalnum,1)

atomsSepaNext<-matrix(0,curterminalnum,1)   #pointers to value,index
atomsSepaAtom<-matrix(0,curterminalnum,1)
atomsLBounNext<-matrix(0,curterminalnum,1)
atomsRBounNext<-matrix(0,curterminalnum,1)

leafloc<-findleafs(left,right)

lkm<-0
node<-nodelkm
while (node>=1){ 
          # root is in position 1
    if ((leafloc[node]==1) && (tobehandled[node]==1)){   #we are in leaf 

          tobehandled[parent[node]]<-1

          lkm<-lkm+1
          separy[node]<-lkm
          atomsSepaAtom[lkm]<-infopointer[node]

          begsSepaBegs[lkm]<-lkm

          #begsLeftBoun[lkm]<-lkm
          #begsRighBoun[lkm]<-lkm

          #obs: we need not change
          #begsSepaNext, atomsSepaNext, atomsLBounNext, atomsRBounNext
          #since at the beginning set consist only one member:
          #pointer is always 0, since we do not have followers

    }
    else if (tobehandled[node]==1){   #not a leaf

        tobehandled[parent[node]]<-1

        leftbeg<-left[node]
        rightbeg<-right[node]

        if ((leftbeg==0) || (separy[leftbeg]==0)){
           #if left child does not exist

           #note that since we consider subsets of the
           #terminal nodes of the original tree, it may happen
           #that leftbeg>0 but left child does not exist
          
           separy[node]<-separy[rightbeg]

           #we need that all lists contain as many members
           #left boundary is empty, but we will make it a list
           #of empty lists

           #note<-separy[node]
           #while (note>0){
           #    begsLeftBoun[note]<-0
           #    note<-begsSepaNext[note]
           #}
           # right boundary stays same as for rightbeg
        }
        else{   # eka else
            if ((rightbeg==0) || (separy[rightbeg]==0)){
            #right child does not exist
    
                   separy[node]<-separy[leftbeg]
    
                   # left boundary stays same as for leftbeg
                   # right boundary is empty
    
                   #note<-separy[node]
                   #while (note>0){
                   #    begsRighBoun[note]<-0
                   #    note<-begsSepaNext[note]
                   #}
            }
            else{   #toka else: both children exist
                    #create boundaries 
                    
                    direktiooni<-vec[node]
                    splittiini<-val[node]     
               
                    #left boundary of right child :
                    #create/check whether empty 
                   
                    Lempty<-TRUE
                    note<-separy[rightbeg]
                    while (note>0){
                        thisnoteempty<-TRUE

                        poiju<-begsSepaBegs[note]
                        prevpoiju<-poiju 
                        while (poiju>0){
                           aatto<-atomsSepaAtom[poiju]
                           if (!(splittiini<low[aatto,direktiooni])){
                               #this atom belongs to boundary
                               if (thisnoteempty==TRUE){
                                   #poiju is the 1st non-empty
                                   begsLeftBoun[note]<-poiju
                               }                
                               Lempty<-FALSE   
                               atomsLBounNext[prevpoiju]<-poiju    
                               atomsLBounNext[poiju]<-0
                               prevpoiju<-poiju

                               thisnoteempty<-FALSE
                           }
                           poiju<-atomsSepaNext[poiju]
                        }
                        if (thisnoteempty) begsLeftBoun[note]<-0
                        note<-begsSepaNext[note]
                     }

                     #right boundary of left child

                     Rempty<-TRUE
                     note<-separy[leftbeg]
                     while (note>0){
                        thisnoteempty<-TRUE
                       
                        poiju<-begsSepaBegs[note]
                        prevpoiju<-poiju 
                        while (poiju>0){
                           aatto<-atomsSepaAtom[poiju]
                           if (!(splittiini>upp[aatto,direktiooni])){
                               #this atom belongs to boundary
                               if (thisnoteempty==TRUE){
                                   #poiju is the 1st non-empty
                                   begsRighBoun[note]<-poiju
                               }                
                               Rempty<-FALSE   
                               atomsRBounNext[prevpoiju]<-poiju    
                               atomsRBounNext[poiju]<-0
                               prevpoiju<-poiju 
                           
                               thisnoteempty<-FALSE
                           }
                           poiju<-atomsSepaNext[poiju]
                        }
                        if (thisnoteempty) begsRighBoun[note]<-0
                        note<-begsSepaNext[note]
                     }

                     #check whether one of boundaries is empty

                     if (Lempty || Rempty){
                        #one of boundaries is empty
                        ############
                        #concatenating separate parts
                        ############
                        akku<-separy[leftbeg]
                        #begsRighBoun[akku]<-0 
                        # right boundaries of sets in left child are empty
                        # begsLeftBoun[akku] does not change
                        
                        #while (begsSepaNext[akku]>0){
                        #    akku<-begsSepaNext[akku]
                        #    begsRighBoun[akku]<-0
                        #}
                        begsSepaNext[akku]<-separy[rightbeg] 
                        #concatenate list of separate sets

                        separy[node]<-separy[leftbeg]
                        akku<-separy[rightbeg]
                          
                        #begsLeftBoun[akku]<-0 
                        #left boundaries of sets in right child are empty

                        #while (begsSepaNext[akku]>0){
                        #   akku<-begsSepaNext[akku]
                        #   begsLeftBoun[akku]<-0
                        #}
                        ############# 
                        #end of concatenating
                        #############
                     }
                     else{  #both children exist, both boundaries non-empty
                     direction<-vec[node]
                     jc<-joincongen(leftbeg,rightbeg,separy,
                     begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
                     atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
                     direction,low,upp)   
                      
                     separy<-jc$separy
                     separy[node]<-jc$totbegSepary

                     begsSepaNext<-jc$begsSepaNext
                     begsSepaBegs<-jc$begsSepaBegs
                     #begsLeftBoun<-jc$begsLeftBoun
                     #begsRighBoun<-jc$begsRighBoun

                     atomsSepaNext<-jc$atomsSepaNext
                     atomsSepaAtom<-jc$atomsSepaAtom
                     #atomsLBounNext<-jc$atomsLBounNext
                     #atomsRBounNext<-jc$atomsRBounNext
                     }
            } #toka else
        } #eka else
    }  #else not a leaf
    node<-node-1
}

return(list(
totbegSepary=separy[1],
begsSepaNext=begsSepaNext,
begsSepaBegs=begsSepaBegs,
atomsSepaNext=atomsSepaNext,
atomsSepaAtom=atomsSepaAtom
))
}














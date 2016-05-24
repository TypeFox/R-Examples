joingene<-function(node,leftbeg,rightbeg,separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
direction,index){
#
            if ((leftbeg==0) || (separy[leftbeg]==0)){
#if left child does not exist    
#note that since we consider subsets of the
#terminal nodes of the original tree, it may happen
#that leftbeg>0 but left child does not exist
                separy[node]<-separy[rightbeg]
                #we need that all lists contain as many members
                #left boundary is empty, but we will make it a list
                #of empty lists
                note<-separy[node]
                while (note>0){
                      begsLeftBoun[note]<-0
                      note<-begsSepaNext[note]
                }
                # right boundary stays same as for rightbeg
            }
            else{   # eka else
                if ((rightbeg==0) || (separy[rightbeg]==0)){  
                              #right child does not exist
                   separy[node]<-separy[leftbeg]
                   # left boundary stays same as for leftbeg
                   # right boundary is empty
                   note<-separy[node]
                   while (note>0){
                       begsRighBoun[note]<-0
                       note<-begsSepaNext[note]
                   }
                } 
                else{   #toka else: both children exist
                    #check whether left boundary of right child is empty
                    Lempty<-TRUE
                    note<-separy[rightbeg]
                    while (note>0){
                        if (begsLeftBoun[note]>0){
                              Lempty<-FALSE
                        }
                        note<-begsSepaNext[note]
                     }
                     #check whether right bound of left child is empty     
                     Rempty<-TRUE
                     note<-separy[leftbeg]
                     while (note>0){
                          if (begsRighBoun[note]>0){
                                 Rempty<-FALSE
                          }
                          note<-begsSepaNext[note]
                     }
                     #check whether one of boundaries is empty
                     if (Lempty || Rempty){
                             #one of boundaries is empty
############
#concatenating separate parts
#and updating boundaries for the separate parts
#separy[node]<- concatenate separy[leftbeg],separy[rightbeg]
###########
akku<-separy[leftbeg]
begsRighBoun[akku]<-0 #right boundaries of sets in left child are empty
                      # begsLeftBoun[akku] does not change
while (begsSepaNext[akku]>0){
  akku<-begsSepaNext[akku]
  begsRighBoun[akku]<-0
}                           
begsSepaNext[akku]<-separy[rightbeg] #concatenate list of separate sets
separy[node]<-separy[leftbeg]
akku<-separy[rightbeg]
begsLeftBoun[akku]<-0 #left boundaries of sets in right child are empty
while (begsSepaNext[akku]>0){
  akku<-begsSepaNext[akku]
  begsLeftBoun[akku]<-0
}        
####################
#end of concatenating
###################
                    }
                    else{  #both children exist, both boundaries non-empty  
jc<-joinconne(leftbeg,rightbeg,separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
direction,index)   #direction<-i
#
separy<-jc$separy
separy[node]<-jc$totbegSepary 
#
begsSepaNext<-jc$begsSepaNext
begsSepaBegs<-jc$begsSepaBegs
begsLeftBoun<-jc$begsLeftBoun
begsRighBoun<-jc$begsRighBoun
#
atomsSepaNext<-jc$atomsSepaNext
atomsSepaAtom<-jc$atomsSepaAtom
atomsLBounNext<-jc$atomsLBounNext
atomsRBounNext<-jc$atomsRBounNext
                        #
                    }
                } #toka else
            } #eka else
###########################################
#          end of child joining
###########################################
return(list(separy=separy,
begsSepaNext=begsSepaNext,begsSepaBegs=begsSepaBegs,
begsLeftBoun=begsLeftBoun,begsRighBoun=begsRighBoun,
atomsSepaNext=atomsSepaNext,atomsSepaAtom=atomsSepaAtom,
atomsLBounNext=atomsLBounNext,atomsRBounNext=atomsRBounNext))
}











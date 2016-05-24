declevdya<-function(beg,AtomlistNext,AtomlistAtom,kg,N,nodenumOfDyaker,
terminalnum,d){
#
#beg is pointer to AtomlistAtom
#nodenumOfDyaker is the num of nodes of the _original dyaker_
#terminalnum is the num of terminal nodes of the "current dyaker"
#
#kg=kernel estimate is represented as binary tree with "nodenumOfDyaker" nodes:
#   vectors whose length is "nodenum":
#     -left,right,parent
#     -infopointer: pointer to "value" and "index" (only for terminal nodes)
#   additional data structures  
#     -value, index
#
#to be created:
#-separy, vector with length "nodenumOfDyaker", points to begsSepaBegs
#-begsSepaBegs, begsSepaNext, begsLeftBoun, begsRighBoun
#   vectors of same length as "value" and "index",
#   list of separate sets, index gives starting point for set in atomsSepaNext
#-atomsSepaAtom, atomsSepaNext, atomsLBounNext, atomsRBounNext
#   vector of same length as "value" and "index",
#   list of atoms in separate sets, index gives the atom in value and index
#
#return: 
#begs: list of beginnings of lists
#AtomlistAtoms, AtomlistNext: list of lists of atoms
#
left<-kg$left
right<-kg$right
parent<-kg$parent
index<-kg$index
nodefinder<-kg$nodefinder
#infopointer<-kg$infopointer
#
separy<-matrix(0,nodenumOfDyaker,1)
#
begsSepaNext<-matrix(0,terminalnum,1)
begsSepaBegs<-matrix(0,terminalnum,1)   #pointers to begsSepaAtoms
begsLeftBoun<-matrix(0,terminalnum,1)
begsRighBoun<-matrix(0,terminalnum,1)
#
atomsSepaNext<-matrix(0,terminalnum,1)  #pointers to value,index
atomsSepaAtom<-matrix(0,terminalnum,1)
atomsLBounNext<-matrix(0,terminalnum,1)
atomsRBounNext<-matrix(0,terminalnum,1)
#
nextFloor<-matrix(0,terminalnum,1)
currFloor<-matrix(0,terminalnum,1)
already<-matrix(0,nodenumOfDyaker,1)
#
#############################
# INITIALIZING: "we go through the nodes at depth "depoftree""
# we make currFloor to be one over bottom floor and initialize 
# separy, boundary, atomsSepaAtom, atomsBounAtom
##############################
#
lkm<-0
curlkm<-0
curre<-beg
while(curre>0){
    lkm<-lkm+1
    atom<-AtomlistAtom[curre]
    node<-nodefinder[atom]
    #
    separy[node]<-lkm
    atomsSepaAtom[lkm]<-atom
    #
    exists<-parent[node]
    if (already[exists]==0){
         curlkm<-curlkm+1
         currFloor[curlkm]<-exists 
         already[exists]<-1
    }
    #
    curre<-AtomlistNext[curre]
}     # obs terminalnum=lkm
# initialize the rest
begsSepaBegs[1:terminalnum]<-seq(1,terminalnum)       
begsLeftBoun[1:terminalnum]<-seq(1,terminalnum)
begsRighBoun[1:terminalnum]<-seq(1,terminalnum)
#obs: we need not change 
#begsSepaNext, atomsSepaNext, atomsLBounNext, atomsRBounNext
#since at the beginning set consist only one member:
#pointer is always 0, since we do not have followers 
#
###########################
#
# START the MAIN LOOP
###########################
i<-d
while (i >= 2){
   j<-log(N[i],base=2)   #depth at direction d
   while (j>=1){
        nexlkm<-0
        k<-1
        while (k <= curlkm){
            node<-currFloor[k]  
            # we create simultaneously the upper floor
            exists<-parent[node]
            if (already[exists]==0){
                 nexlkm<-nexlkm+1
                 nextFloor[nexlkm]<-exists
                 already[exists]<-1
            }
####################################################
#            now we join childs            
####################################################
            leftbeg<-left[node]
            rightbeg<-right[node]
            direction<-i
            jg<-joingene(node,leftbeg,rightbeg,separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
direction,index)
#
separy<-jg$separy
begsSepaNext<-jg$begsSepaNext
begsSepaBegs<-jg$begsSepaBegs
begsLeftBoun<-jg$begsLeftBoun
begsRighBoun<-jg$begsRighBoun
atomsSepaNext<-jg$atomsSepaNext
atomsSepaAtom<-jg$atomsSepaAtom
atomsLBounNext<-jg$atomsLBounNext
atomsRBounNext<-jg$atomsRBounNext
           k<-k+1
        }
        j<-j-1
        curlkm<-nexlkm
        currFloor<-nextFloor
   }
   # now we move to the next direction, correct boundaries
   begsLeftBoun<-begsSepaBegs
   begsRighBoun<-begsSepaBegs
   #
   atomsLBounNext<-atomsSepaNext
   atomsRBounNext<-atomsSepaNext
   #
   i<-i-1
}
#########################
# ENO OF MAIN LOOP
#########################
#
#
###################
# LAST DIMENSION WILL BE handled, (because this contains root node
##################
   i<-1
   j<-log(N[i],base=2)   #depth at direction d
   while (j>=2){
        nexlkm<-0
        k<-1
        while (k <= curlkm){
            node<-currFloor[k]  
            # we create simultaneously the upper floor
            exists<-parent[node]
            if (already[exists]==0){
                 nexlkm<-nexlkm+1
                 nextFloor[nexlkm]<-exists
                 already[exists]<-1
            }
####################################################
#            now we join childs            
#if (right[parent[node]]==node)  #if node is right child 
####################################################
            leftbeg<-left[node]
            rightbeg<-right[node]
            direction<-1 
jg<-joingene(node,leftbeg,rightbeg,separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
direction,index)
#
separy<-jg$separy
begsSepaNext<-jg$begsSepaNext
begsSepaBegs<-jg$begsSepaBegs
begsLeftBoun<-jg$begsLeftBoun
begsRighBoun<-jg$begsRighBoun
#
atomsSepaNext<-jg$atomsSepaNext
atomsSepaAtom<-jg$atomsSepaAtom
atomsLBounNext<-jg$atomsLBounNext
atomsRBounNext<-jg$atomsRBounNext
           k<-k+1
        }
        j<-j-1
        curlkm<-nexlkm
        currFloor<-nextFloor
   }
#########################
# ROOT NODE, we do not anymore update boundaries
#########################
            k<-1
            node<-currFloor[k]  
  while (k <= 1){
####################################################
#            now we join childs            
#if (right[parent[node]]==node)  #if node is right child 
####################################################
            leftbeg<-left[node]
            rightbeg<-right[node]
            if ((leftbeg==0) || (separy[leftbeg]==0)){
                              #if left child does not exist    
                separy[node]<-separy[rightbeg]
            }
            else{   #eka else
                if ((rightbeg==0) || (separy[rightbeg]==0)){  
                              #right child does not exist  
                   separy[node]<-separy[leftbeg]
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
#separy[node]<- concatenate separy[leftbeg],separy[rightbeg]
###########
akku<-separy[leftbeg]
while (begsSepaNext[akku]>0){
  akku<-begsSepaNext[akku]
}                           
begsSepaNext[akku]<-separy[rightbeg]
#
separy[node]<-separy[leftbeg]
####################
#end of concatenating, handle next boundaries
###################
                    }
                    else{  #both children exist, both boundaries non-empty  
direction<-i
jc<-joinconne(leftbeg,rightbeg,separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
direction,index)  #direction<-i
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
 k<-k+1
}
###########################################
#          end of child joining
###########################################
###################################
# END of ROOT
###################################
#
###########################
###########################
totbegSepary<-separy[node]
return(list(totbegSepary=totbegSepary,
begsSepaNext=begsSepaNext,begsSepaBegs=begsSepaBegs,
atomsSepaNext=atomsSepaNext,atomsSepaAtom=atomsSepaAtom))
}














vera <-
function(numattable1, crds1, board1,  round1, currentbet,  mychips1, pot1, 
    roundbets, blinds1, chips1, ind1, dealer1, tablesleft){
    ## if any pair, suited anything, or if the smaller card is at least 9, then all in    
    a1 = 0
    if((crds1[1,1] == crds1[2,1]) || (crds1[1,2] == crds1[2,2]) || (crds1[2,1] > 8.5)) a1 = mychips1
    a1
} ## end of vera


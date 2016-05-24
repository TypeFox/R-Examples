timemachine <-
function(numattable1, crds1, board1,  round1, currentbet, 
    mychips1, pot1, roundbets, blinds1, chips1, ind1, dealer1, tablesleft){
    ## any pair 7 or higher
    ## AK, AQ
    ## AJ if nobody's all in yet.
    ## if less than 3 times bb & one card is Ten or higher, then 75%
    a1 = 0
    x = runif(1)
    if((crds1[1,1] == crds1[2,1]) && (crds1[1,1] > 6.5)) a1 = mychips1
    if((crds1[1,1] == 14) && (crds1[2,1] > 11.5)) a1 = mychips1
    if((crds1[1,1] == 14) && (crds1[2,1] == 11) && 
	(currentbet <= blinds1)) a1 = mychips1
    if((mychips1 < 3*blinds1) && (crds1[1,1] >= 10) && (x<.75)) a1 = mychips1
    a1
} ## end of timemachine


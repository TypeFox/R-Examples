xena <-
function(numattable1, crds1, board1,  round1, currentbet, mychips1, pot1,
    roundbets, blinds1, chips1, ind1, dealer1, tablesleft){
    ## if pair of 10s or higher, all in for sure, no matter what.
    ## if AK or AQ, all in with probability 75%.
    ## if pair of 7s or higher and there are 6 or fewer players 
    ##      at your table (including you), then all in. 
    ## if your chip count is less than twice the big blind, go all in with any cards.
    ## if nobody's raised yet:
    ##     ... and if there are 3 or fewer players left behind you,
    ##              then go all in with any pair or any ace.
    ##     ... and there's only 1 or 2 players behind you,
    ##              then go all in with any cards.    
    a1 = 0
    x = runif(1)               ## x is a random number between 0 and 1.
    y = max(roundbets[,1])     ## y is the maximum bet so far.
    big1 = dealer1 + 2
    if(big1 > numattable1) big1 = big1 - numattable1
    z = big1 - ind1
    if(z<0) z = z + numattable1 
    ## the previous 4 lines make it so z is the number of players left to act behind you.
    if((crds1[1,1] == crds1[2,1]) && (crds1[2,1] > 9.5)) a1 = mychips1
    if((crds1[1,1] == 14) && (crds1[1,2]>11.5) && (x<.75)) a1 = mychips1
    if((crds1[1,1] == crds1[2,1]) && (crds1[2,1] > 6.5) && (numattable1 < 6.5)) a1 = mychips1
    if(mychips1 < 2*blinds1) a1 = mychips1
    if(y <= blinds1){
	if((z < 3.5) && ((crds1[1,1] == crds1[2,1]) || (crds1[1,1] == 14))) a1 = mychips1
	if(z < 2.5) a1 = mychips1
    }
    a1
} ## end of xena


william <-
function(numattable1, crds1, board1,  round1, currentbet,  mychips1, pot1, 
    roundbets, blinds1, chips1, ind1, dealer1, tablesleft){
    ## if you only have less than 3 times the big blind, then all in.
    ## if AA, then all in.    
    ## if T9, then all in with 40% prob.
    ## if nobody's gone all in yet, then go all in.
    ## if KK or QQ, and less than 10 players at table, then all in.
    a1 = 0
    if(mychips1 < 3*blinds1) a1 = mychips1
    if((crds1[1,1] == 14) && (crds1[2,1] == 14)) a1 = mychips1
    if((crds1[1,1] == 10) && (crds1[2,1] == 9)){
	u1 = runif(1)
	if(u1 < .4) a1 = mychips1
	if(u1 > .4) a1 = 0
    }
    if(currentbet == blinds1) a1 = mychips1
    if((crds1[1,1] == crds1[2,1]) && (crds1[1,1] > 11.5) && (numattable1<10)) a1 = mychips1
    a1
} ## end of william


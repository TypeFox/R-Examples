yosef <-
function(numattable1, crds1, board1,  round1, currentbet,  mychips1, pot1, 
    roundbets, blinds1, chips1, ind1, dealer1, tablesleft){
    ## if pair of at least 9, then all in  
    ## if nobody's bet (or raised the blinds) yet, then bet (or raise the blinds) the minimum amount possible.
    ## note: if you're big blind, or if you've already called, then
    ## those cases have to be handled separately, since for instance if 
    ## someone's raised you exactly one more big blind, then you should fold, but it will look like 
    ## nobody's raised yet since it'll again be one big blind to you.
    a1 = 0
    bigb = dealer1 + 2
    if(bigb > numattable1) bigb = bigb - numattable1
    if(round1 == 1){
	if((roundbets[ind1,1] < blinds1 - .5) && (currentbet < blinds1+.5)){
	    a1 = min(2*blinds1, mychips1)
	} else if ((ind1 == bigb) && (currentbet < blinds1 - .5)) a1 = min(blinds1, mychips1)
    }
    if((round1 > 1.5) && (currentbet < .5)) a1 = min(blinds1, mychips1)
    if((crds1[1,1] == crds1[2,1]) && (crds1[2,1] > 8.5)) a1 = mychips1
    a1
} ## end of yosef


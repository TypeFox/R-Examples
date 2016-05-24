zelda <-
function(numattable1, crds1, board1,  round1, currentbet,  mychips1, pot1, 
    roundbets, blinds1, chips1, ind1, dealer1, tablesleft){
    a1 = 0 ## how much I'm gonna end up betting. Note that the default is zero.
    a2 = min(mychips1, currentbet) ## how much it costs to call
    
    if(round1 == 1){ ## pre-flop:
	## AK: Make a big raise if nobody has yet. Otherwise call.
	## AQ: call a small raise, or make one if nobody has yet.
	## AJ, AT, KQ, KJ, QJ: call a tiny raise.
	## A9, KT, K9, QT, JT, T9: call a tiny raise if in late position (within 2 of the dealer).
	## Suited A2-AJ: call a small raise.
	## 22-99: call a small raise.
	## TT-KK: make a huge raise. If someone's raised huge already, then go all in.
	## AA: make a small raise. If there's been a raise already, then double how much it is to you.
	
	a3 = 2*blinds1+1 ## how much a tiny raise would be
	a4 = 4*blinds1+1 ## how much a small raise would be
	a5 = max(8*blinds1,mychips1/4)+1 ## how much a big raise would be
	a6 = max(12*blinds1,mychips1/2)+1 ## how much a huge raise would be
	a7 = dealer1 - ind1
	if(a7 < -.5) a7 = a7 + numattable1 ## your position: a7 = how many hands til you're dealer
	
	if((crds1[1,1] == 14) && (crds1[2,1] == 13)){
	    a1 = max(a2,a5)
	}
	if((crds1[1,1] == 14) && (crds1[2,1] == 12)){
	    if(a2 < a4){
		a1 = a4
	    } else if(a2 > a5){
		a1 = 0
	    } else a1 = a2
	}
	if(((crds1[1,1] == 14) && ((crds1[2,1] < 11.5) && (crds1[2,1] > 9.5))) || 
	    ((crds1[1,1] == 13) && (crds1[2,1] > 10.5)) ||
	    ((crds1[1,1] == 12) && (crds1[2,1] == 11))){
	    if(a2 < a3) a1 = a2
	}
	if(((crds1[1,1] == 14) && (crds1[2,1] == 9)) || 
	    ((crds1[1,1] == 13) && ((crds1[2,1] == 10) || (crds1[2,1] == 9))) ||
	    ((crds1[1,1] == 12) && (crds1[2,1] == 10)) ||
	    ((crds1[1,1] == 11) && (crds1[2,1] == 10)) ||
	    ((crds1[1,1] == 10) && (crds1[2,2] == 9))){
	    if((a2 < a3) && (a7<2.5)) a1 = a2
	}
	if((crds1[1,2] == crds1[2,2]) && (crds1[1,1] == 14) && (crds1[2,1] < 11.5)){
	    if(a2<a4) a1 = a2
	    ## Note: this trumps the previous section, since it comes later in the code.
	}
	if((crds1[1,1] == crds1[2,1])){ ## pairs:
	    if(crds1[1,1] < 9.5){
		if(a2 < a4) a1 = a2
	    } else if(crds1[1,1] < 13.5){
		if(a2<a5) a1 = a5 else a1 = mychips1
	    } else {
		if(a2 < blinds1 + .5) a1 = a4 else a1 = min(2*a2,mychips1)
	    }
	}
    }
    if(round1 == 2){ ## post-flop: 
	## If there's a pair on the board and you don't have a set, then check/call up to small bet.
	## Same thing if there's 3-of-a-kind on the board and you don't have a full house or more.
	## If you have top pair or an overpair or two pairs or a set, make a big bet (call any bigger bet). 
	## Otherwise, if nobody's made even a small bet yet, then with prob. 20% make a big bluff bet.
	## If you're the last to decide and nobody's bet yet, then increase this prob. to 50%.
	## If you have an inside straight draw or flush draw then make a small bet (call any bigger bet).
	## If you have a straight or better, then just call.
	## Otherwise fold.
	
	a5 = min(sum(roundbets[,1]),mychips1) ## how much a big bet would be (prev round's pot size)
	a6 = min(.5*sum(roundbets[,1]),mychips1) ## how much a small bet would be	
	x = handeval(c(crds1[1:2,1], board1[1:3,1]), c(crds1[1:2,2], board1[1:3,2])) ## what you have
	x1 = handeval(c(board1[1:3,1]),c(board1[1:3,2])) ## what's on the board
	y = straightdraw1(c(crds1[1:2,1], board1[1:3,1]))
	z = flushdraw1(c(crds1[1:2,2], board1[1:3,2]))
	topcard1 = max(board1[1:3,1])
	a7 = runif(1) ## random number uniformly distributed between 0 and 1
	a8 = (1:numattable1)[roundbets[,1] == roundbets[ind1,1]] ## others who can still bet with you
	## The next 5 lines may seem weird, but the purpose is explained in the next comment:
	a9 = a8 - dealer1
	for(i in 1:length(a9)) if(a9[i]<.5) a9[i] = a9[i] + numattable1
	a10 = ind1 - dealer1
	if(a10 < .5) a10 = a10 + numattable1
	a11 = 2*(a10 == max(a9))   ## So a11 = 2 if you're last to decide; otherwise a11 = 0.
	
	if((x1 > 1000000) && (x < 3000000)){
	    if(a2 < a6) a1 = a2
	} else if((x1 > 3000000) && (x < 6000000)){
	    if(a2 < a6) a1 = a2
	} else if(x > 1000000 + 15^3*topcard1){
	    a1 = max(a5,a2)
	} else if((a2 < a6) && ((a7 < .20) || ((a7 < .50) && (a11>1)))){
	    a1 = a6
	}
	if((y == 4) || (z == 4)) a1 = max(a6, a2)
	if(x > 4000000) a1 = a2
    }
    if(round1 == 3){ ## after turn: 
	## If there's a pair on the board and you don't have a set, then check/call up to small bet.
	## Same thing if there's 3-of-a-kind on the board and you don't have a full house or more.
	## Otherwise, if you have top pair or better, go all in.
	## If you had top pair or overpair but now don't, then check/call a medium bet but fold to more.
	## If you have an inside straight draw or flush draw then check/call a medium bet as well.
	## Otherwise check/fold.
	a6 = min(1/3*sum(roundbets[,1:2]),mychips1) ## small bet (1/3 of prev round's pot size)
	a5 = min(.75*sum(roundbets[,1:2]),mychips1) ## medium bet (3/4 of prev round's pot size)
	x = handeval(c(crds1[1:2,1], board1[1:4,1]), c(crds1[1:2,2], board1[1:4,2])) ## what you have
	x1 = handeval(c(board1[1:4,1]),c(board1[1:4,2])) ## what's on the board
	y = straightdraw1(c(crds1[1:2,1], board1[1:4,1]))
	z = flushdraw1(c(crds1[1:2,2], board1[1:4,2]))
	topcard1 = max(board1[1:4,1])
	oldtopcard1 = max(board1[1:3,1])
	if((x1 > 1000000) && (x < 3000000)){
	    if(a2 < a6) a1 = a2
	} else if((x1 > 3000000) && (x < 6000000)){
	    if(a2 < a6) a1 = a2
	} else if(x > 1000000 + 15^3*topcard1){
	    a1 = mychips1
	} else if(x > 1000000 + 15^3*oldtopcard1){
	    if(a2 < a5) a1 = a2
	} else if((y == 4) || (z == 4)){
	    if(a2 < a5) a1 = a2
	}
    }
    if(round1 == 4){ ## after river: 
	## If there's a pair on the board and you don't have a set, then check/call up to small bet.
	## Same thing if there's 3-of-a-kind on the board and you don't have a full house or more.
	## Otherwise, if you have two pairs or better, go all in.
	## If you have one pair, then check/call a small bet.
	## With nothing, go all-in with probability 10%; otherwise check/fold.
	a6 = .45+runif(1)/10  ## random number between .45 and .55
	a5 = min(a6*sum(roundbets[,1:3]),mychips1) ## small bet: around 1/2 of pot size; VARIES RANDOMLY
	x = handeval(c(crds1[1:2,1], board1[1:5,1]), c(crds1[1:2,2], board1[1:5,2]))
	x1 = handeval(c(board1[1:5,1]),c(board1[1:5,2])) ## what's on the board
	if((x1 > 1000000) && (x < 3000000)){
	    if(a2 < a5) a1 = a2
	} else if((x1 > 3000000) && (x < 6000000)){
	    if(a2 < a5) a1 = a2
	} else if(x > 2000000){
	    a1 = mychips1
	} else if(x > 1000000){
	    if(a2 < a5) a1 = a2
	} else if(runif(1)<.10){
	    a1 = mychips1
	}
    }
    round(a1)
} ## end of zelda


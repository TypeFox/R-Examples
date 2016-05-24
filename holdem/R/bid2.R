bid2 <-
function(numattable1, playerseats1, blinds1, dealer1, b3, b4, round1, ntable1, decision1){
   if(b4$all1 > 1) return(b4)
   in1 = b4$i1
   if(sum(in1)<1.5) return(b4)
   betlist1 = c(0)
   indlist1 = c(0)
   outlist1 = c(0)
   chips1 = b4$c1
   bet1 = rep(0,numattable1)   
   ind1 = dealer1+1
   if(ind1>numattable1) ind1 = 1
   currentbet = 0
   better1 = dealer1+1
   if(better1 > numattable1) better1 = 1
   v1 = c(0,3,4,5)[round1] 
   board1 = matrix(rep(0,10),ncol=2)
   board1[1:v1,1] = b3$brdnum1[1:v1]
   board1[1:v1,2] = b3$brdsuit1[1:v1]
   roundbets = b4$rb
   pot1 = b4$p1  
   stp = 0
   while(stp < 1){
       out1 = 0
       if(in1[ind1] > 0.5){   
	   crds1 = matrix(c(b3$plnum1[ind1,],b3$plsuit1[ind1,]),ncol=2,byrow=F)
	   # cat("...",playerseats1[ind1],"'s turn...")
	   bmax1 = max((bet1[-ind1] + chips1[-ind1])[in1[-ind1]>.5])
	   b1 = round(decision1[[playerseats1[ind1]]](numattable1, 
	            crds1, 
		    board1, 
		    round1, 
		    currentbet - bet1[ind1], 
		    chips1[ind1],
		    pot1, 
		    roundbets,
		    blinds1[2], 
		    chips1,
		    ind1,
		    dealer1,
		    ntable1))
	   # cat("\n Seat ", ind1,": b1 was ",b1," and it was ",currentbet-bet1[ind1]," to him.")
	   if(b1 > chips1[ind1]) b1 = chips1[ind1] ## if bet is more than you have, fix that.
	   if(b1 > bmax1 - bet1[ind1]) b1 = bmax1 - bet1[ind1] ## can't bet more than anyone else who's in has left
	   ## if bet is between 0.5 and the amount to you, make it a call.
	   if((b1 > 0.5) && (b1 < currentbet - bet1[ind1])) b1 = min(chips1[ind1],currentbet - bet1[ind1])
	   ## if bet is a raise of less than the big blind, make it a raise of the big blind.
	   raiseamt1 = b1 - (currentbet - bet1[ind1])
	   if((raiseamt1 > 0.5) && (raiseamt1 < blinds1[2])) b1 = min(blinds1[2] + currentbet - bet1[ind1],chips1[ind1])
	   if(b1 > currentbet - bet1[ind1]+.5){ ## raise
	       if(b1 > chips1[ind1]) b1 = chips1[ind1]
	       currentbet = b1 + bet1[ind1]
	       better1 = ind1
	       pot1 = pot1 + b1
	       roundbets[ind1,round1] = roundbets[ind1,round1] + b1
	       bet1[ind1] = roundbets[ind1,round1]
	       in1[ind1] = 1
	       chips1[ind1] = chips1[ind1] - b1
	   } else if(b1 == min(chips1[ind1],currentbet-bet1[ind1])){ ## call/check
	       pot1 = pot1 + b1
	       roundbets[ind1,round1] = roundbets[ind1,round1] + b1
	       bet1[ind1] = roundbets[ind1,round1]
	       in1[ind1] = 1
	       chips1[ind1] = chips1[ind1] - b1
	   } else if((chips1[ind1]>0.5) && (b1 < min(chips1[ind1],currentbet-bet1[ind1]))){ ## fold
	       in1[ind1] = 0
	       out1 = 2
	   }
	   betlist1 = c(betlist1,bet1[ind1])
	   indlist1 = c(indlist1,ind1)
	   outlist1 = c(outlist1, out1)
       }
       ind1 = ind1 + 1
       if(ind1 > numattable1) ind1 = 1
       if(better1 == ind1) stp = 2
       if(sum(in1) < 1.5) stp = 2  
   }
   z3 = 0 ## now see if all the betting is over: if so, let z3 = 2.
   if(sum(chips1[c(1:numattable1)[in1 > .5]] > .5) < 1.5) z3 = 2
   list(i1=in1,p1=pot1,c1=chips1,rb=roundbets,all1=z3,bl1 = betlist1,il1 = indlist1,out1=outlist1)
} ## end of bid2


bid1 <-
function(numattable1, playerseats1, chips1, blinds1, dealer1, b3, ntable1, decision1){
   round1 = 1
   in1 = 1*(chips1>0.5)
   bet1 = rep(0,numattable1)   
   betlist1 = c(0)
   indlist1 = c(0)
   outlist1 = c(0)
   ind1 = dealer1+1
   if(ind1>numattable1) ind1 = 1
   j = 0
   prevbet = 0
   currentbet = 0
   better1 = 0
   board1 = matrix(rep(0,10),ncol=2)
   board1[1:3,1] = b3$brdnum1[1:3]
   board1[1:3,2] = b3$brdsuit1[1:3]
   roundbets = matrix(rep(0,(4*numattable1)),ncol=4)
   pot1 = 0  
   stp = 0
   while(stp < 1){
       myout1 = 0
       j = j+1
       if(j==1){
	   bet1[ind1] = min(chips1[ind1],blinds1[1])
	   prevbet = bet1[ind1]
	   better1 = ind1
	   pot1 = bet1[ind1]
	   roundbets[ind1,1] = bet1[ind1]
	   chips1[ind1] = chips1[ind1] - bet1[ind1]
	   betlist1 = c(betlist1,bet1[ind1])
	   indlist1 = c(indlist1,ind1)
	   outlist1 = c(outlist1,myout1)
       }
       if(j==2){
	   bet1[ind1] = min(chips1[ind1],blinds1[2])
	   currentbet = max(prevbet,blinds1[2])
	   better1 = ind1
	   pot1 = pot1 + bet1[ind1]
	   roundbets[ind1,1] = bet1[ind1]
	   chips1[ind1] = chips1[ind1] - bet1[ind1]
	   betlist1 = c(betlist1,bet1[ind1])
	   indlist1 = c(indlist1,ind1)
	   outlist1 = c(outlist1,myout1)
       }
       if(j==3){better1 = ind1}  # so that big blind can raise!
       if((j > 2.5) && (in1[ind1] > 0.5)){
	   crds1 = matrix(c(b3$plnum1[ind1,],b3$plsuit1[ind1,]),ncol=2)
	   bmax1 = max((bet1[-ind1] + chips1[-ind1])[in1[-ind1]>.5])
	   # cat(".a. Seat ",ind1,"'s turn, bmax is ",bmax1,", lastbetter is ",better1,"...\n")
	   if(bmax1 < currentbet) currentbet = bmax1
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
	   # cat("Now b1 = ",b1,". Currentbet = ",currentbet,".\n")
	   if(b1 > currentbet-bet1[ind1]+.5){ ## raise
	       prevbet = currentbet
	       currentbet = b1+bet1[ind1]
	       # cat("curbet = ",currentbet,".\n")
	       better1 = ind1
	       pot1 = pot1 + b1
	       roundbets[ind1,round1] = roundbets[ind1,round1] + b1
	       bet1[ind1] = roundbets[ind1,round1]
	       in1[ind1] = 1
	       chips1[ind1] = chips1[ind1] - b1
	   } else if(b1 == min(chips1[ind1],currentbet-bet1[ind1])){ ## call
	       pot1 = pot1 + b1
	       roundbets[ind1,round1] = roundbets[ind1,round1] + b1
	       bet1[ind1] = roundbets[ind1,round1]
	       in1[ind1] = 1
	       chips1[ind1] = chips1[ind1] - b1
	   } else if((chips1[ind1]>0.5) && (b1 < min(chips1[ind1],currentbet-bet1[ind1])-.5)){ ## fold
	       in1[ind1] = 0
	       myout1 = 2
	   }
	   betlist1 = c(betlist1,bet1[ind1])
	   indlist1 = c(indlist1,ind1)
	   outlist1 = c(outlist1,myout1)
       }
       ind1 = ind1 + 1
       if(ind1 > numattable1) ind1 = 1
       if(better1 == ind1) stp = 2
       # cat("\n\n who:",numattable1,playerseats1, chips1,".\n")
       if(sum(in1) < 1.5) stp = 2  
   }
   z3 = 0 ## now see if all the betting is over: if so, let z3 = 2.
   if(sum(chips1[c(1:numattable1)[in1 > .5]] > .5) < 1.5) z3 = 2
   # cat("\n Round 1 over.  z3 = ",z3,". in1 = ",in1,".\n")
   list(i1=in1,p1=pot1,c1=chips1,rb=roundbets,all1=z3,bl1 = betlist1,il1 = indlist1,out1=outlist1)
} ## end of bid1


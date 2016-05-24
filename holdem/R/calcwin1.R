calcwin1 <-
function(numattable1,playerseats1, b3, b7){
## First, for everyone who's in, evaluate their hands.
    qual1 = rep(0,numattable1)
    inmuch1 = rep(0,numattable1)
    for(i in c(1:numattable1)){
	if(b7$i1[i] > .5) {
	    qual1[i] = handeval(c(b3$brdnum1[1:5],b3$plnum1[i,1:2]),
		c(b3$brdsuit1[1:5],b3$plsuit1[i,1:2]))
	}
	inmuch1[i] = sum(b7$rb[i,])
    }
    z1 = order(unique(inmuch1[inmuch1 > 0.5]))
    difbets = sum(unique(inmuch1[inmuch1 > 0.5])>0) ## How many unique betting amounts.
    pot1 = b7$p1
    chips1 = b7$c1
    prevamount1 = 0    
    for(j in c(1:difbets)){
	amount1 = unique(inmuch1[inmuch1 > 0.5])[z1[j]] ## Start with the least amount
	n1 = sum(inmuch1 >= amount1)  ## How many bet at least that much.
	pl1 = c(1:numattable1)[inmuch1 >= amount1]  ## Who bet at least that much
	# winner1 = pl1[order(qual1[pl1],decreasing=TRUE)[1]]  
	winner1 = pl1[c(1:n1)[qual1[pl1]==max(qual1[pl1])]] ## Who has best hand (may be vector, if ties)
	more1 = amount1 - prevamount1
	for(k in winner1) chips1[k] = chips1[k] + round(n1*more1/length(winner1))
	pot1 = pot1 - n1*more1
	prevamount1 = amount1
    } 
    chips1
} ## end of calcwin1


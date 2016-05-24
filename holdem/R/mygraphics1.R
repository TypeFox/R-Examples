mygraphics1 <-
function(numattable1, playerseats1, chips1, blinds1, dealer1,b3,b4,b5,b6,b7,chips2,
		       ntable1,myfast1,chipstart1,name1,lowercut1){
    nplayers1 = length(name1)
    cardname1 = c(as.character(1:9),"T","J","Q","K","A")
    suitname1 = c(2,3,4,6)
    plot(c(0,10+10*numattable1),c(0,100),type="n",xlab="",ylab="",xaxt="n",yaxt="n")
    a1 = paste("Blinds are",blinds1[1],"and",blinds1[2],".")
    text(10,95,a1)
    if(ntable1>1.5) a2 = paste(ntable1,"tables left.") else a2 = "1 table left."
    text(10*numattable1,95,a2)
    text(10*dealer1,90,"D")    
    text(5,5,"click to continue",cex=.8)
    for(j in c(1:numattable1)){
	i = 2 + dealer1 + j
	if(i > numattable1 + .5) i = i - numattable1
	if(i > numattable1 + .5) i = i - numattable1
	text(10*i,80,as.character(name1[playerseats1[i]]),cex=1+.1*b7$i1[i],col=1)
	a1 = paste("(",chips1[i],")")
	text(10*i,75,a1)
	text(10*i,50,"BETS:")
    }
    ## pre-flop bets:
    writebets1(b4,1,numattable1,b3,playerseats1,chips1,chips2,0,myfast1,name1)
	## Now flop:
    for(j in c(1:3)){
	text(j/6*(10+10*numattable1),15,cardname1[b3$brdnum1[j]],col=suitname1[b3$brdsuit1[j]],cex=2)
    }
    if(myfast1<1) locator(1)
    over1 = b4$all1
    writebets1(b5,2,numattable1,b3,playerseats1,chips1,chips2,over1,myfast1,name1)
    ## Now turn:
    over1 = b5$all1
    j = 4
    text(j/6*(10+10*numattable1),15,cardname1[b3$brdnum1[j]],col=suitname1[b3$brdsuit1[j]],cex=2)
    if(myfast1<1) locator(1)
    writebets1(b6,3,numattable1,b3,playerseats1,chips1,chips2,over1,myfast1,name1)   
    ## Now river:
    over1 = b6$all1
    j = 5
    text(j/6*(10+10*numattable1),15,cardname1[b3$brdnum1[j]],col=suitname1[b3$brdsuit1[j]],cex=2)
    if(myfast1<1) locator(1)
    writebets1(b7,4,numattable1,b3,playerseats1,chips1,chips2,over1,myfast1,name1)
    for(i in c(1:numattable1)) if(chips1[i] != chips2[i]) text(10*i, 70, chips2[i])
    for(i in c(1:numattable1)){
	if(chips2[i]>chips1[i]){
	    text(10*i,85,"$",cex=3)
	    text(10*i,80,as.character(name1[playerseats1[i]]),cex=1+.1*b7$i1[i],col=5)
	}
	if(chips2[i] < .5){
	    text(10*i,80,as.character(name1[playerseats1[i]]),cex=1+.1*b7$i1[i],col=2)
	}
    }
    if(myfast1<1) locator(1)
    plot(c(0,nplayers1+1),c(lowercut1,chipstart1*nplayers1),type="n",xlab="player number",ylab="chips",log="y")
} ## end of mygraphics1


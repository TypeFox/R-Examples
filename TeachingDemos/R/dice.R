"dice" <-
function(rolls=1, ndice=2, sides=6, plot.it=FALSE, load=rep(1,sides))
# Simulate the tossing of some dice.
# rolls is the number of times to roll the dice
# ndice is the number of dice to roll each time
# sides is the number of sides to the dice
# load is how the dice are loaded, can be though of as odds
{
	temp <- matrix( sample(sides, ndice*rolls, TRUE, load), ncol=ndice )
	temp <- as.data.frame(temp)
	
	names(temp) <- c("Red","Green","Blue","Black","Yellow","Purple",
                         "Orange","Brown","Grey","White")[1:ndice]
	
	#if(ndice==1) return(temp$Red)
	
	oldClass(temp) <- c("dice","data.frame")
	
	if(plot.it){
		plot.dice(temp)
		return(invisible(temp))
	}
	temp	
}


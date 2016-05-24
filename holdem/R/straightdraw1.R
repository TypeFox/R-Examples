straightdraw1 <-
function(x){
    ## returns 4 is there are 2 possibilities for a straight. 
    ## returns 2 for a gutshot straight draw.
    ## returns 0 otherwise
    ## Note: returns 26 if you already have a straight!
    a1 = 0
    for(i in c(2:14)){
	a2 = straight1(c(i,x))
	if(a2 > .5) a1 = a1 + 1
    }
    a1 * 2 
}	## end of straightdraw1


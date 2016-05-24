`D2ACWmat` <-
function(J, filter.number = 10., family = "DaubLeAsymm", switch = "direction",
OPLENGTH = 100000.)
{
#
# Program to generate an array which stores the 2-D autocorrelation waveles
#
# This program does not play a role in the main LS2W suite, but has been included in case others
# find it useful!
#
# It works by calling the programs D2ACWmat.l and D2ACWmat.d
#
##############################
#
# We strart with some checks:
#
if(J >= 0.) stop("Sorry! J must be a negative integer.\n")
####
if(switch != "direction" && switch != "level") {
stop("Sorry, but switch can only take the values level and direction.\n")
}
if(switch == "direction") {
return(D2ACWmat.d(J, filter.number = filter.number, family = 
family, OPLENGTH = OPLENGTH))
}
else {
return(D2ACWmat.l(J, filter.number = filter.number, family = 
family, OPLENGTH = OPLENGTH))
}
}


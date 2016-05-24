"formatA" <-
function(x,digits=2, FUN=round,...){
	noquote(format(FUN(x, digits=digits),...))
}


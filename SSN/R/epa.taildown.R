epa.taildown <-
function(dist.hydro, a.mat, b.mat, parsil = parsil, range1 = range1, 
	useTailDownWeight, weight = NULL)
{
	flow.connect <- b.mat == 0
	V <- parsil/(16*range1^5)*(dist.hydro - range1)^2*
		(16*range1^3*(dist.hydro*0 + 1) +
		17*range1^2*dist.hydro -
		2*range1*dist.hydro^2 - 
		dist.hydro^3)* 
		(dist.hydro < range1)*flow.connect +
	parsil/(16*range1^5)*(a.mat - range1)^2*
		(16*range1^3*(dist.hydro*0 + 1) +
		17*range1^2*a.mat -
		15*range1^2*b.mat -
		20*range1*b.mat^2 -
		2*range1*a.mat^2 +
		10*range1*a.mat*b.mat +
		5*a.mat^2*b.mat -
		a.mat^3 -
		10*a.mat*b.mat^2)*
		(a.mat < range1)*(1 - flow.connect)
	if(useTailDownWeight == TRUE) V <- V*weight
	V
}


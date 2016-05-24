fii2fi <-
function(DESIGN,fii,fi){
	Dsup <- fastEucCalc(fii,fi)
	minD <- apply(Dsup,1,min)
	Group_Assigned=Re(Dsup==repmat(minD,1,ncol(DESIGN)))
	return(list(distances=Dsup,assignments=Group_Assigned,confusion=t(Group_Assigned) %*% DESIGN))
}

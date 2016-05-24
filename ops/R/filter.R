filter <-
function(x,ia,ib,t=-1) {
	y=x[x[,ia]>t & x[,ib]>t,c(ia,ib)];
	y;	
}


belongtoint <-
function(bin.vals.vect,int.x,int.y,vect){

	xval<-bin.vals.vect[1];
	yval<-bin.vals.vect[2];	



	sum(vect[,1] >= xval & vect[,1] < (xval+int.x) & vect[,2] >= yval & vect[,2] < (yval+int.y));
 
}


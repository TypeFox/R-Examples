addem <-
function(startingvals,xmatlen,ymatlen,data){

	x<-startingvals[1];
	y<-startingvals[2];


	#print("here");
	#print(data[x:(x+xmatlen-1),y:(y+ymatlen-1)]);

	sum(data[x:(x+xmatlen-1),y:(y+ymatlen-1)])^2;




}


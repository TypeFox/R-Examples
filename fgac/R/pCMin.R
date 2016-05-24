"pCMin" <-
function(theta,delta,s,t)
{if(missing(theta) | missing(delta)){theta<-1;delta<-1};
	pc<-pmax(0,(s+t-1));
	resu<-pc
}


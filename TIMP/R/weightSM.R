"weightSM" <-
function (weight, typ, sm, pt) 
{
#weightSM
	if(typ==0){
		if(pt < weight[1]){ 
			lina<-approxfun(c(weight[1],weight[1]-sm),
                                        c(weight[3],1))
			weight[3]<-lina(pt)
		}
		else{ 
			lina<-approxfun(c(weight[1],weight[2]+sm),
                                        c(weight[3],1))
			weight[3]<-lina(pt)
			
		}
	}

	weight[3]	
}


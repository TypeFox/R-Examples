#### tarjet function to be minimized for IRT CCM parameter-linking

target = function(z,parm,common,model,icc,meth,D,...){

grid.theta = seq(-4,4, 0.05)

	if(model=="1PL"){
		A<-1
		B<-z
		if(icc=="cloglog"){
			Pij <- Pij.cloglog
					}
		else{
			Pij <- Pij.logistic
		    }
	}
	else{
		A<-z[1]
		B<-z[2]
			Pij <- Pij.logistic
		}


parm$aIj = parm$aIj/A
parm$bIj = A*parm$bIj + B
		
		
		Hdiff.theta = 0 
		SL1 = 0
		SL2 = 0
		for(j in 1:length(common)){
				Hdiff.theta = Hdiff.theta + (Pij(as.vector(t(parm[j, 1:3])), grid.theta,D) - Pij(as.vector(t(parm[j, 4:6])), grid.theta,D))^2
				SL1 = SL1 + Pij(as.vector(t(parm[j, 1:3])), grid.theta,D)
				SL2 = SL2 + Pij(as.vector(t(parm[j, 4:6])), grid.theta,D)

		Hcrit = sum(Hdiff.theta)
		SLdiff = (SL1 - SL2)^2
		SLcrit = sum(SLdiff)
		}
if(meth=="H") res<-Hcrit
else if(meth=="SL") res<-SLcrit

return(res)

}


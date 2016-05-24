##' Decision Rules
##'
##' Different rules for making decisions in the framework of belief functions
##' @export
##' @param mass The matrix containing the masses. Each column represents a piece of mass.
##' @param criterion The decision baseline:
##'
##'     criterion=1 maximum of the plausibility
##'
##'		criterion=2 maximum of the credibility
##'
##'		criterion=3 maximum of the credibility with rejection
##'
##'		criterion=4 maximum of the pignistic probability
##'
##'		criterion=5 Appriou criterion (decision onto \eqn{2^\Theta})
##' @param r The parameter in BayesianMass function. If criterion 5 is used, it should be given. 
##' Otherwise it will be set to the default value 0.5.
##' @return The decision vector. E.g., in classification problem, class labels. 
##' @examples
##' m1=c(0,0.4, 0.1, 0.2, 0.2, 0, 0, 0.1);
##' m2=c(0,0.2, 0.3, 0.1, 0.1, 0, 0.2, 0.1);
##' m3=c(0.1,0.2, 0, 0.1, 0.1, 0.1, 0, 0.3);
##' 
##' m3d=discounting(m3,0.95);
##' 
##' M_comb_Smets=DST(cbind(m1,m2,m3d),1);
##' M_comb_PCR6=DST(cbind(m1,m2),8);
##' 
##' class_fusion=decisionDST(M_comb_Smets,1)
##' class_fusion=decisionDST(M_comb_PCR6,1)
##' class_fusion=decisionDST(M_comb_Smets,5,0.5)
##' class_fusion=decisionDST(cbind(M_comb_Smets,M_comb_PCR6),1)
##' 
decisionDST <- function (mass,criterion,r=0.5){

  if (is.vector(mass) || (is.matrix(mass) && nrow(mass) == 1)) {
    mass = matrix(mass,, 1)
  }
	lm=nrow(mass);
	nbvec_test=ncol(mass);
	nbclasses=round(log2(lm));



	class_fusion=c();
	for(k in 1:nbvec_test){
		masstmp=mass[,k];

		if(criterion==1){
			# case 1
			plau=mtopl(masstmp);
			ii=1:nbclasses;
			plau_singl=plau[1+2^(ii-1)];
			indice=which.max(plau_singl);
			class_fusion=c(class_fusion,indice);
		}else if(criterion==2||criterion==3){
			# case {2, 3}
      # browser()
			croy=mtobel(masstmp);
			ii=1:nbclasses;
			croy_singl=croy[1+2^(ii-1)];
	        valuemax=max(croy_singl);
			indice=which.max(croy_singl);
				if(criterion==3){
					indice_complementaire=0;
					for (i in seq(nbclasses,indice,by=-1)){
						indice_complementaire=indice_complementaire+2^(nbclasses-(nbclasses-i+1));
				    }
			    if (valuemax>=croy[indice_complementaire]){
						class_fusion=c(class_fusion,indice);
				}else{
						class_fusion=c(class_fusion,0);
				}
				}else{
				  class_fusion=c(class_fusion,indice);
        }
		}else if(criterion==4){
			# case 4
				pign=mtobetp(t(masstmp));
				indice=which.max(pign);
				class_fusion=c(class_fusion,indice);
		}else if(criterion==5){
			# case 5
			  

				plau=mtopl(t(masstmp));
				lambda=1;
				md=BayesianMass(lambda,r,nbclasses);
				indice=which.max(plau*md);
				class_fusion=c(class_fusion,indice);

		 }else{
				stop('ACCIDENT: The critertion given is not right\n')
		}
  }
	return(class_fusion)
}


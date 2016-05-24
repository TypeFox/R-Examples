##' Combination rules
##' 
##' Different rules to combine masses
##' 
##' @export
##' 
##' @param MassIn The matrix containing the masses. Each column represents a
##' piece of mass.
##' @param criterion The combination criterion:
##' 
##' criterion=1 Smets criterion (conjunctive combination rule)
##' 
##' criterion=2 Dempster-Shafer criterion (normalized)
##' 
##' criterion=3 Yager criterion
##' 
##' criterion=4 Disjunctive combination criterion
##' 
##' criterion=5 Dubois criterion (normalized and disjunctive combination)
##' 
##' criterion=6 Dubois and Prade criterion (mixt combination), only for Bayesian masses whose focal elements are singletons
##' 
##' criterion=7 Florea criterion
##' 
##' criterion=8 PCR6
##' 
##' criterion=9 Cautious Denoeux Min for functions non-dogmatics
##' 
##' criterion=10 Cautious Denoeux Max for separable masses
##' 
##' criterion=11 Hard Denoeux for functions sub-normal
##' 
##' criterion=12 Mean of the bbas
##' 
##' criterion=13 LNS rule, for separable masses
##' 
##' criterion=131 LNSa rule, for separable masses
##' @param TypeSSF If TypeSSF = 0, it is not a SSF, the general case. If TypeSSF = 1, a SSF with a singleton as a focal element. If TypeSSF = 2, a SSF with any subset of \eqn{\Theta} as a focal element. 
##' @return The combined mass vector. One column. 
##' @examples
##' 
##' m1=c(0,0.4, 0.1, 0.2, 0.2, 0, 0, 0.1);
##' m2=c(0,0.2, 0.3, 0.1, 0.1, 0, 0.2, 0.1);
##' m3=c(0.1,0.2, 0, 0.1, 0.1, 0.1, 0, 0.3);
##' 
##' m3d=discounting(m3,0.95);
##' 
##' M_comb_Smets=DST(cbind(m1,m2,m3d),1);
##' M_comb_Smets
##' M_comb_PCR6=DST(cbind(m1,m2),8);
##' M_comb_PCR6
##' M_comb_LNS = DST(cbind(m1,m2),13);
##' M_comb_LNS 
##' M_comb_LNSa = DST(cbind(m1,m2),131);
##' M_comb_LNSa 
##' 
##' n1 = 5
##' ThetaSize = 3
##' mass_mat = matrix(0, 2^ThetaSize, n1 + 1);
##' mass_mat[2, 1 : n1] = c(0.12, 0.16, 0.15, 0.11, 0.14) 
##' mass_mat[3, n1 + 1] = 0.95;
##' mass_mat[8, ] = 1 - colSums(mass_mat)
##' mass_ssf_mat = mass_mat[c(2^(1:ThetaSize-1)+1, 8), ]
##' # the following three functions could produce the same results
##' DST(mass_mat, 13)
##' DST(mass_mat, 13, TypeSSF = 2)
##' DST(mass_ssf_mat, 13, TypeSSF = 1)


DST <- function(MassIn, criterion, TypeSSF = 0){

	 
     n = nrow(MassIn);
	 m = ncol(MassIn);

     if(criterion %in% c(4, 5, 6, 7)){
		 b_mat = apply(MassIn, 2, mtob);
		 b = apply(b_mat, 1, prod)
	 }
	 
	 if(criterion %in% c(1, 2, 3, 6, 7)){
		 q_mat = apply(MassIn, 2, mtoq);
		 q = apply(q_mat, 1, prod)
	 }


		if(criterion==1){
			#Smets criterion
			Mass=qtom(q);
			Mass[1]=1-sum(Mass[2:length(Mass)]); #In case: with very high conflict Mass(1) could be >1 !
		}else if (criterion==2){
			#Dempster-Shafer criterion (normalized)
			Mass=qtom(q);
			Mass=Mass/(1-Mass[1]);                 
			Mass[1]=0;
		}else if (criterion==3){
			#Yager criterion
			Mass=qtom(q);
			Mass[n]=Mass[n]+Mass[1];
			Mass[1]=0; 
		}else if (criterion==4){
			# disjunctive combination criterion
			Mass=btom(b);
		}else if (criterion==5){
			#Dubois criterion (normalized and disjunctive combination)
			Mass=btom(b);
			Mass=Mass/(1-Mass[1]);
			Mass[1]=0;
		}else if (criterion==6){
			#Dubois and Prade criterion (mixt combination). Only if the focal
			#element are the singletons
			Mass=qtom(q);
			Mass[n]=0;
	   
			Mass_disjonc=btom(b);
			Mass_disjonc[n]=0;
			for(i in 1:floor(log2(n))){
				Mass_disjonc[1+2^(i-1)]=0;
     		}
			Mass=Mass+Mass_disjonc;
			Mass[1]=0;
			Mass[n]=1-sum(Mass);
		}else if (criterion==7){
			#Florea criterion
			Mass_conjonc=qtom(q);
			Mass_disjonc=btom(b);
			
			k=Mass_conjonc[1];
			x=0.9;
			alpha=log((1+x)/(k+x))/log((1+x)/x);
			beta=log((1+x)^k*(k+x)^(1-k)/x)/log((1+x)/x);
			
			Mass=alpha*Mass_disjonc+beta*Mass_conjonc;
		
			Mass[1]=0;
		}else if (criterion==8){
			#PCR6 combination Martin & Osswald Criteria
		     re=PCR6(MassIn);	
			 Conf=re$Conf;
			 Mass=re$Mass;
		 	 Mass=t(Mass);
			
		}else if (criterion==9){
			#Cautious Denoeux min for fonctions non-dogmatic
			wtot = apply(MassIn, 2, mtow)
			w=apply(wtot,1,min)
			Mass=wtom(w);
		}else if (criterion==10){
			#Cautious Denoeux max only for separable fonctions
			wtot = apply(MassIn, 2, mtow)
			w=apply(wtot,1,max)
			Mass=wtom(w);
		}else if (criterion==11){ 
			# Hard Denoeux for fonctions sub-normal
			vtot = apply(MassIn, 2, mtov)
			v=apply(vtot,1,min)
			Mass=vtom(v);
		}else if (criterion==12){
			# mean of the masses
			 Mass = apply(MassIn,1,mean);
		 
		}else if (criterion==13){
	        # LNS rule	
			if(TypeSSF == 0){
			  Mass = tCombine(MassIn, mygamma = 1)
			}else if(TypeSSF == 1){
			  Mass = tCombine_SSF(MassIn, mygamma = 1, singleton = TRUE)
			}else if(TypeSSF == 2){
			  Mass = tCombine_SSF(MassIn, mygamma = 1)
			}
		}else if (criterion == 131){
	        # LNSa rule	
			if(TypeSSF == 0){
			  Mass = tCombine(MassIn, approximate = TRUE)
			}else if(TypeSSF == 1){
			  Mass = tCombine_SSF(MassIn, mygamma = 1, approximate = TRUE, singleton = TRUE)
			}else if(TypeSSF == 2){
			  Mass = tCombine_SSF(MassIn, mygamma = 1, approximate = TRUE)
			}
		}else{
			stop('Accident: in DST choose of criterion: uncorrect\n');
        }
	return(matrix(Mass,,1))
}


tCombine <- function(MassIn, mygamma = 1, ifnormalize = FALSE, ifdiscount = TRUE, eta = 0, approximate = FALSE){
    
    ## LNS rule.  Also can be used for conjunctive rule, cautious rule and DS rule
	## only for seperable mass
	## if want to use LNS rule, run like: tCombine(MassIn, mygamma = 1)
    ## if want to use LNSa rule, run like: tCombine(MassIn, approximate = TRUE)
 
	
    ## mygamma is the parameter of the family of conjunctive and disjunctive rules using triangular norms by Denoeux.
	## mygamma = 1, with ifnormalize = FALSE, smets conjunctive rule
	## mygamma = 1, with ifnormalize = TRUE, Dempster rule
	## mygamma = 0, cautious rule
	## mygamma between 0 and 1, the generalized case of the rules using triangular by Denoeux
	## mygamma = 1, ifnormalize = FALSE, ifdiscount = TRUE, LNS rule
	## approximate = TRUE, LNSa rule, the approximation method for LNS rule
	## eta, the parameter in LNS rule, control the specificity of the decision

    nf = nrow(MassIn)
	n = ncol(MassIn)
	ThetaSize = log2(nf)
	w_mat = apply(MassIn, 2, mtow)
    if(approximate){

		  num_eff = apply(w_mat, 1, function(x){sum(abs(x - 1) > 1e-6)}) 
		  id_eff = which(num_eff > 0)
		  num_group_eff = num_eff[id_eff];
		  beta_vec = rep(1, length(id_eff));
		  if(eta != 0){
   	   		myc = sapply(1:nf -1, function(xx){
								sum(dec2bin(xx, ThetaSize))
						}) 
			 beta_vec = (ThetaSize/myc[id_eff])^eta
		  }
		  alpha_vec = beta_vec * num_group_eff / sum(beta_vec * num_group_eff)
		  w_eff =  1 - alpha_vec 
		  w_vec = rep(1, nf)
		  w_vec[id_eff] = w_eff 
        
	}else{ 
		if(mygamma == 1){
		  ## conjunctive rule or dempster rule
		  w_vec = apply(w_mat, 1, prod)	
		}else if(mygamma == 0){
		  w_vec = apply(w_mat, 1, min) 
		}else{
		  ## I donot know if it is right for more than 2 masses
		  w_vec = apply(w_mat, 1, function(x){prod(x)/prod(max(c(x, mygamma))[1:(length(x)-1)])})	
		}

		if(ifdiscount){
		  ## find the masses that are not total ignorance
		  num_eff = apply(w_mat, 1, function(x){sum(abs(x - 1) > 1e-6)}) 
		  id_eff = which(num_eff > 0)
		  w_eff = w_vec[id_eff];
		  num_group_eff = num_eff[id_eff];
		  beta_vec = rep(1, length(id_eff));
		  if(eta != 0){
   			myc = sapply(1:nf -1, function(xx){
								sum(dec2bin(xx, ThetaSize))
						}) 
			 beta_vec = (ThetaSize/myc[id_eff])^eta
		  }
		  alpha_vec = beta_vec * num_group_eff / sum(beta_vec * num_group_eff)
		  w_eff =  1 - alpha_vec + alpha_vec * w_eff 
		  w_vec[id_eff] = w_eff 
		}
	}

    out = wtom(w_vec)
	if(ifnormalize && mygamma == 1){
	   ## dempster rule
	   out[1] = 0;
	   out = out/sum(out)
	}
	return(t(out))
}



tCombine_SSF <- function(MassIn, mygamma, ifnormalize = FALSE, ifdiscount = TRUE, approximate = FALSE, eta = 0, singleton = FALSE){

	## Problem find: for conjunctive rule (and also for DS rule), when the number of masses if quite large, 
	##             more than one focal elements will have mass one
    
    ## LNS rule.  Also can be used for conjunctive rule, cautious rule and DS rule

	## only for SSF and non-dogmatic masses
	## singleton = TRUE, all the focal elements are singletons
	## singleton = FALSE, any kind of focal elements, including singletons certainly

	## useful for the EKNN and BeliefKNN, when singleton = TRUE

 
    ## MassIn: if singleton = TRUE, it is a matrix of (ThetaSize+1) * n. Each column is a bba. 
	##           The first ThetaSize rows are for masses on the singleton, the last row is for total ignorance, Theta
	##         if singleton = FALSE, it is a matrix of (2^ThetaSize) * n. Each column is a bba

	## we note that, when all the focal elements are singletons, we can also use singleton = FALSE to get the same results


    ## mygamma is the parameter of the family of conjunctive and disjunctive rules using triangular norms by Denoeux.
	## mygamma = 1, with ifnormalize = FALSE, smets conjunctive rule
	## mygamma = 1, with ifnormalize = TRUE, Dempster rule
	## mygamma = 0, cautious rule
	## mygamma between 0 and 1, the generalized case of the rules using triangular by Denoeux
	## mygamma = 1, ifnormalize = FALSE, ifdiscount = TRUE, LNS rule
	## approximate = TRUE, LNSa rule, the approximation method for LNS rule
	## eta, the parameter in LNS rule, control the specificity of the decision, only singleton = TRUE is useful

    if(singleton){
		ThetaSize = nrow(MassIn) - 1;
		nf = 2^ThetaSize
		n = ncol(MassIn)
		w_mat = MassIn[1:ThetaSize, ];
		w_mat = 1 - w_mat
		eta = 0
    }else{
		nf = nrow(MassIn);
		ThetaSize = log2(nf);
		w_mat = MassIn[1: (nf - 1), ]
		w_mat = 1 - w_mat
	}
    if(approximate){

		  num_eff = apply(w_mat, 1, function(x){sum(abs(x - 1) > 1e-6)}) 
		  id_eff = which(num_eff > 0)
		  num_group_eff = num_eff[id_eff];
		  if(eta != 0){
			  beta_vec = rep(1, length(id_eff));
			  myc = sapply(1:nf -1, function(xx){
									sum(dec2bin(xx, ThetaSize))
							}) 
			  beta_vec = (ThetaSize/myc[id_eff])^eta
		      alpha_vec = beta_vec * num_group_eff / sum(beta_vec * num_group_eff)
		  }else{
		      alpha_vec =  num_group_eff / sum(num_group_eff)
		  }
		  w_eff =  1 - alpha_vec 
		  if(singleton){
		    w_vec = rep(1, ThetaSize)
		  }else{
		    w_vec = rep(1, nf - 1)
		  }
		  w_vec[id_eff] = w_eff 
	}else{ 
		if(mygamma == 1){
		  ## conjunctive rule or dempster rule
		  w_vec = apply(w_mat, 1, prod)	
		}else if(mygamma == 0){
		  w_vec = apply(w_mat, 1, min) 
		}else{
		  ## I donot know if it is right for more than 2 masses
		  w_vec = apply(w_mat, 1, function(x){prod(x)/prod(max(c(x, mygamma))[1:(length(x)-1)])})	
		}

		if(ifdiscount){
		  ## find the masses that are not total ignorance
		  num_eff = apply(w_mat, 1, function(x){sum(abs(x - 1) > 1e-6)}) 
		  id_eff = which(num_eff > 0)
		  w_eff = w_vec[id_eff];
		  num_group_eff = num_eff[id_eff];
		  if(eta != 0){
			  beta_vec = rep(1, length(id_eff));
			  myc = sapply(1:nf -1, function(xx){
									sum(dec2bin(xx, ThetaSize))
							}) 
			  beta_vec = (ThetaSize/myc[id_eff])^eta
		      alpha_vec = beta_vec * num_group_eff / sum(beta_vec * num_group_eff)
		  }else{
		     alpha_vec =  num_group_eff / sum(num_group_eff)
		  }
		  w_eff =  1 - alpha_vec + alpha_vec * w_eff 
		  w_vec[id_eff] = w_eff 
		}
	}
	w_vec_complete = rep(1, nf);
    if(singleton){
		w_vec_complete[2^(1:ThetaSize-1) + 1] = w_vec
	}else{
		w_vec_complete[1: (nf - 1)] = w_vec
	}
	if(min(w_vec_complete)>0){
	   out = wtom(w_vec_complete)
	}else{
	   id = which(w_vec_complete == 0)
       out = rep(0, nf) 	
	   out[id] = 1
	}
	if(ifnormalize && mygamma == 1){
	   ## dempster rule
	   out[1] = 0;
	   out = out/sum(out)
	}
	return(t(out))
}





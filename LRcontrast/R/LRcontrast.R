#SIMULATING QUANTILES:

#define variables:

#not changeable for the user:

env1<-new.env()													#environment in which scal, off, i is saved
SCALE1<-c(1,1,1,1)												#scaling of the parameters in fct optim (parscale)
LOW1<-c(0,0,0,0)												#lower bounds for (phi1, phi2, theta1, theta2)
UP1<-c(pi,pi/2,6,6)												#upper bounds for (phi1, phi2, theta) (LOW and UP needed in nlminb)
INTERVAL<-c(0, pi)												#interval for phi (needed in optimize)


#starting values for nlminb in qLRcontrast and pLRcontrast:
#for the models linear, linlog, quadratic, emax, exponential:

initVal1<-matrix(nrow=3,ncol=3)									#columns of the matrix are different starting values for (phi1, phi2, theta) in nlminb
initVal1[3,]<-rep(2,3)
initVal1[2,]<-rep(1,3)
initVal1[1,]<-c(0,pi/2,pi)

#for the models logistic and betaMod

initVal2<-matrix(nrow=4,ncol=8)	
initVal2[4,]<-rep(2,8)
initVal2[3,]<-rep(2,8)
initVal2[2,]<-rep(1,8)
initVal2[1,]<-c(0,0.5,1,1.5,2,2.5,3,pi)

#for the model sigEmax

initVal3<-matrix(nrow=4,ncol=4)			
initVal3[4,]<-rep(2,4)
initVal3[3,]<-rep(2,4)
initVal3[2,]<-rep(1,4)
initVal3[1,]<-c(0,1,2,pi)


#define the stochastic process:
#factor (-1), because nlminb minimizes functions instead of maximizing ( max{f(x)} = -min{-f(x)} )):
#parametrization with polar coordinates:

#here: beta0 = cos(phi), sigma = sin(phi)

WS1<-function(phi){return(sum(sqrt(env1$weight)*(cos(phi)*env1$Z[,env1$i]+sqrt(2)*sin(phi)*env1$Y[,env1$i]))/sqrt(sum(env1$weight*(cos(phi)^2+2*sin(phi)^2))))}	


#here: beta0 = cos(phi1)*cos(phi2), beta1 = sin(phi2), sigma = sin(phi1)*cos(phi2)

WS2_linear<-function(x){return(.C("WS2linear", phi1=as.double(x[1]), phi2=as.double(x[2]), weight=as.double(env1$weight), dose=as.double(env1$dose), Z=as.double(env1$Z[,env1$i]), Y=as.double(env1$Y[,env1$i]), doselength = as.integer(length(env1$dose)))$phi1)}	
WS2_emax<-function(x){return(.C("WS2emax", phi1=as.double(x[1]), phi2=as.double(x[2]), ed50=as.double(x[3]), weight=as.double(env1$weight), dose=as.double(env1$dose), Z=as.double(env1$Z[,env1$i]), Y=as.double(env1$Y[,env1$i]), doselength = as.integer(length(env1$dose)))$phi1)}
WS2_exponential<-function(x){return(.C("WS2exponential", phi1=as.double(x[1]), phi2=as.double(x[2]), delta=as.double(x[3]), weight=as.double(env1$weight), dose=as.double(env1$dose), Z=as.double(env1$Z[,env1$i]), Y=as.double(env1$Y[,env1$i]), doselength = as.integer(length(env1$dose)))$phi1)}
WS2_linlog<-function(x){return(.C("WS2linlog", phi1=as.double(x[1]), phi2=as.double(x[2]), off=as.double(env1$off), weight=as.double(env1$weight), dose=as.double(env1$dose), Z=as.double(env1$Z[,env1$i]), Y=as.double(env1$Y[,env1$i]), doselength = as.integer(length(env1$dose)))$phi1)}
WS2_sigEmax<-function(x){return(.C("WS2sigEmax", phi1=as.double(x[1]), phi2=as.double(x[2]), ed50=as.double(x[3]), h=as.double(x[4]), weight=as.double(env1$weight), dose=as.double(env1$dose), Z=as.double(env1$Z[,env1$i]), Y=as.double(env1$Y[,env1$i]), doselength = as.integer(length(env1$dose)))$phi1)}
WS2_quadratic<-function(x){return(.C("WS2quadratic", phi1=as.double(x[1]), phi2=as.double(x[2]), b=as.double(x[3]), weight=as.double(env1$weight), dose=as.double(env1$dose), Z=as.double(env1$Z[,env1$i]), Y=as.double(env1$Y[,env1$i]), doselength = as.integer(length(env1$dose)))$phi1)}
WS2_betaMod<-function(x){return(.C("WS2betaMod", phi1=as.double(x[1]), phi2=as.double(x[2]), delta1=as.double(x[3]), delta2=as.double(x[4]), scal=as.double(env1$scal), weight=as.double(env1$weight), dose=as.double(env1$dose), Z=as.double(env1$Z[,env1$i]), Y=as.double(env1$Y[,env1$i]), doselength = as.integer(length(env1$dose)))$phi1)}
WS2_logistic<-function(x){return(.C("WS2logistic", phi1=as.double(x[1]), phi2=as.double(x[2]), ed50=as.double(x[3]), delta=as.double(x[4]), weight=as.double(env1$weight), dose=as.double(env1$dose), Z=as.double(env1$Z[,env1$i]), Y=as.double(env1$Y[,env1$i]), doselength = as.integer(length(env1$dose)))$phi1)}




#maximizing the functions:

#maximize WS1 over INTERVAL=[0,2*pi] (M1) (equal for all models)

M1max<-function(x){return(max(optimize(WS1, interval=INTERVAL, maximum=TRUE)$objective, 0)^2)} 


#maximize WS2 for the different models over M2. Maximize with nlminb with different starting values and take the maximum.  

MAX_linear<-function(x){max(max(0,(-1)*nlminb(initVal1[1:2,1], WS2_linear, lower=LOW1[1:2], upper=UP1[1:2], control=list(iter.max=5000,eval.max=5000))$objective)^2,
							max(0,(-1)*nlminb(initVal1[1:2,2], WS2_linear, lower=LOW1[1:2], upper=UP1[1:2], control=list(iter.max=5000,eval.max=5000))$objective)^2,
							max(0,(-1)*nlminb(initVal1[1:2,3], WS2_linear, lower=LOW1[1:2], upper=UP1[1:2], control=list(iter.max=5000,eval.max=5000))$objective)^2)}

MAX_emax<-function(x){max(max(0,(-1)*nlminb(initVal1[,1], WS2_emax, lower=LOW1[1:3], upper=UP1[1:3], control=list(iter.max=5000,eval.max=5000))$objective)^2,
						  max(0,(-1)*nlminb(initVal1[,2], WS2_emax, lower=LOW1[1:3], upper=UP1[1:3], control=list(iter.max=5000,eval.max=5000))$objective)^2,
						  max(0,(-1)*nlminb(initVal1[,3], WS2_emax, lower=LOW1[1:3], upper=UP1[1:3], control=list(iter.max=5000,eval.max=5000))$objective)^2)}

MAX_exponential<-function(x){max(max(0,(-1)*nlminb(initVal1[,1], WS2_exponential, lower=LOW1[1:3], upper=UP1[1:3], control=list(iter.max=5000,eval.max=5000))$objective)^2,
								 max(0,(-1)*nlminb(initVal1[,2], WS2_exponential, lower=LOW1[1:3], upper=UP1[1:3], control=list(iter.max=5000,eval.max=5000))$objective)^2,
								 max(0,(-1)*nlminb(initVal1[,3], WS2_exponential, lower=LOW1[1:3], upper=UP1[1:3], control=list(iter.max=5000,eval.max=5000))$objective)^2)}

MAX_linlog<-function(x){max(max(0,(-1)*nlminb(initVal1[1:2,1], WS2_linlog, lower=LOW1[1:2], upper=UP1[1:2], control=list(iter.max=5000,eval.max=5000))$objective)^2,
							max(0,(-1)*nlminb(initVal1[1:2,2], WS2_linlog, lower=LOW1[1:2], upper=UP1[1:2], control=list(iter.max=5000,eval.max=5000))$objective)^2,
							max(0,(-1)*nlminb(initVal1[1:2,3], WS2_linlog, lower=LOW1[1:2], upper=UP1[1:2], control=list(iter.max=5000,eval.max=5000))$objective)^2)}


MAX_quadratic<-function(x){max(max(0,(-1)*nlminb(initVal1[,1], WS2_quadratic, lower=LOW1[1:3], upper=UP1[1:3], control=list(iter.max=5000,eval.max=5000))$objective)^2,
							   max(0,(-1)*nlminb(initVal1[,2], WS2_quadratic, lower=LOW1[1:3], upper=UP1[1:3], control=list(iter.max=5000,eval.max=5000))$objective)^2,
							   max(0,(-1)*nlminb(initVal1[,3], WS2_quadratic, lower=LOW1[1:3], upper=UP1[1:3], control=list(iter.max=5000,eval.max=5000))$objective)^2)}

MAX_betaMod<-function(x){max(max(0,(-1)*nlminb(initVal2[,1], WS2_betaMod, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							 max(0,(-1)*nlminb(initVal2[,2], WS2_betaMod, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							 max(0,(-1)*nlminb(initVal2[,3], WS2_betaMod, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							 max(0,(-1)*nlminb(initVal2[,4], WS2_betaMod, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							 max(0,(-1)*nlminb(initVal2[,5], WS2_betaMod, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							 max(0,(-1)*nlminb(initVal2[,6], WS2_betaMod, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							 max(0,(-1)*nlminb(initVal2[,7], WS2_betaMod, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							 max(0,(-1)*nlminb(initVal2[,8], WS2_betaMod, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2)}

MAX_logistic<-function(x){max(max(0,(-1)*nlminb(initVal2[,1], WS2_logistic, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							  max(0,(-1)*nlminb(initVal2[,2], WS2_logistic, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							  max(0,(-1)*nlminb(initVal2[,3], WS2_logistic, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							  max(0,(-1)*nlminb(initVal2[,4], WS2_logistic, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							  max(0,(-1)*nlminb(initVal2[,5], WS2_logistic, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							  max(0,(-1)*nlminb(initVal2[,6], WS2_logistic, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							  max(0,(-1)*nlminb(initVal2[,7], WS2_logistic, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							  max(0,(-1)*nlminb(initVal2[,8], WS2_logistic, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2)}

MAX_sigEmax<-function(x){max(max(0,(-1)*nlminb(initVal3[,1], WS2_sigEmax, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							 max(0,(-1)*nlminb(initVal3[,2], WS2_sigEmax, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							 max(0,(-1)*nlminb(initVal3[,3], WS2_sigEmax, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2,
							 max(0,(-1)*nlminb(initVal3[,4], WS2_sigEmax, lower=LOW1, upper=UP1, control=list(iter.max=5000,eval.max=5000))$objective)^2)}
							  
							 

#with that function the user can simulate the quantiles:

qLRcontrast<-function(dose, probs, models, weight = rep(1/length(dose), length(dose)), off = 0.01 * max(dose), scal = 1.2 * max(dose), nsim = 10000, info = TRUE)
             {
              env1$weight<-weight;
              if(sum(env1$weight)!=1){stop("The sum of 'weight' needs to be 1")};
              if(min(env1$weight)<=0){stop("Only positive values in 'weight' allowed")};
              env1$dose<-dose;
              if(min(env1$dose) < 0){stop("Only dose-levels >= 0 allowed")};
              if(length(env1$dose)!=length(env1$weight)){stop("Lengths of the vectors 'dose' and 'weight' need to be equal")};
              env1$probs<-sort(probs); 
              if(env1$probs[1]<0 | 1<max(env1$probs)){stop("'probs' needs to be a numeric vector of probabilities with values in [0,1]")};
              env1$off<-off;
              if(env1$off <= 0){stop("'off' parameter needs to be positive")};
              env1$scal<-scal;
              if(env1$scal < max(dose)){stop("'scal' parameter needs to be >= max(dose)")};
              if(length(models)!=length(unique(models))){stop("Only one list entry allowed for each model")};
              if(nsim<=0 | nsim%%1!=0){stop("'nsim' needs to be a positive integer")};
              env1$Z<-matrix(data=rnorm(length(dose)*nsim),nrow=length(dose));
              env1$Y<-matrix(data=rnorm(length(dose)*nsim),nrow=length(dose));
              env1$L<-matrix(nrow=length(models),ncol=nsim);
              L_max<-c();
              QUANTILE<-matrix(nrow=length(models)+1, ncol=length(env1$probs)); colnames(QUANTILE)<-env1$probs; rownames(QUANTILE)<-c(models, "max");
              maxFkt<-list()

              for(j in 1:length(models))
              {
               switch(models[j],
               linear={maxFkt[[j]]<-MAX_linear},
               emax={maxFkt[[j]]<-MAX_emax},
               exponential={maxFkt[[j]]<-MAX_exponential},
               linlog={maxFkt[[j]]<-MAX_linlog},
               sigEmax={maxFkt[[j]]<-MAX_sigEmax},
               quadratic={maxFkt[[j]]<-MAX_quadratic},
               betaMod={maxFkt[[j]]<-MAX_betaMod},
               logistic={maxFkt[[j]]<-MAX_logistic},
               stop("The vector 'models' contains an incorrect entry")
               )
               };
			   
			  start<-Sys.time();
			  if(info==TRUE){cat(paste("% done    \t", "nsim done\t", "secs passed\t", "secs remaining ~", collapse=""), sep="", fill=TRUE)};
              for(j in 1:nsim){env1$i<-j; M1_max<-M1max(1);
                                for(k in 1:length(models)){env1$L[k,j]<-maxFkt[[k]](1)-M1_max};
                                L_max[j]<-max(env1$L[,j])
								
								if(info==TRUE){
								passed <- round(difftime(Sys.time(),start,units="secs")[[1]],0);
								percent<- round(j/nsim*100,2);
								if(j%%10==0){env1$remaining<- round(passed/j*(nsim-j),0)}; 
								cat(paste(rep.int("\b", 60),collapse=""),sep="");
								cat(paste(percent, "\t\t", j, "\t\t", passed, "\t\t", env1$remaining, collapse=""), sep="");
								flush.console()}
								};
              for(j in 1:length(models)){QUANTILE[j,]<-quantile(env1$L[j,],env1$probs)};
              QUANTILE[length(models)+1,]<-quantile(L_max,env1$probs);

              cat("\n");
			  return(QUANTILE)}
			  
			  

#test statistic (2.9):


sLRcontrast<-function(dose, resp, models, off = 0.01 * max(dose), scal= 1.2 * max(dose))
			{
			if(length(dose)!=length(resp)){stop("The length of the vectors 'dose' and 'resp' needs to be equal")};
			if(min(dose) < 0){stop("Only dose-levels >= 0 allowed")};
			if(off <= 0){stop("'off' parameter needs to be positive")};
			if(scal < max(dose)){stop("'scal' parameter needs to be >= max(dose)")};
			if(length(models)!=length(unique(models))){stop("Only one list entry allowed for each model")};
			
			result <- matrix(nrow=length(models)+1, ncol= 1); rownames(result)<-c(models, "max"); colnames(result)<-c("statistic");
			
			RSS <- c();
			addArgs <- list(); addArgs$off <- off; addArgs$scal <- scal;
			bnds <- DoseFinding::defBnds(mD = max(dose));
			
			minT0<-sum((resp-1/length(resp)*sum(resp))^2);
			
			for(j in 1:length(models))
			{
			switch(models[j],
			linear = {fit <- DoseFinding::fitMod(dose, resp, model = "linear", addArgs = addArgs, bnds = bnds$linear); 
			          if(fit$coefs[[2]] > 0){RSS[j] <- fit$RSS} else{RSS[j] <- minT0}},
			
			emax = {fit <- DoseFinding::fitMod(dose, resp, model = "emax", addArgs = addArgs, bnds = bnds$emax); 
			        if(fit$coefs[[2]] > 0){RSS[j] <- fit$RSS} else{RSS[j] <- minT0}},

			exponential = {fit <- DoseFinding::fitMod(dose, resp, model = "exponential", addArgs = addArgs, bnds = bnds$exponential); 
			               if(fit$coefs[[2]] > 0){RSS[j] <- fit$RSS} else{RSS[j] <- minT0}},

			linlog = {fit <- DoseFinding::fitMod(dose, resp, model = "linlog", addArgs = addArgs, bnds = bnds$linlog); 
			          if(fit$coefs[[2]] > 0){RSS[j] <- fit$RSS} else{RSS[j] <- minT0}},
					  
			sigEmax = {fit <- DoseFinding::fitMod(dose, resp, model = "sigEmax", addArgs = addArgs, bnds = bnds$sigEmax); 
			           if(fit$coefs[[2]] > 0){RSS[j] <- fit$RSS} else{RSS[j] <- minT0}},
					  
			quadratic = {fit <- DoseFinding::fitMod(dose, resp, model = "quadratic", addArgs = addArgs, bnds = bnds$quadratic); 
			             if(fit$coefs[[2]] > 0){RSS[j] <- fit$RSS} else{RSS[j] <- minT0}},
					  
			betaMod = {fit <- DoseFinding::fitMod(dose, resp, model = "betaMod", addArgs = addArgs, bnds = bnds$betaMod); 
			           if(fit$coefs[[2]] > 0){RSS[j] <- fit$RSS} else{RSS[j] <- minT0}},
					   
			logistic = {fit <- DoseFinding::fitMod(dose, resp, model = "logistic", addArgs = addArgs, bnds = bnds$logistic); 
			            if(fit$coefs[[2]] > 0){RSS[j] <- fit$RSS} else{RSS[j] <- minT0}},
						
			stop("The vector 'models' contains an incorrect entry")
			)
			};
			
			
			for(j in 1:length(models)){result[j,1]<-round(length(resp)*log(minT0/RSS[j]),4)};
			result[length(models)+1,1]<-max(result[1:length(models),1]);
			return(result)
			}



			
			
# SIMULATING the p-value

pLRcontrast<-function(dose, resp, models, off = 0.01 * max(dose), scal = 1.2 * max(dose), nsim = 1000, info = TRUE)
             {
			  env1$data<-resp[order(dose)];
			  if(length(dose)!=length(resp)){stop("The length of the vectors 'dose' and 'resp' needs to be equal")};
			  env1$dose<-sort(unique(dose));
			  env1$weight<-c(); for(i in 1:length(env1$dose)){env1$weight[i]<-sum(env1$dose[i]==dose)/length(dose)};
			  env1$info=info;
              if(min(env1$dose) < 0){stop("Only dose-levels >= 0 allowed")};
              env1$off<-off;
              if(env1$off <= 0){stop("'off' parameter needs to be positive")};
              env1$scal<-scal;
              if(env1$scal < max(dose)){stop("'scal' parameter needs to be >= max(dose)")};
              if(length(models)!=length(unique(models))){stop("Only one list entry allowed for each model")};
              if(nsim<=0 | nsim%%1!=0){stop("'nsim' needs to be a positive integer")};

              env1$statistic<-sLRcontrast(dose=dose, resp=resp, models=models, off=env1$off, scal=env1$scal);
  
              env1$Z<-matrix(data=rnorm(length(env1$dose)*nsim),nrow=length(env1$dose));
              env1$Y<-matrix(data=rnorm(length(env1$dose)*nsim),nrow=length(env1$dose));
              env1$L<-matrix(nrow=length(models),ncol=nsim);
              env1$L_max<-c();
              maxFkt<-list()
              result<-matrix(nrow=length(models), ncol= 2); rownames(result)<-models; colnames(result)<-c("unadj-p", "adj-p");
              
              for(j in 1:length(models))
              {
               switch(models[j],
               linear={maxFkt[[j]]<-MAX_linear},
               emax={maxFkt[[j]]<-MAX_emax},
               exponential={maxFkt[[j]]<-MAX_exponential},
               linlog={maxFkt[[j]]<-MAX_linlog},
               sigEmax={maxFkt[[j]]<-MAX_sigEmax},
               quadratic={maxFkt[[j]]<-MAX_quadratic},
               betaMod={maxFkt[[j]]<-MAX_betaMod},
               logistic={maxFkt[[j]]<-MAX_logistic},
               stop("The vector 'models' contains an incorrect entry")
               )
               };
			  
			  start<-Sys.time();
			  if(info==TRUE){cat(paste("% done    \t", "nsim done\t", "secs passed\t", "secs remaining ~", collapse=""), sep="", fill=TRUE)};
              for(j in 1:nsim){env1$i<-j; M1_max<-M1max(1);
                                for(k in 1:length(models)){env1$L[k,j]<-max(maxFkt[[k]](1)-M1_max,0)};
                                env1$L_max[j]<-max(env1$L[,j]);
								
								if(info==TRUE){
								passed <- round(difftime(Sys.time(),start,units="secs")[[1]],0);
								percent<- round(j/nsim*100,2);
								if(j%%10==0){env1$remaining<- round(passed/j*(nsim-j),0)}; 
								cat(paste(rep.int("\b", 60),collapse=""),sep="");
								cat(paste(percent, "\t\t", j, "\t\t", passed, "\t\t", env1$remaining, collapse=""), sep="");
								flush.console()}
								};
              
              for(k in 1:length(models)){result[k,1]<-round(sum(env1$L[k,] >= env1$statistic[k])/nsim,4)};
			  for(k in 1:length(models)){result[k,2]<-round(sum(env1$L_max >= env1$statistic[k])/nsim,4)};
               

              cat("\n");
              return(result)}

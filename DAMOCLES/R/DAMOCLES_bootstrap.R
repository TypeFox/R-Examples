DAMOCLES_bootstrap = function(
   phy = rcoal(10),
   pa = matrix(c(phy$tip.label,sample(c(0,1),Ntip(phy),replace = T)),
      nrow = Ntip(phy),ncol = 2),
   initparsopt = c(0.1,0.1),
   idparsopt = 1:length(initparsopt),
   parsfix = 0,
   idparsfix = (1:3)[-idparsopt],
   pars2 = c(1E-3,1E-4,1E-5,1000),
   pchoice = 0,
   runs = 999,
   estimate_pars = FALSE,
   conf.int = 0.95)
{
				
	# CALCULATE OSBERVED N, MNTD AND MPD
	
	obs.samp = matrix(as.numeric(pa[,2]),nrow = 1)
	colnames(obs.samp) = phy$tip.label

	n.obs = sum(obs.samp)
	mntd.obs = picante::mntd(obs.samp, cophenetic(phy))
	mpd.obs = picante::mpd(obs.samp, cophenetic(phy))
	
	# ESTIMATE COLONISATION AND LOCAL EXTINCTION RATES FROM EMPIRICAL DATA		
	
	pars = DAMOCLES_ML(phy,pa,initparsopt = initparsopt,idparsopt = idparsopt,parsfix = parsfix,idparsfix = idparsfix,pars2 = pars2,pchoice = pchoice)
	
  # SIMULATE NULL COMMUNITY UNDER THESE ESTIMATED RATES	
  
  runs2 = runs
  if(pchoice == 0)
  {
     root.state = c(0,1)
     # MAKE NUMBER OF RUNS EVEN AND THEN DIVIDE BY HALF
     runs2 = ((runs + runs %% 2))/2
  }
  if(pchoice == 1)
  {
     root.state = 0
  }
  if(pchoice == 2)	
  {
     root.state = 1
  }
  
	# SIMULATE NULL COMMUNITIES UNDER RANDOM DRAW MODEL AND DAMOCLES USING ESTIMATED PARAMETERS AND CALCULATE PHYLOGENETIC STRUCTURE
	nullCommStats = expand.grid(run = 1:runs2,root.state = root.state,n = sum(as.numeric(pa[,2])),n.RD = sum(as.numeric(pa[,2])),mntd.RD = NA,mpd.RD = NA,n.DAMOCLES = NA,mntd.DAMOCLES = NA,mpd.DAMOCLES = NA,loglik.DAMOCLES = NA,mu.DAMOCLES = NA,gamma_0.DAMOCLES = NA)

	for(run in 1:dim(nullCommStats)[1]){     
    	DAMOCLES_community = DAMOCLES_sim(phy,mu = pars[1],gamma_0 = pars[2],gamma_td = 0,sigma = 0,psiBranch = 0,psiTrait = 0,z = 0,phi = 0,traitOpt = 0,br0 = 0,br_td = 0,nTdim = 1,root.state = nullCommStats$root.state[run],root.trait.state = 0,plotit = FALSE,keepExtinct = FALSE)	
  		DAMOCLES.samp = matrix(as.numeric(DAMOCLES_community[[1]]$state),nrow = 1)
  		colnames(DAMOCLES.samp) = phy$tip.label
  		nullCommStats$n.DAMOCLES[run] = sum(DAMOCLES.samp)
  		nullCommStats$mntd.DAMOCLES[run] = picante::mntd(DAMOCLES.samp,cophenetic(phy))
  		nullCommStats$mpd.DAMOCLES[run] = picante::mpd(DAMOCLES.samp,cophenetic(phy))
  		RD.samp = matrix(as.numeric(sample(pa[,2])),nrow = 1)
  		colnames(RD.samp) = phy$tip.label
  		nullCommStats$mntd.RD[run] = picante::mntd(RD.samp,cophenetic(phy))
  		nullCommStats$mpd.RD[run] = picante::mpd(RD.samp,cophenetic(phy))
  		
  		if(estimate_pars == TRUE)
      {
	  	   DAMOCLES.pa = pa
	  	   DAMOCLES.pa[,2]= DAMOCLES_community[[1]]$state	
  		   parsNull = DAMOCLES_ML(phy,DAMOCLES.pa,initparsopt = c(pars[[1]],pars[[2]],0),idparsopt = idparsopt,parsfix = parsfix,idparsfix = idparsfix,pars2 = pars2,pchoice = pchoice)
		     nullCommStats$mu.DAMOCLES[run] = parsNull[[1]]
	       nullCommStats$gamma_0.DAMOCLES[run] = parsNull[[2]]
		     nullCommStats$loglik.DAMOCLES[run] = parsNull[[4]]
      }
  	}
	# CALCULATE MEAN AND STANDARD DEVIATION OF N, MNTD AND MPD FOR SIMULATED DATA  	
	n.mean.DAMOCLES = mean(nullCommStats$n.DAMOCLES)
	mntd.mean.DAMOCLES = mean(nullCommStats$mntd.DAMOCLES)
	mntd.sd.DAMOCLES = sd(nullCommStats$mntd.DAMOCLES)
	mpd.mean.DAMOCLES = mean(nullCommStats$mpd.DAMOCLES)
	mpd.sd.DAMOCLES = sd(nullCommStats$mpd.DAMOCLES)
	loglik.mean.DAMOCLES = mean(nullCommStats$loglik.DAMOCLES)
	loglik.sd.DAMOCLES = sd(nullCommStats$loglik.DAMOCLES)
	
	mntd.mean.RD = mean(nullCommStats$mntd.RD)
	mntd.sd.RD = sd(nullCommStats$mntd.RD)
	mpd.mean.RD = mean(nullCommStats$mpd.RD)
	mpd.sd.RD = sd(nullCommStats$mpd.RD)
	
	# CALCULATE STANDARDISED EFFECTS SIZES OF OBSERVED DATA
	mntd.obs.z.RD = -1 * (mntd.obs - mntd.mean.RD)/mntd.sd.RD
	mntd.obs.rank.RD = rank(c(mntd.obs,nullCommStats$mntd.RD))[1]
	mntd.obs.q.RD = (mntd.obs.rank.RD/(dim(nullCommStats)[1] + 1))
	
	mpd.obs.z.RD = -1 * (mpd.obs - mpd.mean.RD)/mpd.sd.RD
	mpd.obs.rank.RD = rank(c(mpd.obs,nullCommStats$mpd.RD))[1]
	mpd.obs.q.RD = (mpd.obs.rank.RD/(dim(nullCommStats)[1] + 1))
		
	mntd.obs.z.DAMOCLES = -1 * (mntd.obs - mntd.mean.DAMOCLES)/mntd.sd.DAMOCLES
	mntd.obs.rank.DAMOCLES = rank(c(mntd.obs, nullCommStats$mntd.DAMOCLES))[1]
	mntd.obs.q.DAMOCLES = mntd.obs.rank.DAMOCLES/(dim(nullCommStats)[1] + 1)
	
	mpd.obs.z.DAMOCLES = -1 * (mpd.obs - mpd.mean.DAMOCLES)/mpd.sd.DAMOCLES
	mpd.obs.rank.DAMOCLES = rank(c(mpd.obs, nullCommStats$mpd.DAMOCLES))[1]
	mpd.obs.q.DAMOCLES = mpd.obs.rank.DAMOCLES/(dim(nullCommStats)[1] + 1)
	
	if(estimate_pars == TRUE)
  {
  	loglik.obs.z.DAMOCLES = (pars[[4]] - loglik.mean.DAMOCLES)/loglik.sd.DAMOCLES
  	loglik.obs.rank.DAMOCLES = rank(c(pars[[4]], nullCommStats$loglik.DAMOCLES))[1]
  	loglik.obs.q.DAMOCLES = loglik.obs.rank.DAMOCLES/(dim(nullCommStats)[1] + 1)
  	mu = paste(pars[[1]],"(",quantile(nullCommStats$mu.DAMOCLES,(1 - conf.int)/2),":",quantile(nullCommStats$mu.DAMOCLES,conf.int + (1 - conf.int)/2),")",sep = "")
  	gamma_0 = paste(pars[[2]],"(",quantile(nullCommStats$gamma_0.DAMOCLES,(1 - conf.int)/2),":",quantile(nullCommStats$gamma_0.DAMOCLES,conf.int + (1 - conf.int)/2),")",sep = "")
 	} else {
  	loglik.obs.z.DAMOCLES = NA
  	loglik.obs.rank.DAMOCLES = NA
  	loglik.obs.q.DAMOCLES = NA
  	mu = pars[[1]]
  	gamma_0 = pars[[2]]		
	}	
  if(pchoice == 0)	
  {
   	 root.state.print = paste(0,",",1,sep = "")
  } else {
     root.state.print = root.state
  }
	  	
	value = c(root.state.print,mu,gamma_0,pars[[4]],pars[[5]],pars[[6]],n.obs,mntd.obs,mpd.obs,dim(nullCommStats)[1],
	mntd.mean.RD,mntd.sd.RD,mntd.obs.z.RD,mntd.obs.rank.RD,mntd.obs.q.RD,
	mpd.mean.RD,mpd.sd.RD,mpd.obs.z.RD,mpd.obs.rank.RD,mpd.obs.q.RD,
	n.mean.DAMOCLES,
	mntd.mean.DAMOCLES,mntd.sd.DAMOCLES,mntd.obs.z.DAMOCLES,mntd.obs.rank.DAMOCLES,mntd.obs.q.DAMOCLES,
	mpd.mean.DAMOCLES,mpd.sd.DAMOCLES,mpd.obs.z.DAMOCLES,mpd.obs.rank.DAMOCLES,mpd.obs.q.DAMOCLES,
	loglik.mean.DAMOCLES,loglik.sd.DAMOCLES,loglik.obs.z.DAMOCLES,loglik.obs.rank.DAMOCLES,loglik.obs.q.DAMOCLES)
	
	desc = c("Root_state","mu","gamma_0","loglik.obs","df","conv","n.obs","mntd.obs","mpd.obs","runs",
	"mntd.mean.RD","mntd.sd.RD","mntd.obs.z.RD","mntd.obs.rank.RD","mntd.obs.q.RD",
	"mpd.mean.RD","mpd.sd.RD","mpd.obs.z.RD","mpd.obs.rank.RD","mpd.obs.q.RD",
	"n.mean.DAMOCLES",
	"mntd.mean.DAMOCLES","mntd.sd.DAMOCLES","mntd.obs.z.DAMOCLES","mntd.obs.rank.DAMOCLES","mntd.obs.q.DAMOCLES",
	"mpd.mean.DAMOCLES","mpd.sd.DAMOCLES","mpd.obs.z.DAMOCLES","mpd.obs.rank.DAMOCLES","mpd.obs.q.DAMOCLES",
	"loglik.mean.DAMOCLES","loglik.sd.DAMOCLES","loglik.obs.z.DAMOCLES","loglik.obs.rank.DAMOCLES","loglik.obs.q.DAMOCLES")
		
	summaryResults = data.frame(desc,value)
	
	# SAVE SUMMARY RESULTS AND ALSO RESULTS FOR EACH NULL COMMUNITY
	results = list()
	results[[1]] = summaryResults
	results[[2]] = nullCommStats
  names(results) = c("summary_table","null_community_data")
	
	return(results)
}	


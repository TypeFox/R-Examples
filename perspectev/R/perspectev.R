
perspectev.read <-
function(data,upper,lower,traits,survProb,traitfun,extinctionAge=NULL,occurrenceAge=NULL,vlist=NULL,trim=TRUE,projection=FALSE,quiet=FALSE){
	p = data
	
	cn = colnames(p)
	if(missing(survProb)){
		if(missing(occurrenceAge) || missing(extinctionAge)){stop("No supplied way of determining survivorship probability!")}
		if(!occurrenceAge %in% cn){stop(paste("Occurrence age column",occurrenceAge,"not in dataframe!"))}
	}else{
		if(!missing(occurrenceAge) || !missing(extinctionAge)){
			if(!quiet){print("Warning: occurrenceAge/extinctionAge specified but will be ignored.")}
		}
	}
	if(missing(traitfun) && trim==TRUE){stop("Trait function (traitfun) must be specified if trim=TRUE!")}
	if(missing(traits)){stop("Trait column names must be specified!")}
	if(!upper %in% cn){stop(paste("Upper level column",upper,"not in dataframe!"))}
	if(!lower %in% cn){stop(paste("Lower level column",lower,"not in dataframe!"))}
	for(t in traits){
		if(!t %in% cn){stop(paste("trait column",t,"not in dataframe!"))}
	}

	if(projection == TRUE){ #convert coordinates to mollweide projection
		proj = mapproject(p[,traits[2]],p[,traits[1]],projection="mollweide")
		p[,traits[2]] = proj$x * 100
		p[,traits[1]] = proj$y * 100
	} 

	p[,upper] = as.character(p[,upper])
	p[,lower] = as.character(p[,lower])

	before = NULL
	if(missing(survProb)){
		if(!quiet){print("Assigning survivorship = 1 to taxa with ages < extinctionAge (backward time)")}
		before = p[p[,occurrenceAge] >= extinctionAge,]
		after = p[p[,occurrenceAge] < extinctionAge,]
		beforesp = unique(paste(before[,upper],before[,lower]))
		aftersp = paste(after[,upper],after[,lower])
	}else{
		if(is.vector(survProb) && is.atomic(survProb)){
			if(length(survProb) == 0){stop("No survivors!")}
			if(!quiet){print("Assigning survivorship = 1 to taxa in survivorship vector")}
			before = p
			beforesp = unique(paste(before[,upper],before[,lower]))
			aftersp = survProb
			if(sum(!aftersp %in% beforesp) > 0){
				notina = aftersp[!aftersp %in% beforesp]
				print(notina)
				stop("Above are found in survProb but not in data set. All surviving taxa must be in data set.")
			}
		}
		if(class(survProb) == "data.frame" && ncol(survProb) == 2){
			if(nrow(survProb) == 0){stop("No survivors!")}
			if(!quiet){print("Assigning survivorship = 1 to taxa in survivorship dataframe")}
			before = p
			beforesp = unique(paste(before[,upper],before[,lower]))
			aftersp = paste(survProb[,1],survProb[,2])
			if(sum(!aftersp %in% beforesp) > 0){
				notina = aftersp[!aftersp %in% beforesp]
				print(notina)
				stop("Above are found in survProb but not in data set. All surviving taxa must be in data set.")
			}
		}
		if(class(survProb) == "data.frame" && ncol(survProb) == 3){
			if(class(survProb[,3]) != "numeric"){stop("Supplied survivorship probability must be numeric")}
			if(sum(survProb[,3] < 0) > 0 || sum(survProb[,3] > 1) > 0){stop("Supplied survivorship probability must be between 0 and 1")}
			if(nrow(survProb) == 0){stop("No survivors!")}
			if(!quiet){print("Assigning continuous survivorship probability by survivorship dataframe")}
			before = p
			beforesp = unique(paste(before[,upper],before[,lower]))
			aftersp = unique(paste(survProb[,1],survProb[,2]))
			if(sum(!aftersp %in% beforesp) > 0){
				notina = aftersp[!aftersp %in% beforesp]
				print(notina)
				stop("Above are found in survProb but not in data set. Need 1:1 mapping.")
			}
			if(sum(!beforesp %in% aftersp) > 0){
				notinb = beforesp[!beforesp %in% aftersp]
				print(notinb)
				stop("Above are found in data set but not in survProb. Need 1:1 mapping.")
			}
		}
	}
	if(length(before) == 0){stop("Survivorship probability input improperly formatted (perhaps not a dataframe?)")}
		
	
	survsp = unique(aftersp[aftersp %in% beforesp]) #units found across extinction age
	final = data.frame(matrix(NA,nrow=length(beforesp),ncol=3))
	colnames(final) = c("Upper","Lower","Survivorship")
	rownames(final) = beforesp
	beforenames = paste(before[,upper],before[,lower])
	before = cbind(before,beforenames)
	traitl = list() #list containing trait values
	for(spec in beforesp){
		temp = before[before$beforenames == spec,]
		traitl[[spec]] = list()
		for(tn in 1:length(traits)){	
			tname = paste0("T",tn)
			traitl[[spec]][tn] = list(temp[,traits[tn]])
		}
		final[spec,]$Lower = spec
		final[spec,]$Upper = temp[1,upper]
		final[spec,]$Survivorship = sum(spec %in% survsp)
	}

	if(!missing(survProb)){
	if(class(survProb) == "data.frame" && ncol(survProb) == 3){
		final$Survivorship == -1
		survspecs = paste(survProb[,1],survProb[,2]) 
		finalspecs = rownames(final)
		if(sum(!finalspecs %in% survspecs) > 0){
			print(finalspecs[!finalspecs %in% survspecs])
			stop("Above units not found in survivorship data frame!")
		}
		if(sum(!survspecs %in% finalspecs) > 0){
			print(survspecs[!survspecs %in% finalspecs])
			stop("Above units not found in occurrence data frame!")
		}
		final[paste(survProb[,1],survProb[,2]),]$Survivorship = survProb[,3]
		if(sum(final$Survivorship < 0) > 0){
			print(survspecs[!survspecs %in% finalspecs])
			stop("Above units not found in survivorship data frame!")
		}
	}
	}

	if(trim == TRUE){ #remove lower units that give NA for trait calculation
		nfinal = data.frame()
		ntraitl = list()
		for(i in 1:nrow(final)){
			traitdata = traitl[[final[i,]$Lower]]
			if(!is.na(traitfun(traitdata,vlist))){
				nfinal = rbind(nfinal,final[i,])
				ntraitl[[final[i,]$Lower]] = traitl[[final[i,]$Lower]]
			}
		}
		final = nfinal
		traitl = ntraitl
		if(nrow(final) == 0){stop("No taxa left after trimming!")}
		colnames(final) = c("Upper","Lower","Survivorship")
	}

	finall = list()
	finall$Survivorship = final
	finall$Traits = traitl
	finall$Key = data.frame(cbind(1:length(traits),traits))
	finall$Key[,1] = as.numeric(finall$Key[,1])
	colnames(finall$Key) = c("traitdata_index","traitdata_value")
	if(!quiet){print(finall$Key)}
	return(finall)
}



perspectev.calc <-
function(data,traitfun,vlist=list(),na.rm=FALSE){ 
	low = data.frame()
	lowns = c()
	for(i in 1:nrow(data$Survivorship)){
		td = data$Traits[[data$Survivorship[i,]$Lower]] 
		value = traitfun(td,vlist)
		if(!na.rm){
			low = rbind(low,c(data$Survivorship[i,]$Survivorship,value))
			lowns = c(lowns,data$Survivorship[i,]$Lower)
		}else{
			if(!is.na(value)){
				low = rbind(low,c(data$Survivorship[i,]$Survivorship,value))
				lowns = c(lowns,data$Survivorship[i,]$Lower)
			}
		}
	}
	up = data.frame()
	uppers = unique(data$Survivorship$Upper)
	upperns = c()
	for(g in uppers){
		temp = data$Survivorship[data$Survivorship$Upper == g,]
		ul = data$Traits[[temp[1,]$Lower]]
		if(nrow(temp) > 1){
			for(ln in 2:nrow(temp)){
				for(tln in 1:length(ul)){
					ul[[tln]] = c(ul[[tln]],data$Traits[[temp[ln,]$Lower]][[tln]])
				}
			}
		}	
		value = traitfun(ul,vlist)
		survivorship = 1-prod(1-temp$Survivorship)#probability that any species survives
		if(!na.rm){
			up = rbind(up,c(survivorship,value))
			upperns = c(upperns,g)
		}else{
			if(!is.na(value)){
				up = rbind(up,c(survivorship,value))
				upperns = c(upperns,g)
			}
		}
	}
	q = quantile(up[,2])[2:4] 		 #get quantiles 
	up[,2] = scale(up[,2]) 			 #scale value
	if(var(up[,1]) == 0 || var(up[,2]) == 0){
		correlation = NA
	}else{
		correlation = cor(up[,1],up[,2]) #calculate correlation
	}
	if(nrow(low)!=0){
		rownames(low) = lowns
		colnames(low) = c("Survivorship","Trait")
		low[,2] = scale(low[,2])
	}
	rownames(up) = upperns
	colnames(up) = c("Survivorship","Trait")
	total = list()
	total$lower = low
	total$upper = up
	total$stats = c(correlation,q)
	return(total)
}


perspectev.permute <-
function(data,permutations,traitfun,vlist,na.rm){
	simulations = c(); count = 0; fcount = 0
	while(count < permutations){
		temp = data 		
		allg = temp$Survivorship$Upper
		temp$Survivorship$Upper = allg[sample(1:length(allg),replace=FALSE)] 
		p = perspectev.calc(temp,traitfun,vlist,na.rm)
		if(!is.na(p$stats[1])){
			simulations = rbind(simulations,p$stats)
			count = count + 1
		}
		fcount = fcount + 1
		if(fcount >  100 & (fcount-count)/fcount > 0.9){
			stop('>90% of iterations yield either no victims or no survivors. There is likely too high a ratio of survivors or victims in the data.')
		}
	}
	rownames(simulations)=NULL
	return(simulations)
}

perspectev.test <-
function(data,iterations=1000,cores=1,traitfun=mcpRange,vlist=NULL,na.rm=FALSE){
	sim = c(); 	i = NULL
	real = perspectev.calc(data,traitfun,vlist,na.rm)
	if(cores > 1){ 
		div = floor(iterations/cores)
		plist = rep(div,length=cores)
		plist[length(plist)] = plist[length(plist)] + iterations %% cores
		cl = makeCluster(cores)
		registerDoParallel(cl)
		sim = foreach(i = 1:cores,.packages=c('perspectev'),.combine='rbind') %dopar% perspectev.permute(data,plist[i],traitfun,vlist,na.rm)
		stopCluster(cl)
	}else{
		sim = perspectev.permute(data,iterations,traitfun,vlist,na.rm)
	}
	total = list() 
	total$correlation_permuted = sim[,1]
	total$correlation_observed= real$stats[1]
	total$pvalue = sum(sim[,1] >= real$stats[1])/nrow(sim)
	permutedqs = cbind(sim[,2],sim[,3],sim[,4])
	colnames(permutedqs) = c("25%","50%","75%")
	total$permuted_quantiles = permutedqs
	return(total)
}


perspectev.simulation.permute <-
function(data,simulations,intercept,slope,level,traitfun,vlist,binary,noise){
	sim = c(); count = 0; fcount = 0
	pc = perspectev.calc(data,traitfun,vlist)
	if(level=="lower"){
		while(count < simulations){
			temp = data
			survProb = inv.logit(slope * pc$lower[,2] + rnorm(nrow(temp$Survivorship), sd=noise) + intercept)
			if(binary==TRUE){
				temp$Survivorship$Survivorship = as.numeric(survProb >= runif(nrow(temp$Survivorship)))
			 	while(sum(temp$Survivorship$Survivorship) == 0 || sum(temp$Survivorship$Survivorship) == nrow(temp$Survivorship)){
			 		temp$Survivorship$Survivorship = as.numeric(survProb >= runif(nrow(temp$Survivorship)))}
			}else{temp$Survivorship$Survivorship = survProb}
			p = perspectev.test(temp,1,1,traitfun,vlist) 
			if(!is.na(p$correlation_permuted) & !is.na(p$correlation_observed) ){
				sim = rbind(sim,c(p$correlation_permuted,p$correlation_observed))
				count = count + 1}
			fcount = fcount + 1
			if(fcount >  100 & (fcount-count)/fcount > 0.9){
				stop('>90% of simulations yield either no victims or no survivors. There is likely too high a ratio of survivors or victims in the data.')}
		}
	}else{
		while(count < simulations){
			temp = data; addOrder = c()
			survProb = inv.logit(slope * pc$upper[,2] + rnorm(nrow(pc$upper), sd=noise) + intercept)
			names(survProb) = rownames(pc$upper)
			Espec = sum(data$Survivorship$Survivorship)
			temp$Survivorship$Survivorship = 0
			uppernames = unique(temp$Survivorship$Upper)
			for(u in uppernames){
				st = data$Survivorship[data$Survivorship$Upper == u,]
				if(sum(st$Survivorship) > 0){
					surv = sample(st$Lower,size=1,prob=st$Survivorship)
				}else{
					surv = sample(st$Lower,size=1)
				}
				temp$Survivorship[surv,]$Survivorship = survProb[u]
			}
			Ospec = sum(temp$Survivorship$Survivorship)
			leftnames = temp$Survivorship[temp$Survivorship$Survivorship == 0,]$Lower
			leftprob = data$Survivorship[temp$Survivorship$Survivorship == 0,]$Survivorship
			if(length(leftnames[leftprob > 0]) > 0){
				addOrder = sample(leftnames[leftprob > 0],size=length(leftnames[leftprob > 0]),prob=leftprob[leftprob > 0])
				addOrder = c(addOrder,sample(leftnames[leftprob == 0]))
			}else{addOrder=sample(leftnames)}
			counter = 1
			while(Ospec < Espec && counter <= length(addOrder)){
				spsurv = addOrder[counter]
				supper = temp$Survivorship[spsurv,]$Upper
				utemp = temp$Survivorship[temp$Survivorship$Upper == supper & temp$Survivorship$Survivorship != 0,]#get existing non-zero species in upper level
				nprob = 1 - (1-survProb[supper])^(1/(nrow(utemp)+1)) #recalculate survivorship probability 
				temp$Survivorship[utemp$Lower,]$Survivorship = nprob
				temp$Survivorship[spsurv,]$Survivorship = nprob
				counter = counter + 1
				Ospec = sum(temp$Survivorship$Survivorship)
			}
			if(binary == TRUE){
				sProb = temp$Survivorship$Survivorship
				temp$Survivorship$Survivorship = as.numeric(sProb >= runif(nrow(temp$Survivorship)))
			 	while(sum(temp$Survivorship$Survivorship) == 0 || sum(temp$Survivorship$Survivorship) == nrow(temp$Survivorship)){
			 		temp$Survivorship$Survivorship = as.numeric(sProb >= runif(nrow(temp$Survivorship)))}
			}
			p = perspectev.test(temp,1,1,traitfun,vlist) #occasionally get NA values from all or no upper levels surviving
			if(!is.na(p$correlation_observed) & !is.na(p$correlation_permuted) ){
				sim = rbind(sim,c(p$correlation_permuted,p$correlation_observed))
				count = count + 1}
			fcount = fcount + 1
			if(fcount >  100 & (fcount-count)/fcount > 0.9){
				stop('>90% of simulations yield either no victims or no survivors. There is likely too high a ratio of survivors or victims in the data.')}
		}
	}
	return(sim)
}



perspectev.simulate <-
function(data,simulations,cores,traitfun=mcpRange,vlist=NULL,binary=NA,intercept=NA,slope=NA,level=NA,noise=0,fit=FALSE,quiet=FALSE){
	if(is.na(intercept) || is.na(slope)){
		if(!quiet){print('Intercept and slope parameters not specified. Assigning model fit through logistic regression.')}
		fit=TRUE
	}
	if(is.na(level)){
		if(!quiet){print("Level parameter not specified (can be 'lower' or 'upper'). Performing upper level simulation.")}
		level = 'upper'
	}
	if(level != 'lower' && level != 'upper'){
		stop("Level parameter can only be 'lower' or 'upper'")
	}
	if(is.na(binary)){
		if(sum(data$Survivorship$Survivorship == 1) + sum(data$Survivorship$Survivorship == 0) == nrow(data$Survivorship)){
			binary=TRUE}else{binary=FALSE}
			if(!quiet){print(paste("Binary/continuous survivorship not specified, setting binary =",binary))}
	}
	g = perspectev.calc(data,traitfun,vlist)
	lower = g$lower
	upper = g$upper

	if(level == 'lower'){
		if(fit == TRUE){ #base parameters on fitting logistic regression to data
			#Note: in continuous case this involves fitting a binomial regression to 
			#non-binary data. 
			gl = suppressWarnings(glm(lower[,1] ~ lower[,2],family='binomial'))
			intercept = gl$coefficients[1]
			slope = gl$coefficients[2]
		}
		if(!quiet){print(paste("Lower Level Selection Intercept",intercept))}
		if(!quiet){print(paste("Lower Level Selection Slope",slope))}
	}else{
		if(fit == TRUE){
			gl = suppressWarnings(glm(upper[,1] ~ upper[,2],family='binomial'))
			intercept = gl$coefficients[1]
			slope = gl$coefficients[2]
		}
		if(!quiet){print(paste("Upper Level Selection Intercept",intercept))}
		if(!quiet){print(paste("Upper Level Selection Slope",slope))}
	}
	names(intercept) = NULL
	names(slope) = NULL
	if(cores > 1){
		div = floor(simulations/cores)
		plist = rep(div,length=cores)
		plist[length(plist)] = plist[length(plist)] + simulations %% cores
		cl = makeCluster(cores)
		registerDoParallel(cl)
		i=NULL
		sim = foreach(i = 1:cores,.packages=c('perspectev'),.combine="rbind") %dopar% perspectev.simulation.permute(data,plist[i],intercept,slope,level,traitfun,vlist,binary,noise)
		stopCluster(cl)
	}else{
		sim = perspectev.simulation.permute(data,simulations,intercept,slope,level,traitfun,vlist,binary,noise)
	}
	total = list()
	total$level = level
	total$intercept = intercept
	total$slope = slope
	total$noise = noise
	total$fitted_model = fit
	total$binary = binary
	total$pvalue = mean(sim[,1] >= sim[,2])
	total$correlation_permuted = sim[,1]
	total$correlation_observed = sim[,2]
	return(total)
}

t1Range <- 
function(traitdata,vlist){ #calculate latitude breadth in decimal
	t1 = traitdata[[1]]
	if(class(t1) != "numeric"){stop("t1Range requires numeric data in first trait position!")}
	if(length(unique(t1)) < 2){return(NA)} #need at least two unique values
	return(max(t1) - min(t1))
}

t2Range <- 
function(traitdata,vlist){ #calculate longitude breadth
	t2 = traitdata[[2]]
	if(class(t2) != "numeric"){stop("t2Range requires numeric data in second trait position!")}
	if(length(unique(t2)) < 2){return(NA)} #need at least two unique values
	return(max(t2) - min(t2))
}

t1Var <- 
function(traitdata,vlist){ #calculate variance of trait 1
	t1 = traitdata[[1]]
	if(class(t1) != "numeric"){stop("t1Var requires numeric data in first trait position!")}
	if(length(unique(t1)) < 2){return(NA)} #need at least two unique values
	return(var(t1))
}

t2Var <- 
function(traitdata,vlist){ #calculate variance of trait 2
	t2 = traitdata[[2]]
	if(class(t2) != "numeric"){stop("t2Var requires numeric data in second trait position!")}
	if(length(unique(t2)) < 2){return(NA)} #need at least two unique values
	return(var(t2))
}

t1t2Covar <- 
function(traitdata,vlist){ #calculate covariance between trait 1 and trait 2
	t1 = traitdata[[1]]
	t2 = traitdata[[2]]
	if(class(t1) != "numeric"){stop("t1t2Covar requires numeric data in first and second trait positions!")}
	if(class(t2) != "numeric"){stop("t1t2Covar requires numeric data in first and second trait positions!")}
	if(nrow(unique(cbind(t1,t2))) < 2){return(NA)} #need at least two unique values
	return(cov(t1,t2))
}

t1Unique <- 
function(traitdata,vlist){ #calculate number of unique trait 1 values
	t1 = traitdata[[1]]
	if(length(t1) < 1){return(NA)} #need at least 1 unique value
	return(unique(t1))
}

t2Unique <- 
function(traitdata,vlist){ #calculate number of unique trait 2 values
	t2 = traitdata[[2]]
	if(length(t2) < 1){return(NA)} #need at least 1 unique value
	return(unique(t2))
}

t1t2Unique <- 
function(traitdata,vlist){ #calculate unique combinations of trait 1 and trait 2
	t1 = traitdata[[1]] 
	t2 = traitdata[[2]] 
	if(nrow(unique(cbind(t1,t2))) < 2){return(NA)} #need at least two unique values
	l = unique(data.frame(cbind(t1,t2))) 
	return(nrow(l))
}

mcpRange <- 
function(traitdata,vlist){ #calculate log minimum convex polygon range
	t1 = traitdata[[1]] #latitiude
	t2 = traitdata[[2]] #longitude
	if(class(t1) != "numeric"){stop("mcpRange requires numeric data in first and second trait positions!")}
	if(class(t2) != "numeric"){stop("mcpRange requires numeric data in first and second trait positions!")}
	l = unique(data.frame(cbind(t1,t2))) 
	if(nrow(l) < 3){return(NA)} #mcp needs at least three unique values
	range = attributes(hr.mcp(cbind(l[,1],l[,2]),plot=FALSE,n.min=3))$area
	if(range == 0){return(NA)} #close points can be rounded down to zero
	return(log(range))
}

dnaTheta <-
function(traitdata,vlist){ #calculate mean pairwise difference for set of sequences
	t1 = traitdata[[1]]
	if(length(t1) < 2){return(NA)} #need at least two sequences
	if(is.na(vlist)){stop("dnaTheta requires vlist to be specified")}#must provide global distance matrix
	if(class(vlist) != "list"){stop("dnaTheta requires vlist be a list")} #vlist must be a list
	compMat = vlist[[1]][t1,t1]
	return(mean(compMat[upper.tri(compMat)]))
}


perspectev.plot <-
function(observed,simulated,names,title){
	if(missing(observed)){
		stop("Must provide results from perspectev.test as 'observed'!")
	}
	if(missing(simulated)){
		simulated=list()
		names = c()
	}
	if(class(simulated) != 'list'){
		stop("'simulated' must be of type 'list'")
	}
	if(missing(names)){
		stop("Must provide vector of names to simulations using 'names' parameter!")
	}
	if(missing(title)){
		title = 'Difference between observed and permutated genera'
	}
	if(length(names) < 4){
		cbbPalette <- c("#2c7bb6","#abd9e9","#d7191c","#fdae61")
	}else if(length(names) < 9 && length(names) < 6){
		cbbPalette = c("#4575b4","#91bfdb","#e0f3f8","#fee090","#fc8d59","#d73027")
	}else if(length(names) < 9){
		cbbPalette <- c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090","#fdae61","#f46d43","#d73027")
	}else{
		stop("Cannot graph more than eight simulations at a time")
	}
	obs = observed$correlation_observed - observed$correlation_permuted
	tm = cbind('Observed',obs)
	for(l in 1:length(names)){
		temp1 = simulated[[l]]$correlation_observed -  simulated[[l]]$correlation_permuted
		temp2 = cbind(names[l],temp1)
		tm = rbind(tm,temp2)
	}
	tm = data.frame(tm)
	names(tm) = c('Simulation','RT')
	Simulation=RT=NULL
	tm$Simulation = factor(tm$Simulation,levels=c(names,'Observed'))
	tm$RT = as.numeric(as.character(tm$RT))

	ggplot(tm,aes(x=RT,fill=Simulation)) + 
		geom_density(alpha=0.8) + 
		xlab("R - S") +  ylab("Density") + 
		ggtitle(title) +
		guides(fill=guide_legend(title=NULL)) +
		scale_fill_manual(values = cbbPalette) +
		geom_vline(xintercept=0)
}


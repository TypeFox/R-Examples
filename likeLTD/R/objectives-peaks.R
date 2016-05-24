# Creates a likelihood function from the input hypothesis
# Documentation is in man directory.
create.likelihood.vectors.peaks <- function(hypothesis, addAttr=FALSE, 
				likeMatrix=FALSE, diagnose=FALSE, ...) 
	{
	# add arguments to hypothesis
	hypothesis = add.args.to.hypothesis(hypothesis, ...)
	# check hypothesis has the right type
	sanity.check.peaks(hypothesis)
	# convert hypothesis to locus specific
	locusCentric = transform.to.locus.centric.peaks(hypothesis)
	# functions to perform on each locus
	functions <- mapply(create.likelihood.per.locus.peaks, locusCentric,
                      MoreArgs=list(addAttr=addAttr, likeMatrix=likeMatrix, diagnose=diagnose))
	#create.likelihood.per.locus.peaks(locusCentric[[1]],addAttr=addAttr, 
	#				likeMatrix=likeMatrix, diagnose=diagnose)
	# set penalties if missing
	if(is.null(hypothesis$degradationPenalty)) hypothesis$degradationPenalty = 50 
	if(is.null(hypothesis$gradientSshape)) hypothesis$gradientSshape=0.013/(0.01^2/0.013)
	if(is.null(hypothesis$gradientSscale)) hypothesis$gradientSscale=0.01^2/0.013
	if(is.null(hypothesis$gradientAdjustSD)) hypothesis$gradientAdjustSD=0.2
	if(is.null(hypothesis$meanDshape)) hypothesis$meanDshape=0.02/0.018
	if(is.null(hypothesis$meanDscale)) hypothesis$meanDscale=0.018
	if(is.null(hypothesis$meanOshape)) hypothesis$meanOshape=0.02/0.018
	if(is.null(hypothesis$meanOscale)) hypothesis$meanOscale=0.018
	if(is.null(hypothesis$scaleRate)) hypothesis$scaleRate=0.01



	# output function, to be run every iteration
	likelihood.vectors <- function(degradation=NULL, DNAcont=NULL, 
					scale=NULL, gradientS = NULL, 
					gradientAdjust=NULL, #interceptAdjust=NULL, 
#					locusAdjust=NULL,
					dropin=NULL, dropinDeg=NULL, #interceptS=NULL, 
					meanD=NULL, meanO=NULL, repAdjust=NULL, 
					degradationPenalty=hypothesis$degradationPenalty,
					gradientSshape=hypothesis$gradientSshape,
					gradientSscale=hypothesis$gradientSscale,
					gradientAdjustSD=hypothesis$gradientAdjustSD,
					meanDshape=hypothesis$meanDshape,
					meanDscale=hypothesis$meanDscale,
					meanOshape=hypothesis$meanOshape,
					meanOscale=hypothesis$meanOscale,
					scaleRate=hypothesis$scaleRate, 
					...) 
		{
    		# set arguments
		arguments = list(degradation=degradation, DNAcont=DNAcont,
				scale = scale, repAdjust=repAdjust,
				gradientS = gradientS,#interceptS=interceptS, 
				meanD = meanD,meanO=meanO,
				dropin=dropin,dropinDeg=dropinDeg,
				degradationPenalty=degradationPenalty, 
				gradientSshape=gradientSshape,
				gradientSscale=gradientSscale,
				gradientAdjustSD=gradientAdjustSD,
				meanDshape=meanDshape,
				meanDscale=meanDscale,
				meanOshape=meanOshape,
				meanOscale=meanOscale,
				scaleRate=scaleRate)
		# single locus objective function
		callme <- function(objective,grad#,
		#stut,
#					loc
					) 
			{
			args = append(arguments, list(gradientAdjust=grad#,
							#interceptAdjust=stut,
#							locusAdjust=loc
							))
			do.call(objective, args)
			}
		if(length(gradientAdjust) == 1) gradientAdjust = rep(gradientAdjust, length(functions))
		#if(length(interceptAdjust) == 1) interceptAdjust = rep(interceptAdjust, length(functions))
#		if(length(locusAdjust) == 1) locusAdjust = rep(locusAdjust, length(functions))
		# call every objective function (one per locus)
		objectives = mapply(callme, functions, gradientAdjust#,#interceptAdjust,
		#locusAdjust
		)
		arguments = append(arguments, list(...))
		if(diagnose==TRUE) return(objectives)
		# calculate penalties
		pens <- do.call(penalties.peaks, append(arguments,list(gradientAdjust=gradientAdjust,
				#interceptAdjust=interceptAdjust,
#				locusAdjust=locusAdjust,
				nloc=ncol(hypothesis$queriedProfile))))
    		list(objectives=objectives, penalties=pens)
  		}
	if(addAttr) 
		{
		attr(likelihood.vectors, "hypothesis") <- hypothesis
		attr(likelihood.vectors, "functions") <- functions
		}
	return(likelihood.vectors)
	}

alterHypothesis = function(alleleDb,gens)
    {
    nRow = nrow(alleleDb)
    missing = unique(gens[which(!gens%in%as.numeric(rownames(alleleDb)))])
    if(length(missing)>0)
        {
        alleleDb = rbind(alleleDb,matrix(rep(alleleDb[which(rownames(alleleDb)==-1),],times=length(missing)),nrow=length(missing),byrow=TRUE))
        rownames(alleleDb)[(nRow+1):(nRow+length(missing))] = missing
        }
    return(alleleDb)
    }

alterDBvals = function(alleleDb,hypothesis)
	{
	dbVals = as.numeric(rownames(alleleDb))
	outVals = c(dbVals,dbVals-1)
	if(hypothesis$doDoubleStutter) outVals = c(outVals,dbVals-2)
	if(hypothesis$doOverStutter) outVals = c(outVals,dbVals+1)
	return(sort(unique(outVals)))
	}


create.likelihood.per.locus.peaks <- function(hypothesis, addAttr=FALSE, 
				likeMatrix = FALSE, diagnose=FALSE) 
	{
	# Creates a likelihood function for a given hypothesis and locus
	#
	# A hypothesis is given by the number of unknown contributors, whether to model
	# dropin, so on and so forth.
	cons = likelihood.constructs.per.locus.peaks(hypothesis)
	# alter hypothesis so that rare alleles are always different
    hypothesis$alleleDb = alterHypothesis(hypothesis$alleleDb,cons$genotypes)
	# alter dbVals to include any extra rare alleles
	cons$dbVals = alterDBvals(hypothesis$alleleDb,hypothesis)
	doR = !is.null(hypothesis$doR) && hypothesis$doR == TRUE

	result.function <- function(scale,gradientS,gradientAdjust,#interceptAdjust,
#	locusAdjust,
				#interceptS,
				meanD=NULL,meanO=NULL,repAdjust=NULL,
				degradation=NULL, DNAcont=NULL, 
				dropin=NULL, dropinDeg=NULL, ...) 
		{
		# Likelihood function for a given hypothesis and locus
		#
		# This function is specific to the hypothesis for which it was created.
		# It hides everything except the nuisance parameters over which to
		# optimize.
		#
		# Parameters:
		#   degradation: relative degradation from each profiled individual in this
		#                hypothesis
	    	#   ...: Any other parameter, e.g. for penalty functions. genotypeArray=cons$genotypes,
	    	#        These parameters are ignored here.
	    	# Returns: A scalar value giving the likelihood for this locus and
	    	#          hypothesis
		# complete repAdjust vector
		repAdjust = c(1,repAdjust)
		# perform some checks
		if(length(DNAcont) != hypothesis$nUnknowns + ncol(cons$knownPresence))
		stop(sprintf("DNAcont should be %d long.",
	                   hypothesis$nUnknowns + ncol(cons$knownPresence)))
		if(length(degradation) != hypothesis$nUnknowns +
	                              ncol(cons$knownPresence))
		stop(sprintf("degradation should be %d long.",
	                   hypothesis$nUnknowns + ncol(cons$knownPresence)))
		if(any(DNAcont < 0)) stop("found negative DNA contribution.")
		if(any(degradation < 0)) stop("found negative degradation parameter.")
		if(hypothesis$doDropin && is.null(dropin)) 
		stop("Model requires missing argument 'dropin'")
		if(hypothesis$doDropin && is.null(dropinDeg)) 
		stop("Model requires missing argument 'dropinDeg'")
		else if(is.null(dropin)) dropin = 0
		if(hypothesis$doDropin && dropin < 0) 
		stop("Dropin rate should be greater than 0.")
		if(diagnose==TRUE)
			{
			# diagnose result (for computing Z values)
			repRes <- peaks.probabilities(hypothesis=hypothesis, cons=cons, 
							DNAcont=DNAcont,scale=scale,
							gradientS = gradientS, 
							gradientAdjust=gradientAdjust,
#							locusAdjust=locusAdjust,
							#interceptAdjust=interceptAdjust, 
							#interceptS=interceptS,
							meanD = meanD, meanO=meanO,
							degradation=degradation, 
							repAdjust=repAdjust,
							detectionThresh=hypothesis$detectionThresh,
							dropin=dropin,dropinDeg=dropinDeg,
							doR=doR,diagnose=diagnose)
			return(repRes)
			}
		# result
		res <- peaks.probabilities(hypothesis=hypothesis, cons=cons, DNAcont=DNAcont, 
					scale=scale,gradientS = gradientS,gradientAdjust=gradientAdjust, 						#locusAdjust=locusAdjust,					
					#interceptAdjust=interceptAdjust,interceptS=interceptS, 
					meanD = meanD,meanO=meanO,degradation=degradation, 
					repAdjust=repAdjust,dropin=dropin,dropinDeg=dropinDeg,
					detectionThresh=hypothesis$detectionThresh,doR=doR)
		# multiply by allele probabilities
		factorsRes = res*cons$factors
		# Figure out likelihood for good and return.
		if(likeMatrix==FALSE)
			{
			return(sum(factorsRes))
			} else {
			return(list(evProb=res,genProb=cons$factors))
			}
		}
	# add some attributes
	if(addAttr) 
		{
		attr(result.function, "hypothesis") <- hypothesis
		attr(result.function, "constructs") <- cons
		}
	# return the result function to be performed every iteration
	result.function
	}


get.database.values = function(alleleDb, doOverStutter, doDoubleStutter)
	{
	# convert database alleles to numeric
	db = as.numeric(rownames(alleleDb))
	# add stutter allele to db vals
	out = c(db, db-1)
	# add over stutter values
	if(doOverStutter) out = c(out, db+1)
	# add double stutter values
	if(doDoubleStutter) out = c(out, db-2)
	sort(unique(round(out,1)))
	}


likelihood.constructs.per.locus.peaks = function(hypothesis) 
	{
	# Creates the locus-specific data needed by each likehood function.
	#
	# Parameters:
	#   hypothesis: A hypothesis, for instance one returned by
	#               prosecution.hypothesis(...) or defence.hypothesis(...)
	alleles = rownames(hypothesis$alleleDb)
	if(is.null(alleles)) stop("Could not figure out alleles names.")
	alleles.vector = function(n) alleles %in% unlist(n)
	cspPresence     = apply(hypothesis$binaryProfile, 1, alleles.vector)
	knownPresence = apply(hypothesis$knownProfs, 1, alleles.vector)
	if(!is.matrix(cspPresence)) cspPresence = matrix(ncol=0, nrow=length(alleles))
	if(!is.matrix(knownPresence)) knownPresence = matrix(ncol=0, nrow=length(alleles))
	missingReps = apply(hypothesis$binaryProfile, 1, is.na)
	# get intial genotype array
	genotypes = explain.all.peaks(cspPresence,knownPresence,hypothesis$knownProfs,
					alleles,hypothesis$nUnknowns,hypothesis$peaksProfile,
					hypothesis$heightsProfile,hypothesis$doDoubleStutter,
					hypothesis$doOverStutter,hypothesis$doDropin)
	# get database values, including all stutter positions
	dbVals = get.database.values(hypothesis$alleleDb,hypothesis$doOverStutter,
					hypothesis$doDoubleStutter)
	# get index of which alleles are from known contributors
	#do not want population allele probabilities for known contributors
	if(nrow(hypothesis$knownProfs)>0) 
		{
		kIndex = (hypothesis$nUnknowns*2)+(1:(2*nrow(hypothesis$knownProfs)))
		} else {
		kIndex = c()
		}
	# get genotype probabilities, taking into account relatedness
	if(length(kIndex>0)) 
		{
		# Only take into account relatedness between Q and X under Hd 
		if(hypothesis$hypothesis=="defence")
			{
			# exclude known genotypes
		  	factors = genotype.factors(genotypes[-kIndex,,drop=FALSE], 
					hypothesis$alleleDb,hypothesis$nUnknowns, 
					hypothesis$doDropin,hypothesis$queriedProfile,
                             		hypothesis$relatedness) 
			} else {
			# exclude known genotypes
  			factors = genotype.factors(genotypes[-kIndex,,drop=FALSE], 
					hypothesis$alleleDb,hypothesis$nUnknowns, 
					hypothesis$doDropin,hypothesis$queriedProfile,
                             		c(0,0))
			}
		} else {
		# Only take into account relatedness between Q and X under Hd 
		if(hypothesis$hypothesis=="defence")
			{
			factors = genotype.factors(genotypes, hypothesis$alleleDb,
                             		hypothesis$nUnknowns, hypothesis$doDropin,
					hypothesis$queriedProfile,hypothesis$relatedness) 
			} else {
			factors = genotype.factors(genotypes, hypothesis$alleleDb,
					hypothesis$nUnknowns, hypothesis$doDropin,
					hypothesis$queriedProfile,c(0,0))
			}
		}
	# convert genotypes from indices to repeat number
	genotypes = matrix(as.numeric(rownames(hypothesis$alleleDb))[genotypes],ncol=ncol(genotypes))
	# make sure rare combined alleles are always different
	foo = function(x)
    {
    index = which(x<0&x>-100)
    if(length(index)>1)
        {
        x[index] = seq(from=-1,to=-99,by=-4)[1:length(index)]
        }
    return(x)
    }
    genotypes = apply(genotypes,MARGIN=2,foo)
	#genotypes = differentRares(genotypes)
	# output list
	list(cspPresence=cspPresence, knownPresence=knownPresence,
        	missingReps=missingReps,genotypes=genotypes, factors=factors,
       		freqMat=hypothesis$alleleDb[, 1],dbVals=dbVals)
	}


# function to be called at each iteration of maximisation
peaks.probabilities = function(hypothesis,cons,DNAcont,#locusAdjust,
			scale,gradientS,gradientAdjust,#interceptAdjust,
			#interceptS,
			meanD=NULL,meanO=NULL,degradation,
			repAdjust,detectionThresh,dropin=NULL,dropinDeg=NULL,doR=FALSE,diagnose=FALSE)
	{
	# combine mean and adjustment
	locusGradient = gradientS*gradientAdjust
	#locusIntercept = interceptS*interceptAdjust
	locusIntercept = 0
	locusDNAcont = DNAcont#*locusAdjust
	if(hypothesis$doDropin==TRUE)
		{
		# no dropin currently
		if(doR==TRUE|diagnose==TRUE)
			{
			# probabilities for each replicate
			probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) 
				peak.heights.per.locus(genotypeArray=cons$genotypes,
						alleles=hypothesis$peaksProfile[[x]],
						heights=hypothesis$heightsProfile[[x]],
						DNAcont=locusDNAcont,
						gradientS = locusGradient,
						meanD=meanD,meanO=meanO,
						interceptS=locusIntercept,
						scale=scale,degradation=degradation,
						fragLengths=hypothesis$alleleDb[,2],
						fragProbs=hypothesis$alleleDb[,1],
						LUSvals=hypothesis$alleleDb[,3],
						repAdjust=repAdjust[x],
						detectionThresh=detectionThresh,
						dropin=dropin,dropinDeg=dropinDeg,
						diagnose=diagnose))
			} else {
			if(!is.null(meanD)&!is.null(meanO))
				{
				# single, double and over stutter
		    		probs = .Call(.cpp.getProbabilitiesSDO_dropin,
					genotypeArray=cons$genotypes,
					DNAcont=rep(locusDNAcont,each=2), 
					gradientS=locusGradient,
					meanD=meanD,meanO=meanO,
					interceptS=locusIntercept,
					degradation=rep(1+degradation,each=2),
					fragLengths=hypothesis$alleleDb[,2],
					fragNames=as.numeric(rownames(hypothesis$alleleDb)),
					LUSvals = hypothesis$alleleDb[,3],
					alleles=hypothesis$peaksProfile,
					heights=hypothesis$heightsProfile,
					repAdjust=repAdjust,scale=scale,
					detectionThresh=detectionThresh,
					databaseVals = cons$dbVals,
					fragProbs=hypothesis$alleleDb[,1], 
					dropin=dropin,dropinDeg=1+dropinDeg)
				} else if(is.null(meanD)&is.null(meanO)) {
				# single stutter
		    		probs = .Call(.cpp.getProbabilitiesS_dropin,
					genotypeArray=cons$genotypes,
					DNAcont=rep(locusDNAcont,each=2), 
					gradientS=locusGradient,
					interceptS=locusIntercept,
					degradation=rep(1+degradation,each=2),
					fragLengths=hypothesis$alleleDb[,2],
					fragNames=as.numeric(rownames(hypothesis$alleleDb)),
					LUSvals = hypothesis$alleleDb[,3],
					alleles=hypothesis$peaksProfile,
					heights=hypothesis$heightsProfile,
					repAdjust=repAdjust,scale=scale,
					detectionThresh=detectionThresh,
					databaseVals = cons$dbVals,
					fragProbs=hypothesis$alleleDb[,1], 
					dropin=dropin,dropinDeg=1+dropinDeg)
				} else if(!is.null(meanD)&is.null(meanO)) {
				# single, double stutter
		    		probs = .Call(.cpp.getProbabilitiesSD_dropin,
					genotypeArray=cons$genotypes,
					DNAcont=rep(locusDNAcont,each=2), 
					gradientS=locusGradient,
					meanD=meanD,
					interceptS=locusIntercept,
					degradation=rep(1+degradation,each=2),
					fragLengths=hypothesis$alleleDb[,2],
					fragNames=as.numeric(rownames(hypothesis$alleleDb)),
					LUSvals = hypothesis$alleleDb[,3],
					alleles=hypothesis$peaksProfile,
					heights=hypothesis$heightsProfile,
					repAdjust=repAdjust,scale=scale,
					detectionThresh=detectionThresh,
					databaseVals = cons$dbVals,
					fragProbs=hypothesis$alleleDb[,1], 
					dropin=dropin,dropinDeg=1+dropinDeg)
				} else if(is.null(meanD)&!is.null(meanO)) {
				# single, over stutter
		    		probs = .Call(.cpp.getProbabilitiesSO_dropin,
					genotypeArray=cons$genotypes,
					DNAcont=rep(locusDNAcont,each=2), 
					gradientS=locusGradient,
					meanO=meanO,
					interceptS=locusIntercept,
					degradation=rep(1+degradation,each=2),
					fragLengths=hypothesis$alleleDb[,2],
					fragNames=as.numeric(rownames(hypothesis$alleleDb)),
					LUSvals = hypothesis$alleleDb[,3],
					alleles=hypothesis$peaksProfile,
					heights=hypothesis$heightsProfile,
					repAdjust=repAdjust,scale=scale,
					detectionThresh=detectionThresh,
					databaseVals = cons$dbVals,
					fragProbs=hypothesis$alleleDb[,1], 
					dropin=dropin,dropinDeg=1+dropinDeg)
				}
			}
		} else {
		if(doR==TRUE|diagnose==TRUE)
			{
			# probabilities for each replicate
			probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) 
				peak.heights.per.locus(genotypeArray=cons$genotypes,
						alleles=hypothesis$peaksProfile[[x]],
						heights=hypothesis$heightsProfile[[x]],
						DNAcont=locusDNAcont,
						gradientS = locusGradient,
						meanD=meanD,meanO=meanO,
						interceptS=locusIntercept,
						scale=scale,degradation=degradation,
						fragLengths=hypothesis$alleleDb[,2],
						LUSvals=hypothesis$alleleDb[,3],
						repAdjust=repAdjust[x],
						detectionThresh=detectionThresh,
						diagnose=diagnose))
			} else {
			if(!is.null(meanD)&!is.null(meanO))
				{
				# single, double and over stutter
		    		probs = .Call(.cpp.getProbabilitiesSDO,
					genotypeArray=cons$genotypes,
					DNAcont=rep(locusDNAcont,each=2), 
					gradientS=locusGradient,
					meanD=meanD,meanO=meanO,
					interceptS=locusIntercept,
					degradation=rep(1+degradation,each=2),
					fragLengths=hypothesis$alleleDb[,2],
					fragNames=as.numeric(rownames(hypothesis$alleleDb)),
					LUSvals = hypothesis$alleleDb[,3],
					alleles=hypothesis$peaksProfile,
					heights=hypothesis$heightsProfile,
					repAdjust=repAdjust,scale=scale,
					detectionThresh=detectionThresh,
					databaseVals = cons$dbVals)
			} else if(is.null(meanD)&is.null(meanO)) {
		   		# single stutter only
		    		probs = .Call(.cpp.getProbabilitiesS,
					genotypeArray=cons$genotypes,
					DNAcont=rep(locusDNAcont,each=2), 
					gradientS=locusGradient,
					interceptS=locusIntercept,
					degradation=rep(1+degradation,each=2),
					fragLengths=hypothesis$alleleDb[,2],
					fragNames=as.numeric(rownames(hypothesis$alleleDb)),
					LUSvals = hypothesis$alleleDb[,3],
					alleles=hypothesis$peaksProfile,
					heights=hypothesis$heightsProfile,
					repAdjust=repAdjust,scale=scale,
					detectionThresh=detectionThresh,
					databaseVals = cons$dbVals)
			
			} else if(!is.null(meanD)&is.null(meanO)) {
		    		# single and double stutter
		    		probs = .Call(.cpp.getProbabilitiesSD,
					genotypeArray=cons$genotypes,
					DNAcont=rep(locusDNAcont,each=2), 
					gradientS=locusGradient,
					meanD=meanD,
					interceptS=locusIntercept,
					degradation=rep(1+degradation,each=2),
					fragLengths=hypothesis$alleleDb[,2],
					fragNames=as.numeric(rownames(hypothesis$alleleDb)),
					LUSvals = hypothesis$alleleDb[,3],
					alleles=hypothesis$peaksProfile,
					heights=hypothesis$heightsProfile,
					repAdjust=repAdjust,scale=scale,
					detectionThresh=detectionThresh,
					databaseVals = cons$dbVals)

			} else if(is.null(meanD)&!is.null(meanO)) {
		    		# single and over stutter
		    		probs = .Call(.cpp.getProbabilitiesSO,
					genotypeArray=cons$genotypes,
					DNAcont=rep(locusDNAcont,each=2), 
					gradientS=locusGradient,
					meanO=meanO,
					interceptS=locusIntercept,
					degradation=rep(1+degradation,each=2),
					fragLengths=hypothesis$alleleDb[,2],
					fragNames=as.numeric(rownames(hypothesis$alleleDb)),
					LUSvals = hypothesis$alleleDb[,3],
					alleles=hypothesis$peaksProfile,
					heights=hypothesis$heightsProfile,
					repAdjust=repAdjust,scale=scale,
					detectionThresh=detectionThresh,
					databaseVals = cons$dbVals)
			}
		}
	# diagnose
	if(diagnose==TRUE) return(probs)
	# set probability of impossible genotypes to 0
	probs[is.na(probs)] = 0
        }
    return(probs)
    }





# get probability pf each genotype combination
peak.heights.per.locus = function(genotypeArray,alleles,heights,DNAcont,
			gradientS,meanD=NULL,meanO=NULL,interceptS,scale,
			degradation,fragLengths,fragProbs=NULL,LUSvals,repAdjust=NULL,
			detectionThresh,dropin=NULL,dropinDeg=NULL,diagnose=FALSE)
	{
	# get combined probability of all peaks for each genotype combination seperately
	Probs = apply(genotypeArray,MARGIN=2,FUN=function(x) probability.peaks(genotype=x,
	        alleles=alleles,heights=heights,
        	DNAcont=DNAcont,
        	gradientS=gradientS,meanD=meanD,
        	meanO=meanO,interceptS=interceptS,
        	scale=scale,degradation=degradation,
        	fragLengths=fragLengths,fragProbs=fragProbs,LUSvals=LUSvals,
        	repAdjust=repAdjust,
        	detectionThresh=detectionThresh,
		dropin=dropin,dropinDeg=dropinDeg,
	        diagnose=diagnose))
	# diagnose
	if(diagnose==TRUE) return(Probs)
	# give impossible values a probability of 0
	Probs = unlist(Probs)
	Probs[is.na(Probs)] = 0
	return(unlist(Probs))
	}


# get probability of every peak given current parameters
probability.peaks = function(genotype,alleles,heights,DNAcont,
		gradientS,meanD=NULL,meanO=NULL,interceptS,scale,
		degradation,fragLengths,fragProbs=NULL,LUSvals,repAdjust=NULL,
		detectionThresh,dropin=NULL,dropinDeg=NULL,diagnose=FALSE)
	{
	# convert some formats
	heights = unlist(heights)
	alleles = unlist(alleles)
	genotype = as.numeric(genotype)
	# get mean expected peak heights
	gammaMu = peak.height.dose(genotype=genotype,
	        alleles=alleles,heights=heights,
	        DNAcont=DNAcont,gradientS=gradientS,
	        meanD=meanD,meanO=meanO,interceptS=interceptS,
	        degradation=degradation,fragLengths=fragLengths,
		fragProbs = fragProbs,
	        LUSvals=LUSvals,repAdjust=repAdjust,
		dropin=dropin,dropinDeg=dropinDeg)
	names(heights) = unlist(alleles)
	# give peak heights to dropout alleles (height=0)
	peakHeights = unlist(heights)
	gammaMus = gammaMu
	dropoutIndex = which(!names(gammaMus)%in%names(peakHeights))
	if(length(dropoutIndex)!=0)
	    {
        allelesToAdd = unique(names(gammaMus)[dropoutIndex])
        toAdd = rep(0,times=length(allelesToAdd))
        names(toAdd) = allelesToAdd
        peakHeights = c(peakHeights,toAdd)
        peakHeights = peakHeights[order(names(peakHeights))]
	    }
	# sort
	gammaMus = gammaMus[order(names(gammaMus))]
	peakHeights = peakHeights[order(names(peakHeights))]
	# scale = variance/mu
	# scale = (stanDev^2)/mu
	# alpha = (mu^2)/variance
	# alpha = (mu^2)/(stanDev^2)
	# shape = mu/scale
	# create shape and scale
	shapesVec = gammaMus/scale
	scalesVec = rep(scale,times=length(shapesVec)) 
	# diagnose
	if(diagnose==TRUE) return(list(height=peakHeights,
	mu=gammaMus,sigma=sqrt(shapesVec*scalesVec^2)))
	# probability densities
	pdf = vector(length=length(peakHeights))
	dropoutIndex = which(peakHeights==0)
	if(length(dropoutIndex)>0)
		{
		if(length(dropoutIndex)!=length(peakHeights))
		    {
		    # discrete approximation to pmf
		    pdf[-dropoutIndex] = mapply(FUN=function(x,k,t) pgamma(q=x+0.5,shape=k,scale=t)-
		        pgamma(q=x-0.5,shape=k,scale=t), 
		            x=peakHeights[-dropoutIndex], 
		            k=shapesVec[-dropoutIndex], 
		            t=scalesVec[-dropoutIndex])       
		    }
		# dropout = integral from 0 to threshold
		pdf[dropoutIndex] = mapply(FUN=function(k,t) pgamma(q=detectionThresh,shape=k,scale=t), 
		        k=shapesVec[dropoutIndex], 
		        t=scalesVec[dropoutIndex])
		} else {
		# non-dropout = gamma point density
		pdf = mapply(FUN=function(x,k,t) pgamma(q=x+0.5,shape=k,scale=t)-
		    pgamma(q=x-0.5,shape=k,scale=t), 
		        x=peakHeights, 
		        k=shapesVec, 
		        t=scalesVec)
		}
	# set impossible values to 0 likelihood
	pdf[which(is.infinite(pdf))] = 0
	# output
	return(prod(pdf))
	}


# get mean expected peak height at each position
peak.height.dose = function(genotype,alleles,heights,DNAcont,
			gradientS,meanD=NULL,meanO=NULL,
            		interceptS,degradation,fragLengths,fragProbs=NULL,
			LUSvals,repAdjust=NULL,dropin=NULL,dropinDeg=NULL)
	{
	# positions of stutter alleles
	stutterPos = genotype-1
	allPos = c(genotype,stutterPos)
	if(!is.null(meanD)) 
		{
		# positions of double stutter alleles
		doubleStutterPos = genotype-2
		allPos = c(allPos,doubleStutterPos)
		}
	if(!is.null(meanO))
		{
		# positions of over stutter alleles
		overStutterPos = genotype+1
		allPos = c(allPos,overStutterPos)
		}
	allPos = unique(round(allPos,1))
	# index frag lengths and LUS values
	fragLengthIndex = sapply(genotype,FUN=function(x) which(names(fragLengths)==x))
	# base dose
	dose = repAdjust*rep(DNAcont,each=2)*rep(1+degradation,each=2)^-fragLengths[fragLengthIndex]
	# stutter rate
	stutterRate = interceptS+(gradientS*LUSvals[fragLengthIndex])
	# dose from stutters
	muS = dose * stutterRate	
	if(!is.null(meanD))
		{
		# dose from double stutters
		muSd = dose * meanD
		stutterRate = stutterRate + meanD
		}
	if(!is.null(meanO))
		{
		# dose from over stutters
		muSo = dose * meanO
		stutterRate = stutterRate + meanO
		}
	# dose from allelic
	muA = dose * (1 - stutterRate)
	# tidy
	names(muA) = genotype
	names(muS) = stutterPos
	allPos = round(allPos,1)
	stutterPos = round(stutterPos,1)
	genotype = round(genotype,1)
	if(!is.null(meanD))
		{
		names(muSd) = doubleStutterPos
		doubleStutterPos = round(doubleStutterPos,1)
		}
	if(!is.null(meanO))
		{
		names(muSo) = overStutterPos
		overStutterPos = round(overStutterPos,1)
		}
	# combine doses for each position
	if(!is.null(meanD)&!is.null(meanO))
		{
		# S+D+O
		muX = sapply(allPos,FUN=function(x) sum(muA[which(genotype==x)])+
		                                    sum(muS[which(stutterPos==x)])+
		                                    sum(muSd[which(doubleStutterPos==x)])+
		                                    sum(muSo[which(overStutterPos==x)]))
		names(muX) = allPos
		} else if(is.null(meanD)&!is.null(meanO)) {
		# S+O
		muX = sapply(allPos,FUN=function(x) sum(muA[which(genotype==x)])+
		                                    sum(muS[which(stutterPos==x)])+
		                                    sum(muSo[which(overStutterPos==x)]))
		names(muX) = allPos
		} else if(!is.null(meanD)&is.null(meanO)) {
		# S+D
		muX = sapply(allPos,FUN=function(x) sum(muA[which(genotype==x)])+
		                                    sum(muS[which(stutterPos==x)])+
		                                    sum(muSd[which(doubleStutterPos==x)]))
		names(muX) = allPos
		} else if(is.null(meanO)&is.null(meanD)) {
		# S
		muX = sapply(allPos,FUN=function(x) sum(muA[which(genotype==x)])+
		                                    sum(muS[which(stutterPos==x)]))
		names(muX) = allPos
		}
	# add dropin expected dose
	if(!is.null(dropin))
		{
		toAdd = NULL
		for(i in 1:length(fragProbs))
			{
			if(as.numeric(names(fragProbs)[i])>c(-2)|as.numeric(names(fragProbs)[i])<c(-99))
			    {
			    index = which(round(as.numeric(names(muX)),1)==round(as.numeric(names(fragProbs)),1)[i])
			    if(length(index)==0)
				    {
				    toAdd = c(toAdd,fragProbs[i]*dropin*(1+dropinDeg)^-fragLengths[i])
				    names(toAdd)[length(toAdd)] = names(fragProbs)[i]				    
				    } else {
				    muX[index] = muX[index]+fragProbs[i]*dropin*(1+dropinDeg)^-fragLengths[i] 
			    	    }
				}
			}
		if(!is.null(toAdd))
			{
			muX = c(muX,toAdd)
			muX = muX[order(names(muX))]
			}
		}
	return(muX)
	} 





# Penalties to apply to the likelihood.
# Documentation is in man directory.
penalties.peaks <- function(nloc, 
			degradation=NULL, degradationPenalty=50, 
			gradientS=NULL,gradientSshape=0.013/(0.01^2/0.013), gradientSscale=0.01^2/0.013,
                       gradientAdjust=NULL,gradientAdjustSD=0.2,
		       meanD=NULL,meanDshape=0.02/0.018,meanDscale=0.018,
		       meanO=NULL,meanOshape=0.02/0.018,meanOscale=0.018,
                       scale=NULL, scaleRate=0.01,
			dropin=NULL, dropinDeg=NULL, ...) {
    # set initial penalty
    result = 1
    # Normalizes by number of loci so product of penalties same as in old code.
    # Some penalties are per locus and other for all locus. Hence this bit of run
    # around.
    normalization = 1.0 / max(nloc, 1)
    # penalty on degradation
    result = result * exp(-sum(degradation) * degradationPenalty 
                           * normalization )
    # penalty on gradientS
    # mean = 0.013, var=0.1 (actual var=0.01^2) 
    result = result * (dgamma(gradientS,shape=gradientSshape,scale=gradientSscale) * normalization)
    # penalty on interceptS
    # mean = 0.001, var=0.1 (actual mean=-0.076, actual var=0.13^2) 
    #result = result * (dgamma(interceptS,shape=0.001/(0.0001/0.001),scale=0.0001/0.001) * normalization)
    # penalty on gradientAdjust
    result = result * dnorm(log10(gradientAdjust),mean=0, sd=gradientAdjustSD)
    # penalty on interceptAdjust
    #result = result * dnorm(log10(interceptAdjust),mean=0, sd=0.5)
    # penalty on meanD
    if(!missing(meanD) & !is.null(meanD))
        {
        # mean = 0.01, sd = 5e-5
        result = result * (dgamma(meanD,shape=meanDshape,scale=meanDscale) * normalization)
        }
    # penalty on meanO
    if(!missing(meanO) & !is.null(meanO))
	    {
	    # mean = 0.01, sd = 5e-5
	    result = result * (dgamma(meanO,shape=meanOshape,scale=meanOscale) * normalization)
	    }
    # penalty on scale
    result = result * (dexp(scale,rate=scaleRate) * normalization)
    # penalty on dropin
    if(!missing(dropin) & !is.null(dropin))
	    {
	    #result = result * (densigamma(dropin,alpha=1e-30,beta=5)  * normalization)
	    }
    # penalty on dropinDeg
    if(!missing(dropinDeg) & !is.null(dropinDeg))
	    {
    		result = result * exp(-dropinDeg * degradationPenalty 
                           * normalization )
	    }
    # penalty on locusAdjust
    #result = result * dnorm(log10(locusAdjust),mean=0, sd=0.2)
  return(result)
}



    

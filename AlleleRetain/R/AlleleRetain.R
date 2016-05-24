########################################################################
# Predict retention of a rare allele in a small, bottlenecked population
# Predict accumulation of inbreeding in the same population
# Author: Emily Weiser 
# Portions of code drawn from R package 'mohuasim' written by Murray Efford
# 15 May 2012: fixed bug for polygynous systems with juvenile starters
# 1 June 2012: put allele frequency back into output; restructured to use less RAM
# 13 June 2012: added "printplots"
# 16 July 2012: changed abbreviated arguments to full (e.g. "TRUE" instead of "T", "ncol" instead of "nc")
# 20 August 2012:  added "all" option to "startAge" argument
# 2 January 2013:  Fixed pedigree.summary to work with new version of package pedigree
# 30 June 2013: Updated user manual to fix error in pdfinfo.  No change to code or function.

########################################################################

##### Functions needed to run aRetain: ####

## Test for whether a value is a whole number
is.wholenumber<-function(testme, tol=.Machine$double.eps^0.5)  abs(testme - round(testme)) < tol

pad <- function (x, l, fill=0) c(x,rep(fill, l-length(x)))  ## pad to length l

Blcl <- function (p, n, alpha=0.05)
## Binomial confidence interval
## Wilson score interval with continuity correction after Newcombe (1998)
## Newcombe, Robert G. "Two-Sided Confidence Intervals for the Single Proportion:
##    Comparison of Seven Methods," Statistics in Medicine, 17, 857-872 (1998).
{
  z <- qnorm(1-alpha/2)
  L <- (2*n*p + z^2 - 1 - z * (z^2 - 2 -1/n + 4 * p *(n*(1-p) + 1))^0.5) / (2 * (n + z^2))
  ifelse (p == 0, 0, L)
}

Bucl <- function (p, n, alpha=0.05)
## Binomial confidence interval
## Wilson score interval with continuity correction after Newcombe (1998)
## Newcombe, Robert G. "Two-Sided Confidence Intervals for the Single Proportion:
##    Comparison of Seven Methods," Statistics in Medicine, 17, 857-872 (1998).
{
  z <- qnorm(1-alpha/2)
  U <- (2*n*p + z^2 + 1 + z * (z^2 + 2 -1/n + 4 * p *(n*(1-p) - 1))^0.5) / (2 * (n + z^2))
  ifelse (p == 1, 1, U)
}

pairoff <- function (ID, sex, maxpairs)
## ID is vector of ID numbers of single birds
## sex is vector of sex codes, one per ID
# e.g. pairoff (ID = 1:20, sex = sample(1:2, 20, replace=TRUE), 6)
{
    if ((length(ID) < 2) | (maxpairs < 0.1)) {
        mated <- matrix(nrow=0, ncol=2)
        unmated <- ID
    }
    else {
        temp <- split(ID, sex)
        n <- max(sapply(temp, length))
        if (length(temp)==1)
            temp <- cbind(temp[[1]], rep(0,n))  # only one sex
        else
            temp <- sapply(temp, pad, n)
        unmated <- getsingles(temp)
        mated <- dropsingles(temp)
        npair <- nrow(mated)
        if (npair > maxpairs) {
            unmated <- c(unmated, mated[(maxpairs+1):npair,])
            mated <- mated[1:maxpairs,]
        }
    }
    list (mated = mated, unmated = unmated)
}

polypair <- function(males, females, maleLRS)
# All females will gain a mate.
# maleLRS is a list of the average lifetime reproductive success of each corresponding male
{
	if (length(females) < 1 | length(males) < 1) {
		mated <- matrix(nrow=0, ncol=2)
		unmated <- c(males, females)
	}
	else {
		maleprob <- maleLRS / sum(maleLRS)	# percent chance of any male being chosen, relative to other males in the current list
		matedmales <- rep(NA, length(females))
		if(length(males) == 1){
			if(maleLRS > 0) {	# all females will mate with the only male
				mated <- matrix(cbind(females, rep(males, length(females))), ncol=2)
				unmated <- matrix(nrow=0, ncol=2)
			}
			else {		# no males with prob mating > 0; no females mate.
				mated <- matrix(nrow=0, ncol=2)
				unmated <- females		
			}	
		}	
		else {
			for(i in 1:length(females)){ matedmales[i] <- sample(males, size=1, prob=maleprob)}
			mated <- matrix(cbind(females, matedmales), ncol=2)	# list of mated pairs (males may be listed more than once)
			unmated <- males[!(males %in% mated)]
		}	
	}
	list (mated = mated, unmated = unmated)
}


getsingles <- function (prs) {
  if (length(prs)>2) {
    unmated <- prs[apply(prs,1,prod) < 0.1,, drop=FALSE]
    as.numeric(unmated[unmated > 0.1])
  }
  else if ((length(prs)==2) & (prod(prs) < 0.1))
    as.numeric(prs[prs>0])
  else
    numeric(0)
}

dropsingles <- function (prs) {
  if (length(prs)>2) {
    stillpaired <- apply(prs,1,prod) > 0.1  ## otherwise one or other has died
    prs[stillpaired,, drop=FALSE]
  }
  else if ((length(prs)==2) & (prod(prs) > 0.1))
    matrix(prs, ncol=2)
  else
    matrix(nrow=0, ncol=2)
}

breed <- function (prs, genes, jppr, MAXypF, ypFsex, youngSR) {
## Generate juveniles from the birds in prs at the rate of 'jppr' per pair
## parental genes are sampled as appropriate
## sex ratio is proportion (male)
## Input -
## prs	     2-column matrix with ID numbers of pair on each row
## genes     2-column matrix with lookup list of parental alleles; rows indexed by ID
## jppr	     average juveniles per pair... 
## youngSR   proportion male among juveniles (at end of natal breeding season)
	offspring <- matrix(nrow=0, ncol=5)
	dimnames(prs) <- NULL
	prlist <- c(prs)	# listed in same order as genes, jppr	
	for(i in 1:nrow(prs)){
		maleID <- prs[i,2]
		femaleID <- prs[i,1]
		if(ypFsex=="male") jp <- jppr[which(prlist==maleID)]	## avg num offspring/yr produced by this male
		if(ypFsex=="female") jp <- jppr[which(prlist==femaleID)]
		if(ypFsex=="both") jp <- mean(c(jppr[which(prlist==maleID)], jppr[which(prlist==femaleID)]))
		nrecruit <- rpois (1,jp)   ## random Poisson variate, mean = jp --> offspring produced this year		
		if(nrecruit > MAXypF) nrecruit <- MAXypF
		if (nrecruit == 0) temp <- matrix(nrow=0, ncol=5)
		else {
			temp <- cbind (sample( c(1,2), nrecruit,
				prob = c(1-youngSR, youngSR), replace = TRUE),    ## sex, 1=female
				sample(genes[which(prlist==femaleID),], nrecruit, replace = TRUE),  ## maternal allele
				sample(genes[which(prlist==maleID),], nrecruit, replace = TRUE), ## paternal allele
				femaleID, maleID)  							## dam and sire IDs
			dimnames(temp) <- NULL
		}
		offspring <- rbind(offspring, temp)
	}
  matrix(offspring, ncol=5)
}


addnew <- function (N, inisurv, startSR, exactSR, sourceN, q0) {

## Add new individuals (starters, migrants) with sex and genotypes for each
## Input -
## N		number to add
## All other input values are as specified for "aRetain."
	newadded <- matrix(nrow=0, ncol=3)
	N <- rbinom (1, size = N, prob = inisurv)
	if(N > 0) {
		if (exactSR) {
			nfemale <- round(N*(1-startSR))
			nmale   <- N - nfemale
			sex <- rep (c(1,2), c(nfemale, nmale))
		}
		else sex   <- sample (c(1,2), N, prob = c(1-startSR, startSR), replace = TRUE)  ## 1 = female

		if (!is.finite(sourceN)) {
			genes <- sample (c(0,1), N*2, prob = c(1-q0, q0),
				replace = TRUE)  ## 1 = rare allele
			genes <- matrix(genes, ncol=2)
		}
		else {
			n0 <- round(2*sourceN*(1-q0), digits=0)
			n1 <- 2*sourceN - n0
			sourcegenes <- c(rep(0,n0), rep(1,n1))
			sourcegenes <- sourcegenes[sample (2*sourceN)] ## randomly permute genes HWE
			sampled <- sample (1:sourceN, N, replace = FALSE)
			genes <- matrix(sourcegenes,ncol=2)[sampled,]
			genes <- matrix(genes, ncol=2)
		}
		newadded <- cbind(sex, genes[,1], genes[,2])
		dimnames(newadded) <- NULL    ## tidier to leave unlabelled
	}
	newadded
}	

newinfo <- function (new, type, startAge, youngperF, SDypF, MAXypF, meanMLRS, sdMLRS, firstID, y, adsurvival, nonbrsurv, mature, SenesAge, MaxAge, matingSys) {
## Get information for each new individual added to the population
## Input -
## type		type of individual to add (1 = initial founder, 2 = additional founder, 4 = migrant)
## firstID	number at which to start the list of unique IDs for the new individuals
## y		current year
## All other input values are as specified for "aRetain."
	n <- nrow(new)
	info <- matrix(nrow=0, ncol=12)
	if(nrow(new) > 0) {	
		typeL <- rep(type, n)
		if (startAge == "juvenile")	year <- rep(y, n)	# These will be one year old at year 1 (they are added in year 0)
		if (startAge == "young adult") year<-rep((y-mature),times=n)	# birth year for young adults
		if(startAge == "adult" | startAge == "all") {	
			# randomly select age of each individual based on the probability of selecting an adult of each age
			ages1 <- mature:MaxAge								# adult ages
			if(startAge == "all") ages1 <- 1:MaxAge				# all ages
			adsurvivalList <- rep(adsurvival, length(ages1))	# list of survival rate for each adult age
			for(i in 1:length(ages1))	{				# adjust survival rate for age
				if(ages1[i] > SenesAge) adsurvivalList[i] <- adsurvival - (adsurvival / (MaxAge - SenesAge)) * (ages1[i] - SenesAge)
				if(ages1[i] <= mature) adsurvivalList[i] <- nonbrsurv
			}	
			indivList <- rep(NA, length(ages1))	# list to fill of proportion of individuals alive at each age
			indivList[1] <- 1
			for(i in 2:length(indivList)) indivList[i] <- indivList[i-1]*adsurvivalList[i-1]	# proportion of individuals alive at each age.  This is used as the list of probabilities of an individual of each age being selected to add to the population.
			ages <- sample(x=ages1, size=n, prob=indivList, replace=TRUE)	# current ages
			year <- y - ages	# year of birth
		}
		ypF <- rnorm(n = n, mean = youngperF, sd = SDypF)
		indivID <- firstID:(firstID + n - 1)	
		damsire <- rep(NA, length(indivID))
		yrs <- rep(0, length(indivID))
		if(length(ypF) > 0){ for(i in 1:length(ypF)) { 
			if(ypF[i] < 0 ) ypF[i] <- 0				## remove any negatives
			if(ypF[i] > MAXypF ) ypF[i] <- MAXypF	## cannot average more than the maximum annual repro output
			}	}
		if(matingSys == "monogamy") mates <- rep(1, n)	
		else {
			if(sdMLRS == 0) mates <- rep(meanMLRS, n)
			else mates <- round(rgamma(n=n, shape = meanMLRS^2 / sdMLRS^2, scale = sdMLRS^2 / meanMLRS), digits=2)	# round to ensure that some will be straight 0s, not just very tiny numbers
			if(length(mates) > 0){ for(i in 1:length(mates)) { 
				if(mates[i] < 0) mates[i] <- 0			## remove any negatives
			}}	
		}	
		info <- matrix(cbind(indivID, damsire, damsire, new[,1], new[,2], new[,3], typeL, year, ypF, mates, yrs, yrs), ncol=12)
		dimnames(info) <- NULL    ## tidier to leave unlabelled
	}
	info
}
##############################################


###########################################################################################

aRetain <- function (q0 = 0.05, sourceN = Inf, startN = 20, startAge = "juvenile", startSR = 0.5, exactSR= FALSE, inisurv = 0.90, addN = 0, addyrs = c(0),  migrN = 0, migrfreq = 1, mpriority = FALSE, removeL = FALSE, K = 100, Klag = 0, KAdults = FALSE, reprolag = 0, mature = 1, matingSys = "monogamy", matingLength = "seasonal", meanMLRS = 1, sdMLRS = 0, reproAgeM = c(1:200), AgeOnMLRS = "age/age", nMatings = 1, retainBreeders = "male", MaxAge = 25, SenesAge = 10, adsurvivalF = 0.80, adsurvivalM = 0.80, nonbrsurv = 0.80, nonbrsurvK = 0.80, juvsurv = 0.80,  juvsurvK = 0.80, youngperF = 1.5, SDypF = 0.25, ypF1 = 1, ypF1yr = 1, MAXypF = 2, MAXypFK = 2, ypFsex = "female", youngSR = 0.5, trackall = TRUE, GeneCount = "adult", nyears = 50, nrepl = 100, nreplprint = 10, printplots = FALSE)
{

## Simulate frequency of rare neutral allele in a two-allele system 
##
## Calls functions: polypair, pairoff, getsingles, dropsingles, breed, addnew, newinfo, is.wholenumber
##
## Inputs :
##
##  q0	Frequency of rare allele in the source population (range 0-1); defaults to 0.05.
##  sourceN	Size of source population; must be > startN; defaults to Inf (infinite).
##  startN	Number of starters (or size of bottleneck); not all will become genetic founders.  Minimum 2.
##  startAge	Age class (“juvenile”, “young adult”, or “adult”) of starters, supplemental.  If “juvenile”, all individuals added are assigned age 0; if “young adult”, all are assigned age at maturity; if “adult”, ages are selected randomly based on the proportion of individuals in the source population expected to be of each age (based on the survival rates and senescence specified below).
##  startSR	Sex ratio (proportion male) of starters, supplementals, and migrants; defaults to 0.5 (must be between 0 and 1).
##  exactSR	Whether startSR gives the exact sex ratio of individuals released (TRUE) sexes are assigned randomly based on the probability given by startSR (FALSE); defaults to FALSE.
##  inisurv	Initial survival rate, as a proportion (range 0-1), of starters, supplementals, and migrants immediately post-release.  Annual mortality applies after this value is used.  Defaults to 0.90.
##  addN	List of numbers of individuals to release in years soon after population establishment (“supplementals”); defaults to 0.
##  addyrs	List of years in which to release supplementals.  Each year corresponds to the number of individuals in the same position in the addN list.  Defaults to 0.
##  migrN	Number of migrants to add (must be a whole number); defaults to 0.
##  migrfreq	Interval (number of years) at which to add migrN migrants; must be between 1 and nyears, below; defaults to 1.
##  mpriority	TRUE or FALSE:  whether migrants are given priority over locally produced offspring to recruit into any available breeding vacancies; defaults to FALSE.
##  removeL	TRUE or FALSE: whether to remove the corresponding number of locally produced adults to make room for migrants in the population; only necessary if retainBreeders = “both”/”female”/”male”; occurs even when population is below K.  Defaults to FALSE.
##  K	Carrying capacity (population ceiling); defaults to 100.
##  Klag	Number of years for which population is held at or below initial size (breeding still occurs); indicates a prolonged bottleneck.  Defaults to 0.
##  KAdults    	TRUE (K = number of adults) or FALSE (K = total individuals; subadults, nonbreeders, helpers are also subjected to the limit of K).  Defaults to FALSE.
##  reprolag    	Number of years after establishment for which no reproduction occurs.  Defaults to 0.
##  mature	Average age (in years) at sexual maturity (first breeding); defaults to 1.
##  matingSys	Mating system: “monogamy”, “polygyny”, or “polygynandry”. In any case, pairs are formed and breed.  With polygyny, each male can be part of more than one pair.  With polygynandry, females mate and reproduce multiple times each year.  To model a polyandrous system, set to “polygyny” and then input female values for the “male” parameters (and male values for the “female” parameters) in the model.  Defaults to “monogamy”.
##  matingLength	“seasonal” or “lifelong”.  Determines whether individuals retain the same mate from year to year or divorce.  Note that if set to “lifelong” with polygyny/polygynandry, males will not obtain any new mates until ALL of their previous mates have died (probably not realistic in most cases).  Defaults to “seasonal”.
##  meanMLRS	Mean lifetime reproductive success (LRS), in terms of number of matings that produce young (NOT number of offspring) a male gets over his lifetime.  This is a population average for all males, including those that never reproduce, and may be a fraction.  Each male will be assigned an individual average from a gamma distribution with this mean and sdMLRS (the shape of the gamma function = (meanMLRS^2)/(sdMLRS^2); scale = (sdMLRS^2)/meanMLRS; see help for R function "rgamma" for more information).  The gamma distribution was chosen because of its flexibility in shape appropriate to polygynous mating systems (from strongly right-skewed to nearly symmetrical).  The SD:mean ratio is more important than the magnitude of the mean.  This individual mean indicates the male's “quality” and will be used to assess his chance of mating, relative to other males present, each year (does not translate directly into actual LRS experienced by that male).  Not used if matingSys = “monogamy”.  Defaults to 1.
##  sdMLRS	Among-male standard deviation in LRS.  Used with meanMLRS as described above.  If sdMLRS = 0, all males will have the same chance of breeding each year.  Not used if matingSys = “monogamy”.  Defaults to 0.
##  reproAgeM	List of ages at which males are able to mate successfully.  If males may mate at all ages, use 0:maximum possible lifespan.  This maximum cannot be set to Inf (unlike MaxAge) but you can use a really high number if unsure of lifespan.  Defaults to c(1:200).
##  AgeOnMLRS	Expression describing the proportion of LRS achieved by a male at a particular age (for ages contained within reproAgeM).  The user specifies the form of this expression; it must include “age” (the individual's current age) and no other undefined variables.  Example:  “-5.4 + 1.5*age - 0.08*age^2” describes a parabolic relationship between age and mating success (proportion of LRS achieved at each age).  If there is no effect of age, use the default value of “age/age” (equals 1 so all ages will be assigned the same average, given by meanMLRS).  If a given age is not included in reproAgeM, reproductive output at that age will be set to 0 regardless of the value calculated by this equation.
##  nMatings	Average number of matings per female each year.  Only used when matingSys = “polygynandry”.  Each female mates and breeds the corresponding number of times each year.  Repeat matings with the same male are not prohibited and may occur by chance.  Where there are multiple matings by females, the value entered for youngperF (below) is interpreted as the average number of offspring produced from each of these matings.  Defaults to 1; must be a whole number.
##  retainBreeders	Should established breeders retain their breeding status from year to year, and prevent young individuals from recruiting if the population is at K?  Specify which sex should be retained: “none”, “both”, “male”, or “female”.  Only used when matingSys = “monogamy.”  When “none”, all new recruits are added to the breeding population; individuals are randomly removed from that pool to truncate the population at K (so new recruits may randomly replace established breeders).  When adults will likely survive and prevent new individuals from recruiting, e.g. with territorial species, set this at one of the other values as appropriate for your species.  When pairing off widowed or divorced individuals, those of the retained sex(es) that bred previously will be guaranteed a new mate (if available); non-retained adults will compete with new recruits to mate with available adults.  If the population is at K, new recruits will only fill vacancies left by adults that died (they will not replace any surviving adults, including females when retainBreeders = “male” and vice versa; i.e. “both” functions the same as “male” and “female” in this part of the model).  Defaults to “male”.
##  MaxAge	Maximum allowable lifespan (in years); can be Inf.
##  SenesAge	Age (in years) after which annual survival will be reduced by senescence.  Through this age, adult survival values are set according to adsurvivalF and adsurvivalM (below).  After this age, annual survival decreases linearly until MaxAge (at which it is 0):  new survival =  survival - (survival / (MaxAge - SenesAge)) * (age - SenesAge).
##  adsurvivalF	Annual survival rate of adult females.  All survival rates are given as a proportion (e.g. 0.85).
##  adsurvivalM	Annual survival rate of adult males.
##  nonbrsurv	Annual survival rate of nonbreeders (subadults or adults that have never reproduced).  Chronologically mature individuals remain in the nonbreeding pool until they die or recruit and are subject to this survival rate.
##  nonbrsurvK	Annual survival rate of nonbreeders when population is at K (used instead of nonbrsurv).  If given, subadult survival probability in each year depends on density of the population at the beginning of that year, according to the Beverton-Holt function for density dependence in survival (as in Morris & Doak 2002, Quantitative Conservation Biology):  S(E(t)) = S(0)/(1 + beta * E(t)), where S(E(t)) is survival rate at population density E in year t, S(0) is survival when density is near 0, beta is the decline in survival as density increases, and E(t) is population density at time t.  "Density" is defined in our model as the proportion of K that has been filled, as there is no spatial information in the model.  The model solves for beta according to the user-specified values for nonbrsurv (S(0)) and nonbrsurvK (S at carrying capacity, where E = 1), then uses beta and S(0) to calculate density-dependent survival probability in each year.
##  juvsurv	First year survival (from the stage described by youngperF, below, to the beginning of the next breeding season) when population is below K.
##  juvsurvK	First year survival when population is at K (used instead of juvsurv).  Can be equal to juvsurv.  Otherwise, juvenile survival is density-dependent as for nonbrsurvK.
##  youngperF	Average number of offspring produced per mating each year (averaged over all pairs in population).  For a polyandrous female, this is the average number of offspring produced each time she mates with a male (each year): youngperF * nMatings = total average offspring per year.  youngperF can be calculated for any reproductive stage (eggs, chicks, independent juveniles) as long as juvsurv indicates the proportion of individuals that survive from this stage to the beginning of the following breeding season.  Given as offspring per pair, e.g. 1.5 or 0.75.
##  SDypF	Among-individual standard deviation of youngperF, e.g. 0.50 or 2.
##  ypF1    Where younger breeders have reduced reproductive rates, this can be used to define the reproductive success for first reproductive stage (length of that stage is determined by ypF1yr, below).  Given as a proportion of youngperF.  E.g. if youngperF = 2 and ypF1 = 0.5, mean reproductive success during the first stage will be 1 offspring per female.  Will be < 1 if inexperienced females experience lower reproductive success than older females; but can be > 1 to indicate higher reproductive success for younger females than for older females, e.g. with senescence (youngperF always applies after the age indicated by ypF1yr, below).
##  ypF1yr   	Age after which ypF1 changes to youngperF (e.g. 1 if ypF1 applies to one-year-olds only, or 5 if it applies for the first 5 years and then increases to youngperF from age 6 onward).
##  MAXypF	Maximum annual number of offspring per individual (e.g. based on biological constraints such as clutch size/renesting).  Note that if ypF1 > 1, this value may be exceeded by some individuals.
##  MAXypFK	Maximum annual number of offspring per individual when population is at K (if different from MAXypF).
##  ypFsex	Which member of a pair limits the reproductive output for the year, based on the biology of the species of interest.  Can be “male”, “female”, or “both” (the last will average the male's and female's values).
##  youngSR	Proportion of offspring that are male (range 0-1, defaults to 0.5).
##  trackall	Whether to track all individuals from the population through the whole simulation (TRUE or FALSE).  Must be TRUE if you wish to use indiv.summary or pedigree.summary after running the simulation (see Model Output, next page).  The simulations will be noticeably slower (and require more RAM) if this is set to TRUE, especially with larger carrying capacity, higher fecundity, more iterations, and longer simulation periods.  Defaults to FALSE.
##  GeneCount	Which alleles to count as retained:  those in the “adult” population only, or those in “all” individuals (including subadults/nonbreeders).
##  nyears	Number of years to run the simulation.  Consider setting a value that corresponds to 10 generations of your study-species. 
##  nrepl	Number of iterations (replicates) to run.  More iterations will take longer to run, but will produce narrower confidence limits for allele retention, inbreeding coefficients, and population size.  Defaults to 100, but 1000 may be more useful.
##  nreplprint	Interval (number of replicates) at which to print a message with the current system time.  Allows the user to gauge model progress and to estimate time to completion.  Defaults to 10.
##  printplots	logical: whether to plot the population growth (number of individuals, as defined by KAdults, present each year) and allele frequency (in the pool defined by GeneCount) as they change over time. (TRUE or FALSE)  One line will be plotted for each replicate immediately after it runs.  Can be used to immediately gauge the demographics of the population (e.g. if it will grow as expected); will slow down the simulation by ~ 10-20%.  Defaults to FALSE.
##  

## CHECK INPUT VALUES
if (q0<=0 | q0>=1) stop ('Initial gene frequency must be between 0 and 1')
if (startN > sourceN) stop ("'startN' cannot exceed 'sourceN'")
if (startN < 2) stop ("'startN' must be at least 2")
if (startAge != "juvenile" & startAge != "adult" & startAge != "young adult" & startAge != "all") stop ("'startAge' must be 'juvenile', 'young adult', 'adult', or 'all'")
if (startSR<=0 | startSR>=1 | youngSR<=0 | youngSR>=1) stop ("'startSR' and 'youngSR' must be between 0 and 1")
if (inisurv <=0 | inisurv > 1 | nonbrsurv <=0 | nonbrsurv > 1 | nonbrsurvK <=0 | nonbrsurvK > 1 | juvsurv <=0 | juvsurvK > 1 | adsurvivalF <=0 | adsurvivalF > 1) stop ("Survival rates must be between 0 and 1")
for(i in 1:length(addN)) { if (!is.wholenumber(addN[i])) stop ("'addN' must be a list of whole numbers")}
for(i in 1:length(addyrs)){if (!is.wholenumber(addyrs[i])) stop ("'addyrs' must be a list of whole numbers")}
if(length(addN) != length(addyrs)) stop ("'addN' and 'addyrs' must be vectors of the same length")
if (!is.wholenumber(migrN)) stop ("'migrN' must be a whole number")
if (!is.wholenumber(migrfreq)) stop ("'migrfreq' must be a whole number")
if (migrN > 0) if (migrfreq > nyears) stop ("'migrfreq' > 'nyears' so no migration will occur during this simulation; 'migrfreq' must be < 'nyears' or equal to 0.")
if (!is.wholenumber(K)) stop ("'K' must be a whole number")
if (!is.wholenumber(Klag)) stop ("'Klag' must be a whole number (# years)")
if (!is.wholenumber(mature)) stop ("'mature' must be a whole number (# years)")
if(MaxAge != Inf) if (!is.wholenumber(MaxAge)) stop ("'MaxAge' must be a whole number")
if(SenesAge != Inf) if (!is.wholenumber(SenesAge)) stop ("'SenesAge' must be a whole number")
if (!is.wholenumber(reprolag)) stop ("'reprolag' must be a whole number")
if (!is.wholenumber(nyears)) stop ("'nyears' must be a whole number")
if (matingLength != "seasonal") if (matingLength != "lifelong") stop("'matingLength' must be 'seasonal' or 'lifelong'")
if (retainBreeders != "none") if (retainBreeders != "both") if (retainBreeders != "male") if (retainBreeders != "female") stop ("'retainBreeders' must be 'none', 'both', 'male', or 'female'")
if (ypFsex != "female") if (ypFsex != "male") if (ypFsex != "both") stop ("'ypFsex' must be 'male', 'female', or 'both'")
if (trackall != TRUE) if(trackall != FALSE) stop("'trackall' must be TRUE or FALSE")
if (matingSys != "polygyny") if (matingSys != "monogamy") if (matingSys != "polygynandry") stop ("'matingSys' must be 'polygyny', 'polygynandry', or 'monogamy'")
if (meanMLRS < 0) stop ("'meanMLRS' must be > 0")
if (sdMLRS < 0) stop ("'sdMLRS' must be a positive value")
if (nMatings > 1) if (matingSys != "polygynandry" | matingLength == "lifelong") stop ("'matingSys' must be set to 'polygynandry' and 'matingLength' must be 'seasonal' if 'nMatings' > 1.  Pairs are re-shuffled for each of nMatings matings per season.")
if (!is.wholenumber(nMatings)) stop ("'nMatings' must be whole number")
if (GeneCount != "adult") if (GeneCount != "all") stop ("'GeneCount' must be 'all' or 'adult'")

##########################################
## SET UP REPLICATE FUNCTION
run.one.repl <- function (r) {

	## SET UP EMPTY OBJECTS NEEDED LATER
	out <- list()	# empty list for output for this replicate (one matrix will be added to the list for each year)
	census <- matrix(nrow=nyears, ncol=9)
	indivdata <- matrix(nrow=0, ncol=7)
	numpresent <- 0
	population <- matrix(nrow=0, ncol=13)
	nonbreeders <- matrix(nrow = 0, ncol=13)
    migrants <- matrix(nrow = 0, ncol=13)
	startersinfo <- matrix(nrow=0, ncol=13)
	juveniles <- matrix(nrow=0, ncol=13)
	oldpairs <- matrix(nrow=0, ncol=2)	
	firstpairs <- matrix(nrow=0, ncol=2)
	newadded <- matrix(nrow=0, ncol=13)
	indivID <- 0
	pairs <- matrix(nrow=0, ncol=2)
	singles <- numeric(0)

	## ADD STARTERS TO ESTABLISH POPULATION
    ## Individuals are added as if at the end of a breeding season (autumn).  They will not breed this year, and will be subject to annual mortality (on top of initial mortality) before breeding next year.
	starters <- addnew(startN, inisurv, startSR, exactSR, sourceN, q0)
	startersinfo <- matrix(cbind(newinfo(starters, type = 1, startAge, youngperF, SDypF, MAXypF, meanMLRS, sdMLRS, firstID = 1, y=1, adsurvival = mean(adsurvivalF, adsurvivalM), nonbrsurv=nonbrsurv, mature, SenesAge, MaxAge, matingSys), rep(0, nrow(starters))), ncol=13)
	indivID <- startersinfo[,1]

	
	## LOOP OVER YEARS
	for (y in 1:nyears) {	
		if(nrow(startersinfo) > 0){
			if(startAge == "juvenile") juveniles <- matrix(startersinfo, ncol=13) # put initial starters into juveniles matrix
			else nonbreeders <- matrix(startersinfo, ncol=13)
			indivID <- startersinfo[,1]
			startersinfo <- matrix(nrow=0, ncol=13)
		}		
		else juveniles <- matrix(nrow=0, ncol=13)	# reset each year (after year 1)
		newID2 <- numeric(0)
		oldpairs <- pairs		# last year's pairs to assess with retainBreeders this year
		numpresent <- length(unique(c(pairs, singles))) # num live adults to be used in density dependence later

		## BREEDING
		if (nrow(pairs) > 0) {
			if(y > reprolag) {			# do not breed if within reprolag period		
				prs <- c(pairs)		
				parentypp <- rep(NA, length(prs))				
				parentgenes <- matrix(cbind(rep(NA, length(prs)), rep(NA, length(prs))), ncol=2)
				for(i in 1:length(prs)){
					parentypp[i] <- population[which(population[,1]==prs[i]),9]	# indiv average repro output listed in same order as prs
					parentgenes[i,1] <- population[which(population[,1]==prs[i]),5]	# genes listed in same order as prs	
					parentgenes[i,2] <- population[which(population[,1]==prs[i]),6]
				}
				if(ypF1 != 1){		
					for(i in 1:length(parentypp)){
						if((y - population[which(population[,1]==prs[i]),8]) <= ypF1yr) parentypp[i] <- parentypp[i]*ypF1      # reduce young per pair for individuals in first repro stage
					}
				}			
				if(length(unique(c(pairs, singles))) > (K-1))	Mjpp <- MAXypFK
				else Mjpp <- MAXypF
				for(i in 1:nMatings){
					new <- breed (pairs, parentgenes, parentypp, Mjpp, ypFsex, youngSR)  ## offspring 	
					newinf <- matrix(cbind(newinfo(matrix(new[,1:3], ncol=3), type=3, startAge = "juvenile", youngperF, SDypF, MAXypF, meanMLRS, sdMLRS, firstID=(max(indivID)+1), y, adsurvival = mean(adsurvivalF, adsurvivalM), nonbrsurv, mature, SenesAge, MaxAge, matingSys), rep(y, nrow(new))), ncol=13)
					newinf[,2] <- new[,4]	# dam ID known for offspring
					newinf[,3] <- new[,5]	# sire ID known for offspring				
					juveniles <- matrix(rbind(juveniles, newinf), ncol=13)
					indivID <- c(indivID, newinf[,1])
					if (i < nMatings) {	# still more breeding left for this season; re-form all pairs
						singles <- sort(unique(c(singles, pairs)))	# split up all pairs if seasonal monogamy
						pairs <- matrix(nrow=0, ncol=2)	
						popalive <- matrix(subset(population, population[,1] %in% singles), ncol=13)
						males <- subset(popalive, popalive[,4]==2)
						males <- subset(males, males[,8] %in% c(y - reproAgeM))	# males that are of reproductive age
						females <- subset(popalive, popalive[,4]==1)
						maleID <- males[,1] # everything below is sorted by this order
						femaleID <- females[,1]
						malemates <- rep(NA, nrow(males))
						maleages <- y - males[,8]
						# Apply user-specified expression to determine each male's chance of breeding this year, based on his age:
						nMatesAtAge <- function(age) {}
						body(nMatesAtAge) <- parse(text = AgeOnMLRS)
						for(i in 1:length(maleages)){
							malemates[i] <- nMatesAtAge(maleages[i])*males[i,9]
							if(malemates[i] < 0) malemates[i] <- 0
						}
						if(sum(malemates) > 0){		# if all males have 0 chance of breeding, no pairs formed 
							newpairs <- polypair(maleID, femaleID, malemates)
							pairs <- matrix(newpairs$mated, ncol=2)
							singles <- singles[!(singles %in% pairs)]
						}
					}
				}	
			}
		}

		## ADD ANY ADDITIONAL FOUNDERS (will not breed until next year)
		if (sum(addN) > 0) { if (y %in% addyrs)	{
			addNy <- addN[which(addyrs==y)]
			addstarters <- addnew(addNy, inisurv, startSR, exactSR, sourceN, q0)
			addinfo <- matrix(nrow=0, ncol=13)
			addinfo <- matrix(cbind(newinfo(addstarters, type = 2, startAge, youngperF, SDypF, MAXypF, meanMLRS, sdMLRS, firstID = (max(indivID) + 1), y, adsurvival = mean(adsurvivalF, adsurvivalM), nonbrsurv, mature, SenesAge, MaxAge, matingSys), rep(y, nrow(addstarters))), ncol=13)
			indivID <- c(indivID, addinfo[,1])
			if (startAge == "juvenile") juveniles <- matrix(rbind(juveniles, addinfo), ncol=13)
			else nonbreeders <- matrix(rbind(nonbreeders, addinfo), ncol=13)	# adults go into nonbreeders matrix until they recruit to breed (or die)
		}}

		## IF y IS MIGRATION YEAR, ADD MIGRANTS.
		if (migrN > 0) if(migrfreq > 0) {
		if (is.wholenumber(y/migrfreq)){	# migration occurs this year
			migrants <- addnew(migrN, inisurv, startSR, exactSR, sourceN, q0)
			migrinfo <- matrix(cbind(newinfo(migrants, type = 4, startAge, youngperF, SDypF, MAXypF, meanMLRS, sdMLRS, firstID = (max(indivID) + 1), y, adsurvival = mean(adsurvivalF, adsurvivalM), nonbrsurv, mature, SenesAge, MaxAge, matingSys), rep(y, nrow(migrants))), ncol=13)
			indivID <- c(indivID, migrinfo[,1])
			if (startAge == "juvenile") juveniles <- matrix(rbind(juveniles, migrinfo), ncol=13)
			else nonbreeders <- matrix(rbind(nonbreeders, migrinfo), ncol=13)	# adults go into nonbreeders matrix until they recruit to breed (or die)
			if (removeL == TRUE) { # remove local adults to make room for migrants
				adultslive <- c(pairs, singles)
				adultslive <- adultslive[!(adultslive == 0)]
				popalive <- matrix(population[which(population[,1] %in% adultslive),], ncol=13)
				poplocals <- subset(popalive, popalive[,7]==3)
				if(nrow(poplocals) > migrN) Lremove <- sample(poplocals[,1], size=migrN)	# locals to remove			
				else Lremove <- poplocals[,1]	# remove ALL locals
				pairs[pairs %in% Lremove] <- 0	# set flag for later removal
				singles <- singles[!(singles %in% Lremove)]
			}
		}}
		
		## SURVIVAL
		# Adults.  Some may have been removed by removeL above, but that will not affect the probability of other individuals surviving.	
		adultslive <- unique(c(pairs, singles))
		adultslive <- adultslive[!(adultslive == 0)]
		popalive <- matrix(population[which(population[,1] %in% adultslive),], ncol=13)
		if (nrow(popalive) > 0) {
			fem <- subset(popalive, popalive[,4]==1)		# females
			mal <- subset(popalive, popalive[,4]==2)		# males
			survF <- rep(adsurvivalF, nrow(fem))	# female survival
			survM <- rep(adsurvivalM, nrow(mal))	# male survival		
			ageF <- y - fem[,8]						# female ages
			ageM <- y - mal[,8]						# male ages
			if(length(ageF) > 0) for(i in 1:length(ageF)) { 	# Individual survival probabilities based on age for females
				if (ageF[i] >= MaxAge) survF[i] <- 0
				else if (ageF[i] >= SenesAge) survF[i] <- adsurvivalF - (adsurvivalF / (MaxAge - SenesAge)) * (ageF[i] - SenesAge)
			}
			if(length(ageM) > 0) for(i in 1:length(ageM)) { 		# Individual survival probabilities based on age for males
				if (ageM[i] >= MaxAge) survM[i] <- 0
				else if (ageM[i] >= SenesAge) survM[i] <- adsurvivalM - (adsurvivalM / (MaxAge - SenesAge)) * (ageM[i] - SenesAge)
			}			
			survs <- c(survF, survM)
			IDs <- c(fem[,1], mal[,1])
			ages <- c(ageF, ageM)
			for (i in 1:length(survs)) if (survs[i] < 0) survs[i]<-0
			live <- rep(NA, length(survs))
			for (n in 1:length(live)) live[n] <- sample(c(0,1), size=1, prob=c(1-survs[i], survs[i]))
			for (i in 1:length(live)) if(ages[i] >= MaxAge) live[i] <- 0	# force mortality of old individuals			)
			liveID <- IDs * live
			liveID <- liveID[liveID != 0]
			pairs[!(pairs %in% liveID)] <- 0						# set flag for later removal		
			singles <- singles[singles %in% liveID]
		}
		if(sum(pairs)==0) pairs <- matrix(nrow=0, ncol=2)	
		
		# Subadults/nonbreeders.  Those that die are put into "population" to be recorded as dead individuals.
		# Note that nonbreeders were previously referred to as "subadults" so some of the nomenclature below will reflect that (SA/subad = subadult)
		if (nrow(nonbreeders) > 0) { 
			subadindex <- (1:nrow(nonbreeders))
			if(nonbrsurvK != nonbrsurv){
				beta <- nonbrsurv/nonbrsurvK - 1
				pden = numpresent/K	 
				sas <- nonbrsurv/(1 + beta*pden)
			}
			else sas <- nonbrsurv
			ss<-sample(c(0,1), size=length(subadindex), replace=TRUE, prob=c(1-sas, sas))
			alive <- subadindex * ss
			subadindex <- subadindex[subadindex %in% alive]
			nonbrdead <- matrix(nonbreeders[-subadindex,], ncol=13)
			nonbreeders <- matrix(nonbreeders[subadindex,], ncol=13)
			population <- matrix(rbind(population, nonbrdead), ncol=13)
		}

		# Juveniles.  Those that die are put into "population" to be recorded as dead individuals.
		if(juvsurvK < juvsurv){
			beta <- juvsurv/juvsurvK - 1
			pden = numpresent/K
			jsv <- juvsurv/(1 + beta*pden)
			}
		else jsv <- juvsurv
		mort <- 1-jsv
		if (nrow(juveniles) > 0) {
			juvindex <- 1:nrow(juveniles)
			js <- sample (c(0,1), size=length(juvindex), replace=TRUE, prob=c(mort,jsv))
			alive <- juvindex * js
			juvindex <- juvindex[juvindex %in% alive]
			juvdead <- matrix(juveniles[-juvindex,], ncol=13)
			juveniles <- matrix(juveniles[juvindex,], ncol=13)
			population <- matrix(rbind(population, juvdead), ncol=13)
		}
		nonbreeders <- matrix(rbind(nonbreeders, juveniles), ncol=13) 
		juveniles <- matrix(nrow=0, ncol=13)
	
		## RECRUIT INDIVIDUALS old enough to breed NEXT year.  
		 ## If retaining breeders, add only enough to replace adults that die (the rest remain in nonbreeders).  
		 ## If mpriority=TRUE, choose migrants type "4" over locals "3".
		nonbreeders3 <- matrix(nrow=0, ncol=13)
		if (nrow(nonbreeders) > 0) {
			x<-((y+1)-mature)	# birth year of individuals old enough to breed NEXT year
			if(length(subset(nonbreeders, nonbreeders[,8]<=x)) > 0) nonbreeders3<-matrix(subset(nonbreeders, nonbreeders[,8]<=x), ncol=13)	# nonbreeders old enough to breed			
			staySA <- matrix(subset(nonbreeders, nonbreeders[,8]>x), ncol=13)	 # nonbreeders not old enough
			if(retainBreeders == "none" | !KAdults) san <- nrow(nonbreeders3)	# allow all mature nonbreeders to recruit and potentially breed.  If specified, migrants will be prioritized later (during truncation).
			else {
				allad <- unique(c(pairs, singles))
				allad <- allad[allad != 0]
				san <- K - length(c(allad))	# Num subads needed to make up to K
				if (san < 0) san <- 0		
				if (san > 0) {
					if(migrN > 0) if (mpriority == TRUE) {  # Take maturing nonbreeders from previous migrants:		
						subadmigr <- matrix(subset(nonbreeders3, nonbreeders3[,7]==4), ncol=13)	# nonbreeders old enough to breed that came in as migrants
						subadlocal <- matrix(subset(nonbreeders3, nonbreeders3[,7]!=4), ncol=13)
						if (nrow(subadmigr) > 0) {
							if(san >= nrow(subadmigr)) {	# all mature migrants recruit
								SAbreed <- matrix(subadmigr, ncol=13)
								san <- san - nrow(SAbreed)	# Num subads still needed after migrants added
								nonbreeders3 <- matrix(subadlocal, ncol=13)	# All migrants used so only locals remain as nonbreeders
							}
							else {
								subadMindex <- 1:nrow(subadmigr)
								sel <- sample(subadMindex, san)	# Select from previous migrants if not all needed
								SAbreed <- matrix(subadmigr[sel,], ncol=13)
								SAMstay <- matrix(subadmigr[-sel,], ncol=13)
								nonbreeders3 <- matrix(rbind(subadlocal, SAMstay), ncol=13)	# Remaining mature nonbreeders
								san <- 0							# No more subads needed
							}
							if (nrow(SAbreed) > 0) {
								newID2 <- SAbreed[,1]
								population <- matrix(rbind(population, SAbreed), ncol=13)
								SAbreed <- matrix(nrow=0, ncol=13)
							}
						}
					}
				}
			}
			nonbreeders <- matrix(rbind(staySA, nonbreeders3), ncol=13)	# remaining nonbreeders that did not recruit (before final recruitment phase below; if no additional recruitment below, these will carry over to next year)

			# If no previous migrants available (or not enough), or if mpriority == FALSE, take from remaining mature subadult pool (nonbreeders3).  
			if (san > 0) {	
				if (nrow(nonbreeders3) > 0){
					subadindex3<-1:nrow(nonbreeders3)
					if (san < length(subadindex3)) subadindex3 <- sample(subadindex3, size=san)  
					SAmature <- matrix(nonbreeders3[subadindex3,],ncol=13)				
					if (nrow(SAmature) > 0) nonbreeders3 <- matrix(nonbreeders3[-subadindex3,],ncol=13)  # left over nonbreeders that do not mature
					else nonbreeders3 <- matrix(nrow=0, ncol=13)  # no nonbreeders left over
				}
				else SAmature <- matrix(nrow=0,ncol=13)
				if (nrow(SAmature) > 0){
					newID2 <- c(newID2, SAmature[,1])
					population <- matrix(rbind(population, SAmature), ncol=13)	# add to population
					dimnames(population) <- NULL
					nonbreeders <- matrix(rbind(staySA, nonbreeders3), ncol=13)
					SAmature <- matrix(nrow=0, ncol=13)
				}
			}
		}
		# Add new recruits to singles (already in population):	
		if(length(newID2) > 0) singles  <- c(singles, newID2)	
		
		## TRUNCATE AT CARRYING CAPACITY. This occurs only when KAdults = FALSE and/or when retainBreeders = "none".  When KAdults = true AND retainBreeders != "none", recruitment will have been limited so that K is not exceeded (and truncation does not occur).
		if (y < (Klag+0.1))	k <- startN
		else	k <- K
		adults2 <- unique(c(singles, pairs))
		adults <- adults2[adults2 != 0]		# live adults.  live subadults/nonbreeders are in "nonbreeders".
		## If !KAdults, nonbreeders are included in the grand total.  Subadult migrants will also get priority if indicated.
		## If KAdults, only adults are included.
		ALL <- adults
		if(!KAdults) if(nrow(nonbreeders) > 0 )	ALL <- c(adults, nonbreeders[,1])
		if(nrow(population) > 0) {  if(length(ALL) > k) {
			adultslive <- c(pairs, singles)
			adultslive <- adultslive[!(adultslive == 0)]
			popalive <- matrix(population[which(population[,1] %in% adultslive),], ncol=13)	# only live individuals can be selected	
			mig <- matrix(subset(popalive, popalive[,7]==4), ncol=13)
			miga <- mig[,1]	
			if(KAdults) {	## only adults are limited to K	
				# Truncation with KAdults will be necessary only when retainBreeders == none; so established breeders get no priority here.  Migrants get first priority; then all other adults compete to remain here.
				if(migrN > 0) if (mpriority == TRUE){		
					# migrants that are established breeders get first priority for staying in population
					if (length(miga) > k) alive <- sample(miga, k)		# enough migrants to make up K by themselves
					else {					
						nm <- ALL[!(ALL %in% miga)]		# locally produced adults ID
						n <- k - length(miga)			# number of locals needed to supplement migrants
						if(n > 0) {
							if(n < length(nm)) nma <- sample(nm, n) # keep number needed to make up K
							if(n > length(nm)) nma <- nm
							}
						else nma <- numeric(0)
						alive <- c(miga, nma)			# ID numbers of individuals to keep
					}
				}		
				if(migrN == 0 | mpriority == FALSE) alive <- sample (ALL, k)   ## retain random sample of k adults 	   
				singles  <- singles[singles %in% alive]
				pairs[!(pairs %in% alive)] <- 0
			}
			else {		## all individuals are limited to K (including nonbreeders, if any)
				# Truncation will be used regardless of retainBreeders; so give any retained adults first priority (as indicated).  Then select migrants from those remaining (adults and subadults/nonbreeders).  Then select randomly from any other available individuals to fill K.
				sbk <- numeric(0)
				if(length(ALL) > k) {		
					n <- k	# number needed (to be updated as prioritized individuals are added)
					keep <- numeric(0)
					left <- numeric(0)
					
					# Keep retained breeders in population (and calculate how many more are needed)
					if(retainBreeders == "both"){
						keep <- ALL[ALL %in% pairs]  # K is never reduced, so ALL birds that bred this year will breed next year (don't truncate any).
						left <- ALL[!(ALL %in% keep)]	# individuals not yet chosen (including nonbreeders, if any, and non-retained sex, if any)
						n <- n - length(keep)				
					}
					if(retainBreeders == "female"){
						keep <- ALL[ALL %in% pairs[,1]]
						left <- ALL[!(ALL %in% keep)]
						n <- n - length(keep)
					}
					if(retainBreeders == "male"){
						keep <- ALL[ALL %in% pairs[,2]]
						left <- ALL[!(ALL %in% keep)]
						n <- n - length(keep)
					}
					if(retainBreeders == "none") left <- ALL
				
					# Keep migrants in population, if mpriority = true (and calculate how many more are needed)
					if(migrN > 0) if (mpriority == TRUE) {
						miga <- miga[miga %in% left]	# only adult migrants that have NOT already been selected as retained adults
						if(nrow(nonbreeders) > 0) {
							sbm <- matrix(subset(nonbreeders, nonbreeders[,7]==4), ncol=13)	#  migrants
							sm <- sbm[,1]				# live subadult migrant ID	
							asm <- c(miga, sm)			# all live migrant ID	
						}
						else asm <- miga
						if (length(asm) > n){
							alive <- sample(asm, n)		# enough migrants to make up rest of K by themselves
							n <- 0							# no more needed to meet K
						}	
						else{
							alive <- asm				# all migrants are kept
							n <- n - length(alive)		# number still needed to meet K
						}	
						ak <- miga[miga %in% alive]	# adult migrants to keep
						keep <- c(keep, ak)
						if(nrow(nonbreeders) > 0) sbk <- sm[sm %in% alive]	# subadult migrants to keep
					}
					n <- k - length(c(keep, sbk))	# number still needed to make up K	
					
					# If only some migrants were kept, we have no room for the rest (or any locals).
					# If ALL retained breeders and ALL migrants are kept (or no immigrant priority) and we still need more, choose more to keep from locals:

					if(n > 0){	
						if(length(adults[!(adults %in% keep)]) > 0) anm <- adults[!(adults %in% keep)]		# adults not already selected to keep (local non-retained)
						else anm <- 0
						if(nrow(nonbreeders) > 0) {
							sbn <- matrix(subset(nonbreeders, nonbreeders[,7]!=4), ncol=13) # local nonbreeders
							sbnID <- sbn[,1]	# IDs of local nonbreeders = row number in "nonbreeders"
							nm <- c(anm, sbnID)	# all nonmigrant IDs
						}
						else nm <- anm		
						if(length(nm) > n) nma <- sample(nm, n)	# keep number of locals needed to make up K
						else nma <- nm							# keep all locals if room for them all in K
						anmk <- anm[anm %in% nma]	# adult locals to keep
						keep <- c(anmk, keep)			# adults to keep		
						if(nrow(nonbreeders) > 0) {
							sblk <- sbnID[sbnID %in% nma]	# subadult locals to keep
							sk <- c(sblk, sbk)		# nonbreeders to keep
							nonbreeders <- matrix(nonbreeders[which(nonbreeders[,1] %in% sk),], ncol=13)	# nonbreeders still in the population
							}	
						}
						singles  <- singles[singles %in% keep]
						pairs[!(pairs %in% keep)] <- 0
					} }
			} }
	
		## PAIRING
		# length(unique(c(pairs))) is the number of live paired individuals
		# oldpairs gives the list of pairs from end of last year (before mortality etc this year)
		# reassign all adults to pairs (or leave as singles) when seasonal mating; otherwise retain pairs from last year.
		# population truncation to K has already occurred, so the number of pairs that can form is not limited.
		adultslive <- c(pairs, singles)
		adultslive <- adultslive[!(adultslive == 0)]
		popalive <- matrix(population[which(population[,1] %in% adultslive),], ncol=13)
		if (nrow(popalive) > 0) {
			if (matingSys != "monogamy"){
				singles  <- unique(c(singles, getsingles(pairs)))  # add newly widowed to 'singles'
				pairs    <- dropsingles(pairs)  ## drop widowed birds from pair list
				singles <- sort(singles[!(singles %in% pairs)])	# male that lost one female will be put into 'singles' above, but may still have other females, so is not actually single - take him out here.
				if (matingLength == "seasonal") {
					singles <- sort(unique(c(singles, pairs)))	# split up all pairs if seasonal monogamy
					pairs <- matrix(nrow=0, ncol=2)
				}		# if mating is lifelong, previous pairs remain and only single individuals will be paired off
				if(length(singles) > 1) {	# if not, cannot make new pairs
					livesingles <- matrix(subset(population, population[,1] %in% singles), ncol=13)
					males <- subset(livesingles, livesingles[,4]==2)				
					reproyrs <- (y+1) - reproAgeM	# birth years of males that will be of reproductive age next year		
					males <- subset(males, males[,8] %in% reproyrs)	# males that will be of reproductive age	
					females <- subset(livesingles, livesingles[,4]==1)
					if(nrow(males) > 1) if(nrow(females) > 1){	# if not, cannot make new pairs
						maleID <- males[,1]
						femaleID <- females[,1]
						malemates <- rep(NA, nrow(males))
						maleages <- y+1 - males[,8]	# age NEXT year
						# Apply user-specified expression to determine each male's chance of breeding NEXT year, based on his age NEXT year:
						nMatesAtAge <- function(age) {}
						body(nMatesAtAge) <- parse(text = AgeOnMLRS)
						for(i in 1:length(maleages)){
							malemates[i] <- nMatesAtAge(maleages[i])*males[i,10]
							if(malemates[i] < 0) malemates[i] <- 0
						}				
						if(sum(malemates) > 0){		# if all males have 0 chance of breeding, no new pairs formed this year 
							newpairs <- polypair(maleID, femaleID, malemates)
							pairs <- rbind(pairs, newpairs$mated)
							singles <- singles[!(singles %in% pairs)]
						}
					}	
				}
			}	
			else {	# retainBreeders only comes into play with monogamy.
				if (retainBreeders=="none"){  
					singles  <- sort(c(singles, getsingles(pairs)))  # add newly widowed to 'singles'
					pairs    <- dropsingles(pairs)  ## drop widowed birds from pair list
					if (matingLength == "seasonal") {
						singles <- sort(c(singles, pairs))	# all pairs are split up
						pairs <- matrix(nrow=0, ncol=2)
					}	# otherwise pairs remain as they were, and only single/widowed birds are paired up now.	
					livesingles <- matrix(subset(population, population[,1] %in% singles), ncol=13)
					newpairs <- pairoff(livesingles[,1], livesingles[,4], Inf)	# new pairs from all adults				
					pairs <- rbind(pairs, newpairs$mated)
					singles <- newpairs$unmated
				}
				else {			
					fwid <- matrix(nrow=0, ncol=13)
					mwid <- matrix(nrow=0, ncol=13)
					fwidID <- numeric(0)
					mwidID <- numeric(0)
					mprID <- numeric(0)
					fprID <- numeric(0)
					msinID <- numeric(0)
					fsinID <- numeric(0)
					pd <- c(oldpairs)			
					pd <- sort(pd[pd %in% pairs])	# do not include previously paired adults that have died			
					prs <- matrix(subset(population, population[,1] %in% pd), ncol=13)	# Previously paired adults that are still alive	
					mpr <- subset(prs, prs[,4]==2)		# males previously paired
					fpr <- subset(prs, prs[,4]==1)		# females previously paired
					if(nrow(mpr) > 0) mprID <- c(mpr[,1])	# male paired ID
					if(nrow(fpr) > 0) fprID <- c(fpr[,1])	# female paired ID
					widows <- sort(getsingles(pairs))
					pairs <- dropsingles(pairs)		# drop widowed birds
					newpairs <- matrix(pairs[!(pairs %in% oldpairs)], ncol=2)  # not including indiv paired last year
					if (length(widows) > 0) {
						wid <- matrix(subset(population, population[,1] %in% widows), ncol=13)
						mwid <- matrix(subset(wid, wid[,4]==2), ncol=13)		# Select male widow
						fwid <- matrix(subset(wid,wid[,4]==1), ncol=13)		# Select female widows
						if(nrow(mwid)>0) mwidID <- c(mwid[,1])	# male widow ID
						if(nrow(fwid)>0) fwidID <- c(fwid[,1])	# female widow ID
						}
					newsingles <- sort(singles[!(singles %in% oldpairs)])	# does not included paired indiv from last year that are now single
					sing <- matrix(subset(population, population[,1] %in% newsingles), ncol=13)	
					fsin <- subset(sing,sing[,4]==1)				# single females
					msin <- subset(sing,sing[,4]==2)				# single males
					if(nrow(fsin) > 0) fsinID <- c(fsin[,1])	# single female ID
					if(nrow(msin) > 0) msinID <- c(msin[,1])	# single male ID
					if (retainBreeders == "both"){				
						if (matingLength == "lifelong") {
							toprm <- c(mwidID, fsinID)	# widowed males and single females to pair off
							toprf <- c(fwidID, msinID)	# widowed females and single males to pair off
							toprmpop <- matrix(subset(population, population[,1] %in% toprm), ncol=13)
							toprfpop <- matrix(subset(population, population[,1] %in% toprf), ncol=13)
							mp <- pairoff(toprmpop[,1], toprmpop[,4], Inf)	# pair off single females with established males
							fp <- pairoff(toprfpop[,1], toprfpop[,4], Inf)	# pair off single males with established females
							singles <- c(mp$unmated, fp$unmated)
							singpop <- matrix(subset(population, population[,1] %in% singles), ncol=13)
							mf <- pairoff(singpop[,1], singpop[,4], Inf)  # pair off any remaining singles
							singles <- mf$unmated
							pairs <- matrix(rbind(pairs, newpairs, mp$mated, fp$mated, mf$mated), ncol=2)
						}
						if (matingLength == "seasonal") {
							pdpop <- matrix(subset(population, population[,1] %in% pd), ncol=13)
							pp <- pairoff(pdpop[,1], pdpop[,4], Inf)	# pairoff previous pairs first
							pairs <- pp$mated
							ps <- pp$unmated	# previously paired individuals left single
							pspop <- matrix(subset(population, population[,1] %in% ps), ncol=13)
							psf <- subset(pspop, pspop[,4]==1)
							psfID <- psf[,1]		# previously paired females left single
							psm <- subset(pspop, pspop[,4]==2)
							psmID <- psm[,1]		# previously paired males left single
							if(nrow(pairs) < Inf) {
								if (length(psfID) > 0) {	# pair off previously paired females
									topr <- c(psfID, msinID)
									toprpop <- matrix(subset(population, population[,1] %in% topr), ncol=13)
									pf <- pairoff(toprpop[,1], toprpop[,4], Inf)
									pairs <- rbind(pairs, pf$mated)
								}
								if (length(psmID) > 0) {	# pair off previously paired males
									topr <- c(psmID, fsinID)
									toprpop <- matrix(subset(population, population[,1] %in% topr), ncol=13)
									pm <- pairoff(toprpop[,1], toprpop[,4], Inf)
									pairs <- rbind(pairs, pm$mated)
								}
							}
							fsinID <- fsinID[!(fsinID %in% pairs)]
							msinID <- msinID[!(msinID %in% pairs)]
							sinID <- c(msinID, fsinID)	# remaining single individuals to pair
							sinpop <- matrix(subset(population, population[,1] %in% sinID), ncol=13)
							pn <- pairoff(sinpop[,1], sinpop[,4], Inf)
							if(nrow(pairs) > 0) pairs <- rbind(pairs, pn$mated)
							else pairs <- pn$mated
							singles <- singles[!(singles %in% pairs)]
						}
					}
					if (retainBreeders == "male") {
						if (matingLength == "lifelong") if(nrow(mwid) > 0){
							singles <- sort(c(singles, fwidID))	# Intact pairs remain; widowed females join singles
							sing <- matrix(subset(population, population[,1] %in% singles), ncol=13)	
							fsin <- matrix(subset(sing,sing[,4]==1), ncol=13)
							if(nrow(fsin) > 0) fsinID <- c(fsin[,1])	# single female ID
							topr <- c(mwidID, fsinID)	# widowed males and single females to pair off
							toprpop <- matrix(subset(population, population[,1] %in% topr), ncol=13)
							mp <- pairoff(toprpop[,1], toprpop[,4], Inf)	# pair off single females with established males
							singles <- singles[!(singles %in% mp$mated)]
							if(length(mp$mated)>0){
								if(length(pairs)>0) pairs <- matrix(rbind(pairs, mp$mated),ncol=2)
								else pairs <- matrix(mp$mated,ncol=2)
							}
						}
						if (matingLength == "seasonal"){
							singles <- sort(c(singles, fprID))	# when no long-term mating, all females put into singles to re-pair
							sing <- matrix(subset(population, population[,1] %in% singles), ncol=13)	
							fsin <- matrix(subset(sing,sing[,4]==1), ncol=13)
							if(nrow(fsin) > 0) fsinID <- c(fsin[,1])	# single female ID
							topr <- c(mprID, fsinID)	# previously paired males and all females
							toprpop <- matrix(subset(population, population[,1] %in% topr), ncol=13)
							prs <- pairoff(toprpop[,1], toprpop[,4], Inf)	# pair off females with established males
							singles <- singles[!(singles %in% prs$mated)]
							if(length(prs$mated)>0) pairs <- prs$mated
							else pairs <- matrix(nrow=0, ncol=2)
						}
					}
					if (retainBreeders=="female") {
						if (matingLength == "lifelong") if(length(fwidID)>0){
							singles <- sort(c(singles, mwidID))	# Intact pairs remain; widowed males join singles
							sing <- matrix(subset(population, population[,1] %in% singles), ncol=13)
							msin <- matrix(subset(sing,sing[,4]==2), ncol=13)
							if(nrow(msin) > 0) msinID <- c(msin[,1])	# single male ID
							topr <- c(fwidID, msinID)	# widowed females and single males to pair off
							toprpop <- matrix(subset(population, population[,1] %in% topr), ncol=13)
							fp <- pairoff(toprpop[,1], toprpop[,4], Inf)	# pair off single males with established females
							singles <- singles[!(singles %in% fp$mated)]
							if(length(fp$mated)>0){
								if(length(pairs)>0) pairs <- matrix(rbind(pairs, fp$mated),ncol=2)
								else pairs <- matrix(fp$mated,ncol=2)
							}
						}
						if (matingLength == "seasonal") {
							singles <- sort(c(singles, mprID))	# when no long-term mating, all males put into singles to re-pair
							sing <- matrix(subset(population, population[,1] %in% singles), ncol=13)
							msin <- matrix(subset(sing,sing[,4]==2),ncol=13)
							if(nrow(msin) > 0) msinID <- c(msin[,1])	# single male ID
							topr <- c(fprID, msinID)	# previously paired females and all males
							toprpop <- matrix(subset(population, population[,1] %in% topr), ncol=13)
							prs <- pairoff(toprpop[,1], toprpop[,4], Inf)	# pair off any males with established females
							singles <- singles[!(singles %in% prs$mated)]
							if(length(prs$mated)>0) pairs <- matrix(prs$mated,ncol=2)
							else pairs <- matrix(nrow=0, ncol=2)
						}
					}
					if(length(singles) > 0) {			# If still singles left, make more pairs
						sing <- matrix(subset(population, population[,1] %in% singles), ncol=13)
						more <- pairoff(sing[,1], sing[,4], Inf)
						pairs <- matrix(rbind(pairs, more$mated),ncol=2)
						singles <- singles[!(singles %in% more$mated)]
					}
				}
			}
		}
		
		## UPDATE LIVE INDIVIDUALS; # YEARS ALIVE, # YEARS BRED FOR EACH	
		singles <- singles[!(singles==0)]
		liveadults <- unique(c(pairs, singles)) 
		if(nrow(population) > 0){
			# Drop dead individuals if specified:
			if(!trackall) population <- matrix(population[which(population[,1] %in% liveadults),], ncol=13)
			numpresent <- length(liveadults)
			if(!KAdults) numpresent <- numpresent + nrow(nonbreeders)	# nonbreeders are included in K
			# Add 1 year to col 11 in population/nonbreeders for all individuals alive (in nonbreeders or liveadults):			
			if(nrow(nonbreeders) > 0) {
				tally <- c(nonbreeders[,11] + 1)
				nonbreeders <- matrix(cbind(nonbreeders[,1:10, drop=FALSE], tally, nonbreeders[,12:13, drop=FALSE]), ncol=13)
			}			
			popalive <- matrix(population[which(population[,1] %in% liveadults),], ncol=13)
			tally <- c(popalive[,11] + 1)	
			popalive <- matrix(cbind(popalive[,1:10, drop=FALSE], tally, popalive[,12:13, drop=FALSE]), ncol=13)
			popdead <- matrix(population[-(which(population[,1] %in% liveadults)),], ncol=13)
			population <- rbind(popdead, popalive)	
			# Add 1 year bred to col 12 for individuals in pairs:
			popprs <- matrix(population[which(population[,1] %in% pairs),], ncol=13)
			if(nrow(popprs) > 0) {
				tally <- popprs[,12] + 1
				popprs <- matrix(cbind(popprs[,1:11, drop=FALSE], tally, popprs[,13, drop=FALSE]), ncol=13)	
				popunpr <- matrix(population[-(which(population[,1] %in% pairs)),], ncol=13)
				population <- matrix(rbind(popprs, popunpr), ncol=13)
			}
		}


		## CENSUS	
		adults <- unique(c(pairs,singles))		# list of IDs of all live adults
		allID <- c(adults, nonbreeders[,1])	# unique IDs of all live individuals
		if(nrow(nonbreeders) == 0) nnb <- 0
		else nnb <- nrow(nonbreeders)
		if(length(pairs) == 0) npr <- 0
		else npr <- nrow(pairs)
		if (length(allID)==0) break
		census[y,1] <- length(adults)	          	# number of adults	
		census[y,2] <- length(unique(c(pairs[,1])))  # number of breeding females
		census[y,3] <- length(unique(c(pairs[,2])))	# number of breeding males
		popprs <- matrix(population[which(population[,1] %in% unique(c(pairs))),], ncol=13)
		popalive <- matrix(population[which(population[,1] %in% adults),], ncol=13)	# live adults
		if(GeneCount == "adult"){
			census[y,4] <- sum(popprs[,5:6])  	# number of rare alleles in breeding adult population
			counted <- popprs
		}	
		else{
			census[y,4] <- sum(popalive[,5:6]) + sum(nonbreeders[,5:6]) 	# number of rare alleles in whole population
			counted <- rbind(popalive, nonbreeders)
		}	
		census[y,5] <- npr               			# number of pairs
		census[y,6] <- nnb			# number of nonbreeders
		census[y,7] <- nrow(subset(popalive,popalive[,7]==1 | popalive[,7]==2)) + nrow(subset(nonbreeders, nonbreeders[,7]==1 | nonbreeders[,7]==2))	# number of starters remaining, including supplementals (from addN)
		census[y,8] <- nrow(subset(popalive,popalive[,7]==4)) + nrow(subset(nonbreeders, nonbreeders[,7]==4))	# number of migrants in population
		census[y, 9] <- mean((y+1) - popprs[,8])	# mean age of breeding adults in population NEXT yr	
		
		# translate NAs into 0s (for summary purposes)
		for(i in 1:nrow(census)){
			for(j in 1:ncol(census)){
			if(is.na(census[i,j])) census[i,j] <- 0
			}
		}
		
		if(printplots){
			if(KAdults) nindiv <- census[,1]
			else nindiv <- c(census[,1] + census[,6])
			afreqs <- c(census[,4]/(nrow(counted)*2))	
			plotlim <- q0*5
			if (plotlim > 1) plotlim <- 1
			if(r == 1)if(y==1){	# set up empty plots
				split.screen(c(2,1))
				close.screen(n=c(1,2), all.screens=TRUE)
				split.screen(c(2,1))
				screen(1)
				plot(nindiv ~ c(1:nyears), type="n", ylim=c(0, K + K*0.05), ylab="number of individuals", xlab = "year", main = "Population growth")
				screen(2)
				plot(afreqs ~ c(1:nyears), type="n", ylim=c(0, plotlim), ylab="frequency of rare allele", xlab = "year", main = "Allele frequency")
			}
			if(y == nyears){
				screen(1)
				plot(nindiv ~ c(1:nyears), type="l", ylab=" ", xlab=" ", ylim = c(0, K + K*0.05), xaxt="n", yaxt = "n", main = NULL)		
				screen(2)
				plot(afreqs ~ c(1:nyears), type="l", ylab=" ", xlab=" ", ylim=c(0, plotlim), xaxt="n", yaxt = "n", main = NULL)
			}
		}
		
		
		## INDIVIDUAL INFORMATION - from last year of simulation only (contains all individuals from all years, if trackall = TRUE)
		if(y == nyears) {
			totalpop <- rbind(population, nonbreeders)		
			indivdata <- totalpop[,c(1,2,3,7,8,11,12,13)]	# only keep relevant information (ncol=8)		
		}

	}		#  End year loop
	out[[1]] <- census
	out[[2]] <- indivdata
	rm(census, indivdata, population, nonbreeders)
	if(is.wholenumber(r/nreplprint)) cat("Replicate #", r, "finished at", paste(Sys.time()), sep=" ", "\n")
	flush.console()
	out
  }			#  End replicate function
	
  #############
  
## RUN EACH REPLICATE:
result <- lapply (1:nrepl, run.one.repl)  
# This is a list of sub-lists... one sub-list for each replicate, containing 1) census summary matrix and 2) indiv data from last year of the replicate.
result
}		## End aRetain function

	

######################################################################################

aRetain.summary <- function (adata, GeneCount, alpha=0.05, dropextinct = FALSE) {
## summarize output from aRetain (object called adata) across replicates, for an overall average for each year of the simulation.
## GeneCount should be as specified in aRetain.
## use alpha = 0.05 for 95% confidence limits
	SE <- function(x) sqrt(var(x, na.rm=TRUE)/sum(!is.na(x)))
	# pull out just the census information (drop indivdata), from adata[[r]][[1]] where r = each replicate
	nrepl <- length(adata)
	nyears <- dim(adata[[1]][[1]])[1]
	censusdata1 <- array(dim = c(nyears, 9, nrepl))
	for(r in 1:nrepl){
		censusdata1[,,r] <- matrix(adata[[r]][[1]], ncol=9)
	}
	if (dropextinct){
		# all data will be summarized only for replicates in which the population did not go extinct.  P.extant (proportion of replicates in which population did not go extinct) and P.xLCL/P.xUCL are the only summaries that include replicates in which the population did go extinct.
		censusdata <- censusdata1[,,(censusdata1[nyears,1,] > 0)]
		if (is.matrix(censusdata)){
			print(censusdata)
			stop('Error: The population survived to the end of the period in only one replicate (shown above); this output is not averaged over multiple replicates.  Multiple replicates are needed to obtain output from aRetain.summary.')
		}	
		cat("Note: The population persisted to the end of the simulation (did not go extinct) in", dim(censusdata)[3], "of", nrepl, "replicates. Information is summarized for these replicates.\n")
	}
	else{
		censusdata <- censusdata1	# ALL replicates will be included in all of the above summaries.
		cat("Note: All replicates, including those in which the population went extinct, are included in the summary.  This will affect the estimates of population size/composition, mean age, and probability of allele retention given here.\n")
	}
	nyears2 <- dim(censusdata)[1]
	nrepl2 <- dim(censusdata)[3]
	num <- censusdata[,1,]	# number of adults present
	if(GeneCount == "all"){
		numnb <- censusdata[,6,]	# number of nonbreeders present
		numall <- matrix(nrow=nyears2, ncol=nrepl2)
		for(i in 1:nyears2){
			for(j in 1:nrepl2){
				numall[i,j] <- numnb[i,j] + num[i,j]
			}
		}
		num <- numall	# all individuals were counted toward allele retention
	}
	else if(GeneCount != "adult") stop('GeneCount must equal "adult" or "all".')
	alleles <- censusdata[,4,]
	freqs <- matrix(nrow = nyears2, ncol = nrepl2)
	for(i in 1:nyears2) {
		for(j in 1:nrepl2) {
			freqs[i,j] <- alleles[i,j]/(2*num[i,j])
		}
	}
	
	MeanNAd   <- apply(censusdata[,1,], 1, mean, na.rm = TRUE)
	SEN      	<- apply(censusdata[,1,], 1, SE)
	n        	<- apply(censusdata[,1,], 1, function(x) sum(!is.na(x)))
	P.retain 	<- apply(censusdata[,4,]>0, 1, mean, na.rm = TRUE)
	P.LCL    	<- Blcl(P.retain, n, alpha)
	P.UCL    	<- Bucl(P.retain, n, alpha)
	n1        	<- apply(censusdata1[,1,], 1, function(x) sum(!is.na(x)))
	P.extant 	<- apply(censusdata1[,1,]>0, 1, mean, na.rm = TRUE)
	P.xLCL   	<- Blcl(P.extant, n1, alpha)
	P.xUCL   	<- Bucl(P.extant, n1, alpha)
	MeanBrF  	<- apply(censusdata[,2,], 1, mean, na.rm = TRUE)
	SEBrF    	<- apply(censusdata[,2,], 1, SE)
	MeanBrM  	<- apply(censusdata[,3,], 1, mean, na.rm = TRUE)
	SEBrM    	<- apply(censusdata[,3,], 1, SE)
	MeanNNonbr   <- apply(censusdata[,6,], 1, mean, na.rm = TRUE)
	MeanNFound	<- apply(censusdata[,7,], 1, mean, na.rm = TRUE)
	MeanNMigr	<- apply(censusdata[,8,], 1, mean, na.rm = TRUE)
	MeanAge 	<- apply(censusdata[,9,], 1, mean, na.rm = TRUE)
	A.Freq 		<- apply(freqs, 1, mean, na.rm = TRUE)
	A.SE		<- apply(freqs, 1, SE)
	cbind (MeanNAd, SEN, MeanNNonbr, MeanBrF, SEBrF, MeanBrM, SEBrM, MeanNFound, MeanNMigr, MeanAge, P.extant, P.xLCL, P.xUCL, P.retain, P.LCL, P.UCL, A.Freq, A.SE)
}
########################################################################

indiv.summary <- function(adata, genlength, alpha=0.05){
## Summarize individual data.  Input:
# adata: object returned by aRetain
# genlength: mean age of simulated breeding individuals, after stabilizing from any founder age effects (as derived from census output).
# alpha: significance level for confidence limits
	nrepl <- length(adata)
	nyears <- dim(adata[[1]][[1]])[1]
	indivOUT <- array(dim = c(4, 8, nrepl))
	colnames(indivOUT) <- c("n", "pbreed", "pbreed.LCL", "pbreed.UCL", "YrsBred", "YrsBredBr", "lifespan", "effectivegen")
	rownames(indivOUT) <- c("starters", "supplements", "locals", "migrants")
	
	for(r in 1:nrepl) {
		me <- matrix(adata[[r]][[2]], ncol=8)	# individual data
			# cols: ("ID", "dam", "sire", "origin", "birthyr", "NYrAlive", "NYrBred")
		if(nrow(me) > 0){	# population did not go extinct
			for (i in 1:4){		# for each origin
				chosen <- matrix(me[me[,4]==i,], ncol=8)	# indiv of origin i	
				chosen <- matrix(chosen[which(chosen[,8] < (nyears - genlength)),], ncol=8)	# take out those added during the last generation (unequal chance of breeding)				
				indivOUT[i,1,r] <- nrow(chosen)	# total number of origin i
				bred <- matrix(chosen[chosen[,7] > 0,], ncol=8)
				pbred <- round(nrow(bred) / nrow(chosen), digits=2)	# proportion that bred at some point				
				indivOUT[i,2,r] <- pbred
				indivOUT[i,3,r] <- Blcl(pbred, nrow(chosen), alpha)
				indivOUT[i,4,r] <- Bucl(pbred, nrow(chosen), alpha)
				indivOUT[i,5,r] <- round(mean(chosen[,7]), digits=2)	# mean num breeding seasons each				
				indivOUT[i,6,r] <- round(mean(bred[,7]), digits=2)		# mean num breeding seasons for those that DID breed
				indivOUT[i,7,r] <- round(mean(chosen[,6]), digits=2)	#  mean lifespan
				effgen <- (nrow(chosen)*nrow(bred) / nrow(chosen)) / (nyears/genlength)	# average number of effective individuals per generation
				indivOUT[i,8,r] <- effgen
				if(i < 3) indivOUT[i,8,r] <- NA		# founders not added in every generation so don't calculate this			
			}
		}		
		else indivOUT[,,r] <- matrix(NA, nrow=4, ncol=8)	# population had already gone extinct
	}
	
	### AVERAGE EACH CELL ACROSS REPLICATES
	indivinfo <- matrix(nrow=4, ncol=8)
	for(j in 1:nrow(indivOUT)){ for(c in 1:ncol(indivOUT)){
		indivinfo[j, c] <- mean(as.numeric(indivOUT[j,c,]), na.rm=TRUE)	# averaged across replicates (with any NA values, from replicates in which population went extinct, removed)
	}}
	colnames(indivinfo) <- c("n", "pbreed", "pbreed.LCL", "pbreed.UCL", "YrsBred", "YrsBredBr", "lifespan", "effectivegen")
	rownames(indivinfo) <- c("starters", "supplement", "locals", "migrants")	
	indivinfo
}	


########################################################################


pedigree.summary <- function(adata){
## Get mean and variance of inbreeding coefficient for each year, across replicates
## Requires package pedigree (which requires Matrix, lattice, HaploSim, reshape)
	
	require(pedigree)
	nrepl <- length(adata)
	nyears <- dim(adata[[1]][[1]])[1]
	outF <- array(dim = c(nyears, 2, nrepl))
	for(r in 1:nrepl) {
		me <- matrix(adata[[r]][[2]], ncol=8)	# individual data (ID, dam, sire, origin, birthyr, nyrsalive, nyrsbred)
		if(nrow(me) > 0){	# population exists
			IDdamsire <- data.frame(me[,1], me[,2], me[,3])
			colnames(IDdamsire) <- c("ID", "dam", "sire")
			if(nrow(IDdamsire) > 0){
				ord <- orderPed(IDdamsire) # Order to work with new requirements for calcInbreeding
				IDdamsire <- IDdamsire[order(ord),]
				Fvalues <- calcInbreeding(IDdamsire)	# calculate F of each individual
			}
			else Fvalues <- numeric(0)
			meF <- matrix(cbind(me, c(Fvalues)), ncol=9)	# attach F values to individual information
			avgF <- rep(NA, nyears)
			varF <- rep(NA, nyears)
			for(y in 1:nyears) {	# subset the individuals alive during that year
				added <- subset(meF, meF[,8] <= y)	# individuals added before or during year y
				alive <- subset(added, (added[,8] + added[,6]) > y)	# individuals that were still alive at end of year y			
				if(nrow(alive) > 0){
					avgF[y] <- mean(alive[,9], na.rm=TRUE)	# average Fs across individuals alive now
					varF[y] <- var(alive[,9], na.rm=TRUE) 	# across-individual variance in F for each year
				}
			}
			outF[,,r] <- cbind(avgF, varF)
		}
		else outF[,,r] <- matrix(NA, nrow=1, ncol=2)	# population had already gone extinct
	}
	avgF2 <- rep(NA, nyears)
	varF2 <- rep(NA, nyears)
	indvarF2 <- rep(NA, nyears)
	for(y in 1:nyears){
		avgF2[y] <- mean(outF[y,1,], na.rm=TRUE)	# average F across replicates
		varF2[y] <- var(outF[y,1,], na.rm=TRUE)	# variance of F across replicates
		indvarF2[y] <- mean(outF[y,2,], na.rm=TRUE)	# mean inter-individual variance of F (across replicates)
	}
	
	pedsum <- cbind(avgF2, varF2, indvarF2)
	yrs <- 1:nrow(pedsum)
	pedsum <- cbind(yrs, pedsum)
	colnames(pedsum) <- c("year", "meanF", "varF", "indivVarF")
	pedsum
}
	
########################################################################
	
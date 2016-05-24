###############################################################################
#
# FUNCTION: coalsim, coalsimC
#
# This function utilizes Hudsons ms[1] program to generate samples under the 
# neutral model. It then uses its output to calculate a variety of summary
# statistics by mean of coalescent theory.
# 
# For each loci the ms program is invoked with the user given conditions
# (i.e. mutation rate, recombination, population growth, demography, exchange 
# gene conversion etc.) and its output is then read in and necessary 
# parts of it are saved in order to perform so called summary statistics.
#
# Currently the following stats are collected:
#	- PI
#	- Fay & Wu's ThetaH
#	- Tajima D
#	- Fu and Li's D (1993)
#	- Fu and Li's F (1993)
#	- Fay and Wu H
#
# These computed stats are then used to derive summarized statistics for each 
# loci across all iterations. If a matrix with observed values is passed on as
# an argument, those values are compared to simulated ones giving P-values about
# how those ovserved values could have arised in the past. 
#
# The function returns an object of class csstats, it contains some general 
# values about the condition under which the simulation was performed. It 
# provides some summarized stats for all loci and also a single matrix with all
# p-values across all loci. Furthermore it contains a list of locstats objects,
# which gives detailed statistics about each test in each iteration for each 
# loci
# 
# The schematic structure of the return object:
#   ___________
#  |  csstats  |
#  |  =======  |
#  |  .		   |
#  |  .		   |
#  |___________|		 __________			 __________			 __________
#  |  locus	   | ----> 	| locstats | ---->  | locstats | ---->  | locstats |
#  |___________|		| ======== |		| ======== |		| ======== |
#						|		   |		|		   |		|		   |
#						|		   |		|		   |		|		   |
#						|__________|		|__________|		|__________|
#
#
#
#
# [1] Hudson, R. R. (2002) Generating samples under a Wright-Fisher neutral
# model of genetic variation. Bioinformatics 18: 337-338.
#
#
# FUNCTION CALLS: 	init_coef
#					probabilties
#					neuttest
#					progressBar
#
# PARAMETERS:
#		nsam:		is the number of copies of the locus in each sample. It needs 
#					to be provided as a vector of length nloci 
#					length( c(x,y,z) ) = nloci
#
#		niter:		number of independant samples to generate for each locus
#					single integer value greater than 0
#
#		theta:		mutation parameter theta (4Nmu), where N is the diplod  
#					population size and mu the mutation rate per locus.
#					It needs to be provided as a vector of length nloci 
#					length( c(x,y,z) ) = nloci
#
#		nloci:		number of loci, single integer value greater than 0
#
#		npop:		number of populations from which observed values originate
#					single integer value greater than 0
#
#		nsites:		number of nucleotid sites for each locus, if provided, 
#					must be a vector of length nloci.
#		
#		seeds:		specify 3 random number seeds. a vector of length 3 with 
#					positive values is expected
#
#		obsVal:		a matrix with observed values to test against.
#					needs to be provided as a matrix of size numTests x nloci
#					each row consisting of one locus and each column denote the
#					value for a specific test. The order of the tests are specified 
#					in the vector testNames (below)
#		
#		printtree:	set to 1 to include tree information in class, which can later 
#					be printed using the genetree function
#					printing tree is not available with recombination and geneConv
#
#	fixedSegsites:	usually the number of segregating sites varies in each iteration. 
#					Please provide a  single numeric value if the number of 
#					segregating sites needs to be fixed.
#
#	recombination:	provide a vector of format: c(p, nsites) 
#					p = cross over parameter rate, nsites is the number of sites
#					between recombination occurs
#	
#		geneConv:	in addition to recombination intra-locus non-cross-over 
#					exchange gene conversion can be included in simulation
# 					expected format is c(f, gamma) 
# 					f denote the ratio, g/r, where r is the probability per generation 
#					of crossing-over between adjacent sites. (see Wiuf and Hein 2000)
# 					gamma is the mean conversion tract length			
#
#		growth:		population size is measured by N(t) = N0 exp-^(alpha*t). provide alpha 
#					as integer value. negative values indicate that population was larger 
#					in the past than present
#
#	migration:		specify a vector of length 'npop', each element denoting the 
#					number sampled regarded as to belong to a certain subpopulation.  
#
#
#	demography:		vector of length 3 or 4 with first value denoted as 'type' 
#
#					valid 'types' for vectors of length 3 are as following: 
#					- 1 to set a growth rate change alpha at a certain time t: 
#					  c(1, t, alpha)
#
#					- 2	set all subpop to size x * N_0 and growth rate to zero: 
#					  c(2, t, x)
#
#					- 3 set all elements of migration matrix to x/(npop-1): 
#					  c(3, t, x)					
#
#					valid 'types' for vector of length 4 with the following values:	
#					- 4 set growth rate of subpop i to alpha at time z: 
#					  c(4, t, i, alpha)
#
#					- 5 set subpop i size to x * N_0 at time t and growth rate to zero: 
#					  c(5, t, i, x)
#
#					- 6 split subpopulation i into subpopulation i and a new subpopulation, 
#					  labeled npop + 1. Each ancestral lineage in subpopulation i is randomly  
#					  assigned to subpopulation i with probability p and subpopulation  
#					  npop + 1 with probability 1 - p. The size of subpopulation npop + 1 is  
#					  set to N_0. Migration rates to and from the new subpopulation are assumed  
#					  to be zero and the growth rate of the new subpopulation is set to zero: 
#					  c(6, t, i, p)
#
#					- 7 move all lineages in subpopulation i to subpopulation j at time t. 
#					  Migration rates from subpopulation i are set to zero: 
#					  c(7, t, i, j)
#
#
#
# RETURN VALUES:	an object of class 'csstats', available Slots
#					- probLess: 	Prob. that simulated val. <= to observed val. P(Sim <= Obs)
#					- probEqual:	Prob. that simulated val = to  observed val. P(Sim = Obs)
#					- validIter:	number of valid iteration for each test and loci
#					- obsVal:		observed values for each test (input)
#					- nloci:		number of loci considered for this simulation
#					- niter:		number of iterations for each loci for this simulation
#					- average:		average values of each statistic calculated across all loci
#					- variance:		variance values of each statistic calculated across all loci
#					- locus:		a list of locstats objects, containing detailed stats for each locus
#
# AUTHOR:			Niroshan Nadarajah <niroshan.nadarajah@uni-duesseldorf.de>
#
# LAST MODIFIED:	10/11/04
#
###############################################################################

#------------------------------------------------------------------------------#
#							Definition of some constants					   #
#------------------------------------------------------------------------------#
#specify number of tests performed for each loci
#numTests  <- 11
#specify the names of the column header for matrices
#testNames <- c("Tajima.D","n.segregating.sites","Rozas.R_2","Fu.Li.F","Fu.Li.D","Fu.F_S","Strobeck.S","Fay.Wu.H","Zeng.E","theta_Tajima","theta_Watterson")

#specify the quantile point intervals. Values between 0 and 1 are accepted
quantileProbs = c(0.001, 0.010, 0.025, 0.050, 0.100, 0.250, 0.500, 0.750, 0.900, 0.950, 0.975, 0.990, 0.999)

#change this to 1 to include position information of each halpotypes into the locus stats class, 0 otherwise
preservePositionInfo <- 1

#change this to 1 to include position information of each halpotypes into the locus stats class, 0 otherwise
preserveHaplotypeInfo <- 0


#------------------------------------------------------------------------------#
#                            	 coalsimC									   #
#------------------------------------------------------------------------------#
## This is a convenience function, which is invoked with just an object of class 
## testparams, avoiding to pass on every single parameter to coalsim() directly.  
coalsimC <- function(p,detail=FALSE,testNames,numTests,neutrality=FALSE,linkage=FALSE,F_ST=FALSE,MSMS=FALSE,big.data){
	
	#check validity of arguments
	if (class(p)[1] != "test.params")
		stop("this function needs an object of class 'test.params'")
	
	
	if (p@n.sam > 0 && p@n.iter > 0 && length(p@theta) > 0 ) {
		
		nloci <- p@n.loci
		if ( length(p@n.loci) == 0 )
			nloci = 1
		
		npop <- p@n.pop	
		if ( length(p@n.pop) == 0 )
			npop = 1
		
		nsites <- p@n.sites	
		if ( length(p@n.sites) == 0 )
			nsites = NA
		
		obsVal <- p@obs.val	
		if ( length(p@obs.val) == 0 )
			obsVal = NA
		
		printtree <- p@print.tree
		if ( length(p@print.tree) == 0 )
			printtree = 0
		
		fixedSegsites <- p@fixed.seg.sites
		if ( length(p@fixed.seg.sites) == 0 )
			fixedSegsites = NA
		
		recombination <- p@recombination
		if ( length(p@recombination) == 0 )	
			recombination <- NA
		
		geneConv <- p@gene.conv	
		if ( length(p@gene.conv) == 0 )	
			geneConv <- NA
		
		growth <- p@growth
		if ( length(p@growth) == 0 )		
			growth <- NA
		
		demography <- p@demography
		if ( length(p@demography) == 0 )		
			demography <- NA
		 
                migration <- p@migration
		if ( length(p@migration) == 0 )		
			migration <- NA
		 
                seeds <- p@seeds
		if ( length(p@seeds) == 0 )		
			seeds <- NA
		
		return(coalsim(p@n.sam, p@n.iter, p@theta, nloci, npop, nsites, obsVal, printtree, fixedSegsites, recombination, geneConv, growth,demography,seeds,migration,detail,testNames,numTests,neutrality,linkage,F_ST,MSMS,big.data))
		
	} else {
		stop("not enough parameters set to do the simulation")	
	}
}



#------------------------------------------------------------------------------#
#                            	 coalsim									   #
#------------------------------------------------------------------------------#
coalsim <- function(nsam=c(10), niter=2, theta=c(5.0), nloci=1, npop=1, nsites= NA, obsVal = NA,  printtree = 0, fixedSegsites = NA, recombination=NA, geneConv=NA, growth=NA, demography=NA, seeds=NA, migration=NA,detail,testNames,numTests,neutrality,linkage,F_ST,MSMS,big.data){

	# do some error checking/handling prior to any calculation
	if( !is.numeric(nsam) )
		stop("Please provide a list of integer values for nsam")
		
	#if (length(nsam) != nloci)
	#	stop(paste("The number of loci (", nloci ,") and the number of samples for each loci (", length(nsam) ,") mismatch. Please provide a vector of correct length as the first parameter"))
		
	if (niter < 1)
		stop("The number of iteration is invalid. Please provide a positive value greater 0")
		
		
	if (length(theta) != nloci)	
		stop(paste("The number of loci (", nloci ,") and the number of value of theta for each loci (", length(theta) ,") mismatch. Please provide a vector of correct length as the third      parameter"))
	
	if (nloci < 1 && nloci != -1)
		stop("The number of sites is invalid. Please provide a positive value greater 0.")
	
	
	if ( !is.na(nsites) && nsites < 1)
		stop("The number of sites is invalid. Please provide a positive value greater 0.")
	
	# check for validity of the observed value matrix
	if ( !is.na(obsVal[1][1])  ) {
		if ( dim(obsVal[[1]])[1] != nloci && dim(obsVal[[1]])[2] != numTests )  
			stop(paste("The dimensions of the observed values matrix are incorrect. For (", nloci , ") loci you need to specify a ", numTests ,"x", nloci ," matrix "))
	}
	
	if (printtree != 0 && printtree != 1)
		stop("Accepted values for parameter printtree are 0 or 1")
	
	if (class(nsam)[1] == "testparams")
		stop("Please use function coalsimC in combination with a 'testparams' object")
	
	
	# collate structure for ms parameters according to user input
	param <- ""
	
	# check whether recombination is in format c(f, gamma)
	if ( !is.na(recombination)[1] ) {
		if ( length(recombination) == 2 ) {
			param <- paste(param, " -r", recombination[1], recombination[2])
			
			# printing tree is not available with recombination parameter....omit option
			if (printtree == 1) {
				print("Print tree option is not available with recombination. Option omitted!")
				printtree = 0
			}
			
		} else {
			stop("Please provide a vector of length 2 for recombination!")
		}
	}
	
	# check whether geneConv is in format c(f, gamma)
	if ( !is.na(geneConv)[1]) {  
		if( length(geneConv) == 2 ) {
			param <- paste(param, "-c", geneConv[1], geneConv[2])
			
			# printing tree is not available with recombination parameter....omit option
			if (printtree == 1) {
				print("Print tree option is not available with gene conversion. Option is omitted!")
				printtree = 0
			}
		} else {
			stop("Please provide a vector of length 2 for specifing the gene conversion rate. ")
		}
	}
	
	
	if ( printtree == 1 ) {
		param <- paste(param, "-T")
	}
	
	if ( !is.na(fixedSegsites)) {
		if (is.numeric(fixedSegsites) && fixedSegsites > 0){
			param <- paste(param, "-s ", fixedSegsites)
		}else{
			stop("The number of segregating sites is invalid")
                }	
        }
	
	
	if ( !is.na(growth) ) {
		param <- paste(param, " -G ", growth)
	}
	
	if ( !is.na(demography)[1] ) {
		dgParam<- getDemographyParam(demography)
		param <- paste(param, dgParam)
	}
	
	if ( !is.na(seeds[1])) {
		if ( length(seeds) != 3 ) {
			stop("The 'seeds' vector needs to specify a vector of length 3.")
		} else {
			param <- paste(param, "-seeds", seeds[1], seeds[2], seeds[3])
		}
	}
	

	if (!is.na(migration[1])) {
		
		if (npop < 2) {
			stop("migration assumes al least a population size of 2")
		}
		
		#if (npop != length(migration) - 1) {
		#	stop("Please assign a vector of length 'npop' + 1, specifing the number of loci per population. The last item in the vector is assumed to be the mutation parameter")
		#}
		
		nsubpoploci = c(0)
		mig <- ""
		for(i in 2:length(migration)-1){
			nsubpoploci = nsubpoploci + migration[i]
			
			mig <- paste(mig, migration[i])
		}
		
		#if (nsubpoploci != nsam) {
	#		stop("the sum of the loci of all subpopulation in the migration vector must match the number of samples")
	#	}
		
		param <- paste(param,"-I", npop, mig, migration[length(migration)])
	}
	
	 
	# call function to show a progressbar during computation
	# progr <- progressBar()
	
	
	############################################################## ----------------------------- 
  # Definiere Populationen (nur wenn mehr als eine Population)
  # ----------------------------------------------------------- # ------------------------------ 
  # nsam ist hierbei L\E4nge der Populationen WICHTIG !!!!
  

  for(cloci in 1:nloci) { # -------------------------------------------------------- 
	
	   ## Define Populations
	     mymig <- migration[1:npop]
   if(npop>1){
   
     Populations      <- vector("list",npop)
     ids <- 1:nsam[cloci] #<----------------------------------------------------------------------

      for(xx in 1:npop){
       Populations[[xx]] <- ids[1:mymig[xx]]
       if(npop!=xx){ids  <- ids[(mymig[xx]+1):length(ids)]}
      }

   }else{ Populations <- list(1:nsam[cloci])}
      #################	
	
	         ### Fot the package
                 # p_path <- .path.package("PopGenome")
		
  #              if (.Platform$OS.type == "unix") {
			# mscall <- paste(p_path,"/exec/ms",sep="")
	#		  mscall <- file.path(getwd(),"/ms",sep="")
	#	} else {
	#		if (.Platform$OS.type == "windows") {
				# mscall <- paste(p_path,"\\exec\\ms.exe",sep="")
  #       mscall <- file.path(getwd(),"ms.exe",sep="")
	#		} else {
	#			stop("Your platform is unsupported.")
	#		}
	#	}

	
	
		# visually update progress
		#  progr <- progressBar(cloci, nloci, progr)
		
		# execute ms and get the output
		#if (.Platform$OS.type == "unix") {
		#   	mscall <- paste(getwd(),"/ms",sep="")
                #        "res/ms "	
		#} else {
	#		if (.Platform$OS.type == "windows") {
	#			mscall <- paste(getwd(),"\\ms.exe",sep="")
         #                       # "res\\ms.exe "
	#		} else {
	#			stop("Your platform is unsupported.")
	#		}
	#	}

              
		if (.Platform$OS.type == "unix"){ 
		   mscall <- paste("./ms", nsam[cloci], niter, " -t",theta[cloci], param)#-------------------------------------- call ms !!!!!!
                }
		
                if (.Platform$OS.type == "windows"){
 		   mscall <- paste("ms.exe", nsam[cloci], niter, " -t",theta[cloci], param)
                } 
                
     
         
	   # ---------------------------------------------------------------------------------------------------------       
         # just a dummy insertion to calculate msms
         
		   MSMS_in <- FALSE
		   
                old_workspace  <- getwd()

                if(MSMS[1]!=FALSE){			
	                 ss       <- setwd(file.path(old_workspace,"msms","bin"))

			 if (.Platform$OS.type == "unix"){ 	
 	                  mscall   <- paste("./msms",MSMS,"-ms",nsam[cloci], niter, " -t", theta[cloci], param)
		         }
			
			 if (.Platform$OS.type == "windows"){
  			  mscall   <- paste("msms.exe",MSMS,"-ms",nsam[cloci], niter, " -t", theta[cloci], param)
                         } 
	
			             MSMS_in  <- list(nsam=nsam[cloci])
			                
		}
		             
	#-----------------------------------------------------------------------------------------------------------

	 
	 
	 
	 

       #         output <- system(mscall, intern=TRUE)
  	  	# process/save ms input parameters (i.e. first line)
  #		line1 <- unlist(strsplit(output[1], " "))
		#print(paste("Processing Loci", cloci, "..."))
	
		# process/save seed values (i.e. the second line of output)
#		line2 <- unlist(strsplit(output[2], " "))
#		seeds <- c(as.numeric( line2[1] ), as.numeric( line2[2] ), as.numeric( line2[3] ) )
	  
#		nlines <- length(output)
	
#		# define some variables to accumulate collected data
#		first <- TRUE
#		positions <- matrix()
#		trees <- matrix()
#		
 #   for (i in 1:nlines) {
			# look for next sample data output
#			
#			if (output[i] == "//") {
#				
#				# new sample detected, move to next line
#				i = i+1
#				
##				# check whether gene tree information (in Newick format) are present, i.e. printtree option is active
#				
#				if (substr(output[i], 0, 1) == "(") {
##					
#					if (first){
#						trees <- rbind(output[i])
##					} else {
#						trees <- rbind(trees, output[i])
#					}
#					
##					i = i+1
#				} 
##				
#				# get the line with the number of segregated sites
#				segsitesLine      <- unlist(strsplit(output[i], " "))
#				nsegsites	  <- as.numeric(segsitesLine[2])
##				
#				
 #       # in case of option fixedSegsites another line is omitted giving the propability emission of this number of
##				# segregated sites, this information is omitted for now
#				if (substr(output[i], 0, 1) == "p") {
#					i = i+1
##				} 
#								
#				# haplotypes <- matrix()
#				
##				
   



 #   if ( nsegsites > 0 ) {  	
  #      #       # move along to line with position information
#					i = i+1
#					
#					if (preservePositionInfo) {
#						# get position information
#						# positionsLine <- unlist(strsplit(output[i], " "))
#						# positions <- rbind(positions, as.numeric(positionsLine[2:length(positionsLine)]))
##					
#						if (first){
#							positions <- rbind(substr(output[i], 12, nchar(output[i])-1))
##						} else {
#							positions <- rbind(positions, substr(output[i], 12, nchar(output[i])-1))
##						}
#					}
##					
#					# move on to the next line, i.e the first line of haplotypes
#					i = i+1
##				
#					for (j in 1:nsam[cloci]) {
##
#						row <- as.numeric(unlist(strsplit(output[i], "")))
#				
##						# aggregate all haplotypes in a matrix
#						if (j==1) {
#							haplotypes <- rbind(row)
##						} else {
#							haplotypes <- rbind(haplotypes, row)
##						}
#				
##						# increase i so that this line is not visited again once more
#						i = i+1
##					}
#			
 #                                                  # haplotypes is the biallelic matrix
#					           # analyse sample
##					           print(haplotypes)                           





			 	if (.Platform$OS.type == "unix") {
			  system(paste(mscall," > ms.out"))
		} else {
			if (.Platform$OS.type == "windows") {
			  shell(paste(mscall," > ms.out")) 
			} else {
				stop("Your platform is unsupported.")
			}
		}
			
			
			    if(!big.data){		
  			    msout <- read.ms.output(file.ms.output="ms.out",MSMS=MSMS_in)
                            }
			    if(big.data){
			    msout <- read.big.ms.output("ms.out")
			    }		

			    setwd(old_workspace)	

			for(zz in 1:length(msout$segsites)){

				trees <- matrix()
                                seeds <- NaN	
				if(msout$segsites[zz]>0){ # wenn segsites da

					 if(big.data){
					 open(msout$gametes[[zz]])  
					 }
                                         haplotypes <- msout$gametes[[zz]][,,drop=FALSE] 
					 if(big.data){
					 close(msout$gametes[[zz]])  
					 }						
					 
					 positions  <- msout$positions[[zz]]
 					 #print(haplotypes)

                                             ## STATISTICS ########################################
                                               if(neutrality){

                                                 obj      <- calc_freqstats_FAST(haplotypes,Populations)
                                                 SEG      <- obj$THETA[1,]
                                                 thetaT   <- obj$THETA[3,]
                                                 thetaS   <- obj$THETA[2,]
                                                 TD       <- obj$taj_D
                                                 # R2       <- calcR2(haplotypes,Populations,thetaT,SEG)
                                                 FuLi_F   <- obj$FuLi_F
                                                 FuLi_D   <- obj$FuLi_D
                                                 HnFw     <- obj$HnFw
                                                 Ez       <- obj$Ez
                                                 nix      <- rep(NaN,npop)
                                                 FS       <- nix
                                                 Strobeck <- nix
                                                 R2       <- nix

                                               if(detail){

                                                obj2     <- calc_FS(haplotypes,Populations,thetaT) 
                                                FS       <- obj2$FS
                                                Strobeck <- obj2$Strobeck
                                                
                                               }

                                              }

                                              if(linkage){

                                                res              <- wall99bq(haplotypes,Populations)
                                                WALLB            <- res$B 
                                                WALLQ            <- res$Q
                                                 nix             <- rep(NaN,npop)
                                                 Zns             <- nix
                                                 ZA              <- nix
                                                 ZZ              <- nix 
                                        
						if(detail){
                                                
                                                 res             <- linkdisequ(haplotypes,Populations)
    						 Zns             <- res$Zns
  						 ZA              <- res$ZA
     						 ZZ              <- res$ZZ

                                                }


                                              }
                                              
                                              if(F_ST){
                                              
                                                res             <- calc_hwhafsth(haplotypes,Populations,simulation=TRUE)
                                                hapw            <- res$hapw
                                                Pi              <- res$PIW_nei
                                                FSTH            <- res$fsthALL
                                                FSTN2           <- res$fstnALL
                                                GST             <- res$GstAll
                                                SNN             <- snn(haplotypes,Populations)
						      
                                                #res             <- fstcalc(haplotypes,Populations,data=NULL)
                                                #FSTN2           <- res$FSTALL
                                                #print(hapw)
                                                #print(Pi)
                                                #print(FSTH)
                                                #print(GST)
                                                #print(FSTN2)

                                              }
            
           			
				
                              } else { ### <- keine segregating sites
				
				                                nix       <- rep(NaN,npop)
                                        if(neutrality){
				                               	SEG       <- nix
                                        TD        <- nix
					                              R2        <- nix
					                              FuLi_F    <- nix
					                              FuLi_D    <- nix
                                        FS        <- nix
                                        thetaT    <- nix
                                        thetaS    <- nix
                                        Strobeck  <- nix
                                        HnFw      <- nix
                                        Ez        <- nix

                                        }

                                        if(linkage){
					                              WALLB      <- nix
					                              WALLQ      <- nix
					                              Zns        <- nix
					                              ZA         <- nix
					                              ZZ         <- nix
                                        }
                                        
                                        if(F_ST){
                                        hapw            <- nix
                                        Pi              <- nix
                                        FSTH            <- NaN
                                        FSTN2           <- NaN
                                        GST             <- NaN   
					SNN             <- NaN
 
                                        }
	
			                               		positions <- c(0)
					
		        } # end keine segsites
						
				if (zz==1) {		
					# if the first entry is made to the matrix add speaking column names and
					# csStats <- rbind(c(pi, tH, tD, flD, flF, fwH, nsegsites))
           			csStats <- vector("list",npop)
          for(xx in 1:npop){

                   if(neutrality){csStats[[xx]]  <- rbind(c(TD[xx],SEG[xx],R2[xx],FuLi_F[xx],FuLi_D[xx],FS[xx],Strobeck[xx],HnFw[xx],Ez[xx],thetaT[xx],thetaS[xx]))}
                   if(linkage){csStats[[xx]]     <- rbind(c(WALLB[xx],WALLQ[xx],ZA[xx],ZZ[xx],Zns[xx]))}
                   if(F_ST){hapwX <- hapw[xx];PiX<-Pi[xx];csStats[[xx]]<- rbind(c(hapwX,PiX,FSTH,FSTN2,GST,SNN))} 
                   colnames(csStats[[xx]])      <- c(testNames)
        	        # first = FALSE
	  } # End \FCber alle Populationen


			  		  
       } else {
#				  # add all computed test statistics for each sample to the stats matrix
#					# csStats <- rbind(csStats, c(pi, tH, tD, flD, flF, fwH, nsegsites))				
#	
    for(xx in 1:npop){
              # print(csStats[[xx]])
               if(neutrality){csStats[[xx]]         <- rbind(csStats[[xx]],c(TD[xx],SEG[xx],R2[xx],FuLi_F[xx],FuLi_D[xx],FS[xx],Strobeck[xx],HnFw[xx],Ez[xx],thetaT[xx],thetaS[xx]))}
               if(linkage)   {csStats[[xx]]        <- rbind(csStats[[xx]],c(WALLB[xx],WALLQ[xx],ZA[xx],ZZ[xx],Zns[xx]))}
               if(F_ST)      {hapwX<-hapw[xx];PiX <- Pi[xx];csStats[[xx]]   <- rbind(csStats[[xx]],c(hapwX,PiX,FSTH,FSTN2,GST,SNN))}
          
      }
      
    } #if
		
    
} # end of each simulation biallelic matrix	


		
 		

		# Hier sind die samples abgearbeitet !!!!!! -------------------------------------------------------------- ENDE der SAMPLES f\FCr ein Locus
		csStats           <- as.matrix(csStats)
		rownames(csStats) <- paste("pop",1:npop)
		colnames(csStats) <- "cs.stats"	  
		
    # save all calculated data in slots of the locstats class
		lData <- new("loc.stats", 
			n.sam 		    = nsam[cloci], 
			n.iter 		    = niter,
			theta		    = theta[cloci],
			seeds		    = seeds,
			positions	    = as.matrix(positions),
			trees		    = trees,								
			#haplotypes = haplotypes,
			stats 		    = csStats
		)
		
	
    if (!is.na(obsVal[1][1])) {
			
      obsVal_my <- matrix(,npop,dim(csStats[[1]])[2])
			rownames(obsVal_my) <- paste("pop",1:npop)
			colnames(obsVal_my) <- testNames
			
      ## Konvert obsVal --------------------------

      for(xx in 1:npop){
       
       obsVal_my[xx,] <- obsVal[[xx]][cloci,]
			}
			# for one gene !!!
      #lData@obsVal = obsVal[cloci,]
		  lData@obs.val = obsVal_my
    }
		
		lData <- calc_probabilities(cloci, niter, csStats, obsVal,lData,npop,testNames,numTests) # ------------------- calc_probabilities
		
    # accumulate all locstats objects in one list
		if (cloci == 1) {
			locusData <- c(lData)
		} else {
			locusData <- c(locusData, lData)
		}	
  } # for 1 to nloci  <------------------------------------------------------------- ENDE EINES GEN (LOCI)
  
  

if(is.list(obsVal)){return(locusData)} ### return only when observed Data is present
	
# Prepare the summary class over all loci	
############################################
# Init	
popnames           <- paste("pop",1:npop)
init               <- vector("list",npop)
average            <- as.matrix(init)
variance           <- as.matrix(init)
rownames(average)  <- popnames
rownames(variance) <- popnames
colnames(variance) <- "variance"
colnames(average)  <- "average"
#-----------------------------

for(xx in 1:npop){
	
	# calc average and variance for all loci
	average[[xx]]           <- matrix(ncol=numTests, nrow=niter)
	colnames(average[[xx]]) <- testNames
	
	variance[[xx]]           <- matrix(ncol=numTests, nrow=niter)
	colnames(variance[[xx]]) <- testNames
		
	for (i in 1:niter) {
		for (j in 1:numTests) {
			tmp <- matrix(ncol=1)
			
			for (k in 1:nloci) {	
				tmp <- rbind(tmp, locusData[[k]]@stats[[xx]][i,][j])
			}

			average[[xx]][i,][j]  <- colMeans(tmp, na.rm=TRUE)
			variance[[xx]][i,][j] <- var(tmp, na.rm=TRUE)
		}
	}
	
} # End of for pops	
  # return(as.matrix(average))
	
  # aggregate all calculated stats in an csstats class
	# data of all each single locus is accumulated in a list of class locstats
	lociData <- new("cs.stats", 
		n.loci		= nloci,
		n.iter		= niter,
		n.pop		= npop,
		
		average		= average,
		variance	= variance,
		
		locus		= locusData
	)

	
#	return(lociData)
	
	init      <- as.matrix(vector("list",npop))
	probLess  <- init
	probEqual <- init
	validIter <- init
	rownames(probLess)   <- popnames
	rownames(probEqual)  <- popnames
	colnames(probLess)   <- "prob.less"
	colnames(probEqual)  <- "prob.equal"
	rownames(validIter)  <- popnames
	colnames(validIter)  <- "valid.iter"

	
	# calc probabilities of all loci if observed value has been given to test against
	if ( !is.na(obsVal[1][1]) ) {
	
  for(xx in 1:npop){
	
		probLess[[xx]]           <- matrix(nrow=nloci, ncol=numTests)
		colnames(probLess[[xx]]) <- testNames
		
		probEqual[[xx]] 	  <- matrix(nrow=nloci, ncol=numTests)
		colnames(probEqual[[xx]]) <- testNames
		
		validIter[[xx]] 	<- matrix(nrow=nloci, ncol=numTests)
		colnames(validIter[[xx]]) <- testNames
		
		
		# iterate through all locStats objects and create a matrix with all probabilities
		# this will provide a single matrix with average stats for all loci and is made accessible
		# through the csStats class
		
		# return(locusData)
		for (elem in 1:length(locusData)) {
		
			probLess[[xx]][elem,]  = locusData[[elem]]@loc.prob.less[[xx]]
			probEqual[[xx]][elem,] = rbind(locusData[[elem]]@loc.prob.equal[[xx]])
			validIter[[xx]][elem,] = rbind(locusData[[elem]]@loc.valid.iter[[xx]])
		
		}
		
		#return(probLess)
		
  	        # calculate average and variance values for probabilties
		averageL        <- colMeans(probLess[[xx]], na.rm=TRUE)
		varianceL       <- apply(probLess[[xx]], 2, var, na.rm=TRUE)
		probLess[[xx]]  <- rbind(probLess[[xx]], averageL, varianceL)
		
		# calculate average and variance values for probabilties
                averageE        <- colMeans(probEqual[[xx]], na.rm=TRUE)
		varianceE       <- apply(probEqual[[xx]], 2, var, na.rm=TRUE)
		probEqual[[xx]] <- rbind(probEqual[[xx]], averageE, varianceE)
		
		# calculate average and variance values for the number of valid iterations
                averageV        <- colMeans(validIter[[xx]], na.rm=TRUE)
		varianceV       <- apply(validIter[[xx]], 2, var, na.rm=TRUE)
		validIter[[xx]] <- rbind(validIter[[xx]], averageV, varianceV)
		
		# assign names to rows in format 'locusX', where X is 1,2,3....
		# for each considered loci
                rnames <- paste("locus", 1:nloci)
                rnames <- c(rnames,"average","variance")
                rownames(probEqual[[xx]]) <- rnames
                rownames(probLess[[xx]])  <- rnames
                rownames(validIter[[xx]]) <- rnames
		
	}# of for npops	
	  	 
		# assign all calculated matrices to the created csStats objects
		lociData@prob.less 	= probLess	
		lociData@prob.equal	= probEqual
		lociData@valid.iter	= validIter
 }# End if
		
	# add obsVal to csStats class if present
	if ( !is.na(obsVal[1][1] )) {
		#colnames(obsVal) <- testNames
		lociData@obs.val = obsVal
	}
		
	
		
	# this simply creates a line break before printing the success notice 
	#print("", quote=FALSE);
	
	#print("Completed successfully!")
	
	# object of class csStats is returned
	return(lociData)
	}

#------------------------------------------------------------------------------#
#								getDemographyParam							   #
#------------------------------------------------------------------------------#
## This function verifies the validity of the 'demography' vector and returns a 
## string to invoke ms with changing demographic structure of the population
getDemographyParam <- function(d) {
	
	if ( length(d) == 3 ) {
		
		if (d[1] == 1) {
			
			# set all growth rates to alpha at time t.
			return( paste("-eG", d[2], d[3]) )
		} else {
			if (d[1] == 2) {
				
				# set all subpop\D5s to size x * N_0 and growth rates to zero.
				return( paste("-eN", d[2], d[3]) )
			} else {
				if (d[1] == 3) {
					
					# set all elements of migration matrix to x/(npop-1)
					return( paste("-eM", d[2], d[3]) )
				} else {
					stop("Invalid vector length for 'demography' for this type")
				}
			}
		}
	} else {
		
		if ( length(d) == 4 ) {
			
			if (d[1] == 4) {
				
				# set growth rate of subpop i to alpha at time z
				return( paste("-eg", d[2], d[3], d[4]) )
			} else {
				
				if (d[1] == 5) {
					
					# set subpop i size to x * N_0 at time t and growth rate to zero
					return( paste("-en", d[2], d[3], d[4]) )
				} else {
					
					if (d[1] == 6) {
						
						# subpopulation slit
						return( paste("-es", d[2], d[3], d[4]) )
					} else {
						
						if (d[1] == 7) {
							
							# move between subpopulations
							return( paste("-ej", d[2], d[3], d[4]) )
						} else {
							stop("Invalid vector length for 'demography' for this type.")
						}
						
					}
				}
			}
			
		} else {
			stop("Invalid length of vector 'demography'.")
		}
	}
	
}

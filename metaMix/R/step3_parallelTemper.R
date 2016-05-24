################################################ Perform Parallel Tempering MCMC ###############################################
#' @name parallel.temper
#' @title Parallel Tempering MCMC
NULL
#' @rdname parallel.temper
#' @title parallel.temper
#' @description Performs Parallel Tempering MCMC to explore the species state space. Two types of moves are implemented: a mutation step (within chain) and an exchange step (between neighboring chains).  If working with BLASTn data, use parallel.temper.nucl().
#' @param readSupport The number of reads the user requires in order to believe in the presence of the species. It is used to compute the penalty factor. The default value is 10. We compute the logarithmic penalty value as the log-likelihood difference between two models: one where all N reads belong to the "unknown"  category and one where r reads have a perfect match to some unspecified species and the remaining reads belong to the "unknown"  category.  
#' @param noChains The number of parallel chains to run. The default value is 12.
#' @param seed Optional argument that sets the random seed (default is 1) to make results reproducible.
#' @param step2 list. The output from reduce.space(). Alternatively, it can be a character string containing the path name of the ".RData" file  where step2 list was saved.
#' @return step3: A list with two elements. The first one (result) is a list that records MCMC information from each parallel chain.  The second one (duration) records how much time the MCMC exploration took.
#' @seealso  \code{\link{parallel.temper.nucl}} This function should be used when working with BLASTn data.
#' @keywords parallel.temper
#' @export parallel.temper
#' @import Rmpi Matrix
#' @importFrom gtools rdirichlet
#' @examples
#' ## See vignette for more details
#'
#' \dontrun{
#' # Either load the object created by previous step (i.e from function reduce.space() )
#' data(step2)   ## example output of reduce.space
#' step3<-parallel.temper(step2=step2)
#'
#' # or alternatively point to the location of the step2.RData object
#' step3 <- parallel.temper(step2="/pathtoFile/step2.RData")
#' }
######################################################################################################################
parallel.temper = function(step2, readSupport=10, noChains=12, seed=1){

  if (is.character(step2)) {
    load(step2)
  }
  

  should.be.in.the.list <- c("pij.sparse.mat", "ordered.species", "read.weights", "outDir", "gen.prob.unknown")

  if (sum (!( should.be.in.the.list %in% names(step2)) ) > 0) {
    message('Missing the following arguments')
    print(names(step2)[!(should.be.in.the.list %in% names(step2))] )
    stop()
  }  else {
    
    parallel.temper.wrapped<-function(readSupport.internal=readSupport, noChains.internal=noChains, pij.sparse.mat=step2$pij.sparse.mat, read.weights=step2$read.weights, ordered.species=step2$ordered.species, gen.prob.unknown=step2$gen.prob.unknown, outDir=step2$outDir, seed.internal=seed){

      set.seed(seed.internal);
                                        #print(warnings())
      StartTime<-Sys.time()
      
      sieve <- function(n) {
        n <- as.integer(n)
        if(n > 1e6) stop("n too large")
        primes <- rep(TRUE, n)
        primes[1] <- FALSE
        last.prime <- 2L
        for(i in last.prime:floor(sqrt(n)))
          {
            primes[seq.int(2L*last.prime, n, last.prime)] <- FALSE
            last.prime <- last.prime + min(which(primes[(last.prime+1):n]))
          }
        which(primes)
      }

      list.integers <- sieve(1000)
      
      node.ids <- list.integers[ 1:noChains.internal ]
      mpi.spawn.Rslaves(nslaves = noChains.internal)  #number of slaves to spawn, should be equal to individual chains

# In case R exits unexpectedly, have it automatically clean up
# resources taken up by Rmpi (slaves, memory, etc...)
      .Last <- function(){
        if (is.loaded("mpi_initialize")){
          if (mpi.comm.size(1) > 0){
            print("Please use mpi.close.Rslaves() to close slaves.")
            mpi.close.Rslaves()
          }
          print("Please use mpi.quit() to quit R")
          .Call("mpi_finalize", PACKAGE='metaMix')
        }
      }



      pij.sparse.mat<-cBind(pij.sparse.mat, "unknown"=gen.prob.unknown)

#rm(step2)
      gc()

      allSpecies<-ordered.species[,c("taxonID", "samplingWeight")]
      lenSp<-nrow(allSpecies)

      lpenalty<-computePenalty(readSupport=readSupport.internal, readWeights=read.weights, pUnknown=gen.prob.unknown)  ### penalty for accepting a new species ~ Use pij for the r readSupport perfect reads (=1/median(protein length)). Likelihood jump for 10 reads moving from unknown bin (pij=gen.prob.unknown, default 1e-10) to species (1/1500).

### EMiter
      EMiter<-10

### PT parameters
      exchangeInterval<-1                             ###leave chains run in parallel for that many iterations
      ExternIter<-5*(nrow(ordered.species))                              ### make chains communicate 1 times
      TotalIter<-exchangeInterval * ExternIter

##Tempering Vector --Power Decay
      temper<-vector()
      K<-0.001
      a<-3/2
      for (n in 2:noChains.internal){
        temper[1]<-1
        temper[n]<-(temper[n-1]-K)^a
      }
      for (i in 1:noChains.internal) {names(temper)[i]<- paste("slave",i, sep="")}  ###names

### flag for adding /removing species
      stepAdd<-vector(mode = "logical")

#### broadcast necessary objects/functions to slaves
      mpi.bcast.Robj2slave(exchangeInterval)
      mpi.bcast.Robj2slave(ExternIter)
      mpi.bcast.Robj2slave(TotalIter)
      mpi.bcast.Robj2slave(list.integers)
      mpi.bcast.Robj2slave(node.ids)
      mpi.bcast.Robj2slave(ordered.species)
      mpi.bcast.Robj2slave(read.weights)
      mpi.bcast.Robj2slave(pij.sparse.mat)
      mpi.bcast.Robj2slave(lpenalty)
      mpi.bcast.Robj2slave(EM)
      mpi.bcast.Robj2slave(allSpecies)
      mpi.bcast.Robj2slave(lenSp)
      mpi.bcast.Robj2slave(EMiter)
      mpi.bcast.Robj2slave(noChains.internal)
      mpi.bcast.Robj2slave(temper)
      mpi.bcast.Robj2slave(gen.prob.unknown)
##########################---------------------------------------SINGLE CHAIN FUNCTION -----------------------------------------------------------------------------------------------#######################
      singleChain <- function(TotalIter, exchangeInterval){

        if (file.exists('~/.Rprofile')) source('~/.Rprofile')
        print(.libPaths())
  
       
        ind <- mpi.comm.rank()   # Each slave gets its own copy of ind and chain based on mpi process rank
        print(ind)

        estimNew<-matrix(0, nrow=TotalIter, ncol=1)   ### Create matrix that will hold the estimator of the log-likelihood
        estimNew[1,]<- sum(read.weights[,"weight"])*log(gen.prob.unknown)*temper[ind]  
        print(estimNew[1,])
  
        presentSpecies<-"unknown"   ## create object that receives the species deemed as present, from previous (swapInterv*j) iteration.
        abundUsedSp<-c("unknown"=1)
  
  ### Create 2 lists. One that holds species names and one with their abundances. Each list has 2 elements. Element1: present species  and Element2: tentative species being explored in current iteration.
        speciestoUse<-list("presentSp"=presentSpecies, "tentativeSp"=NULL)
        abundUsedSpecies<-list("presentSp"= abundUsedSp, "tentativeSp"=NULL)
  
        message("\nThis is the temperature for this slave")
        print(temper[ind])
  
### create matrix to hold info on which species was added or removed, which species were present and logL , through iterations
        record<-matrix(0, ncol=(5+lenSp), nrow=(TotalIter))
        record<-as.data.frame(record)
        colnames(record)<-c("Iter", "Move", "Candidate Species", allSpecies[,1],  "unknown", "logL")
        record[,1]=1:(TotalIter)
        record[1, presentSpecies]<-1
        record[1,(5+lenSp)]<-estimNew[1,]             ###last colum is logL
  
        swaps.attempted<-0
        swaps.accepted<-0
        oddFlag<-0

################ ----------------------------------- Begin MCMC (within single chain)  ------------------------------------------------------------------------------------------############  
        for (i in  2:TotalIter)  {    
          cat('\n',i,'\n')                                 ### temporary - debugging purposes

#### Create 3 functions to use repeatedly in adding/removing species steps.
######################### 1a. For add step:  species  I sample my candidate species from.
          potentialSpecies<-function() {                    
            toChooseFrom <-  allSpecies[,1][!(allSpecies[,1] %in% speciestoUse[[1]])]   ### Species remaining after omitting "present species"
            toChooseFrom<-as.character(toChooseFrom)
            potentialSp<- allSpecies[allSpecies$taxonID %in% toChooseFrom,]
            resultPS<-list("potentialSp"=potentialSp, "toChooseFrom"=toChooseFrom)
            return(resultPS)
          }

    ######################## 1b. For add step, sample candidate organism to add and create object for tentative species
          moveAdd<-function() {
            randSp <- sample(as.character(potentialSp[,1]),1, prob=potentialSp[,2])   ### choose random species, weights need not sum to one. Here weights are based on intial read counts.         
            cat('Add species', randSp, '\n')                     ### temporary - debugging purposes
            speciestoUse[[2]]<-c(randSp, speciestoUse[[1]])      ### tentative present species, to use in Gibbs below.
            record[i, 2]<- "Add"                     ### record info on species we are adding
            record[i,3]<- randSp

            result<-list("record"=record, "speciestoUse"=speciestoUse)      
            return(result)
          }

    ######################## 2. For remove step, sample candidate organism to remove and create object for tentative species
          moveRemove<-function() {
            toremoveFrom <- speciestoUse[[1]][!(speciestoUse[[1]]%in%"unknown")]         ### species set from which to remove (i.e present ones bar unknown)
            toremoveFrom<-as.character(toremoveFrom)          
            removeProb<-(1/abundUsedSpecies[[1]][toremoveFrom])/sum(1/abundUsedSpecies[[1]][toremoveFrom])       ###### sampling weight inversely proportional to assigned read counts.
     ###Flatten removal probabilities
            percentiles<-quantile(removeProb,  probs=c(0.2, 0.8))
            removeProb[which(removeProb  >= percentiles["80%"])] <- percentiles["80%"]
            removeProb[which(removeProb  <= percentiles["20%"])] <- percentiles["20%"]

    ### sample candidate species to remove
            randSp<-sample(toremoveFrom, 1, prob=removeProb )   ### weights need not sum to one
            cat('Remove species', randSp, '\n')                 ### temporary - debugging purposes
            speciestoUse[[2]] <- speciestoUse[[1]][!(speciestoUse[[1]] %in% randSp)]    ### tentative present species, to use in Gibbs below.
            record[i, 2]<- "Remove"                     ### record info on species we are adding
            record[i,3]<- randSp
            result<-list("record"=record, "speciestoUse"=speciestoUse)      
            return(result)
          }

          addStep <-0.5      ####proability of doing add step
    
              ######################################################## Add species ########################################      
          if (runif(1)<= addStep) {                             
            resultPS<-potentialSpecies()
            potentialSp<-as.data.frame(resultPS$potentialSp)
            toChooseFrom<-as.character(resultPS$toChooseFrom)
            
            if (length(toChooseFrom) != 0L){  
              result<-moveAdd()
              speciestoUse<-result$speciestoUse
              record<-result$record
                
            } else {                                 ################## ### If pool of species to add empty, go to remove step  
              message("No more species to choose from, all are kept as present.")
              result<-moveRemove()
              speciestoUse<-result$speciestoUse
              record<-result$record        
              addStep<-0             ### so do remove step with prob=1
            }  

          }  ###close if runif(1)<=addStep

      
          else {    ################################################### Remove species #################################################
            toremoveFrom <- speciestoUse[[1]][!(speciestoUse[[1]]%in%"unknown")]         ### species set from which to remove (i.e present ones bar unknown)
            toremoveFrom<-as.character(toremoveFrom)

         #### First check that more than 1 species are present, so don't risk of remaining only with X bin and wasting iteration.
            if (length(toremoveFrom) > 1) {
              result<-moveRemove()
              speciestoUse<-result$speciestoUse
              record<-result$record        
        
            }  else if  (length(speciestoUse[[1]])==1) {                         #### speciestoUse[[1]]==1 when only unknown bin is there / We can only add
              message("We cannot remove, we only have unknown bin")
              potentialSp<-allSpecies
              result<-moveAdd()
              speciestoUse<-result$speciestoUse
              record<-result$record        
              addStep <- 1        ### so do add step with prob=1
        
            } else {                                                ### we are adding species instead. Removing would leave us just with X-bin again.
              message("\nCannot remove a species, so instead we add\n")   ### temporary
              resultPS<-potentialSpecies()
              potentialSp<-as.data.frame(resultPS$potentialSp)
              toChooseFrom<-as.character(resultPS$toChooseFrom)          
              result<-moveAdd()
              speciestoUse<-result$speciestoUse
              record<-result$record        
              addStep <- 1        
            }        
          } ###close else (remove species)


          noSpecies<-length(speciestoUse[[2]])
          print(noSpecies)
          hyperP<-rep(1, noSpecies)
          startW<-rdirichlet(1, hyperP)
          
          output100Tent<-EM(pij=pij.sparse.mat, iter=EMiter, species=speciestoUse[[2]], abund=startW, readWeights = read.weights)  ### EM function
          
          estimator <- output100Tent$logL[EMiter,2] + (noSpecies * lpenalty) #### penalise likelihood                                  
          
          message("EM took: ", output100Tent$RunningTime)       
          mean1<-as.numeric(output100Tent$abundances[EMiter,2:(noSpecies+1)])
          names(mean1)<-names(output100Tent$abundances[EMiter,2:(noSpecies+1)])                        
          abundUsedSpecies[[2]]<-mean1
          
          estimNew[i,]<-estimator*temper[ind]                          ###tempered likelihood
          message('\n', estimNew[i,] , ' ',  estimNew[i-1,])           ### temporary print - which values am I comparing?
    
    ####### flag of adding/removing
          if (record[i,2]=="Add") {
            stepAdd<-TRUE
          } else {stepAdd<-FALSE}

          removeStep <- 1 - addStep    ##remove step


          if (stepAdd) {

            candidateAdd<-allSpecies[allSpecies[,1]==(as.character(record[i,3])),2]
            removeProba<-(1/abundUsedSpecies[[2]])/sum(1/abundUsedSpecies[[2]])       ###### sampling weight inversely proportional to size. Tentative Species

###Flattening the sampling probabilities
            percentiles<-quantile(removeProba,  probs=c(0.2, 0.8))
            removeProba[which(removeProba  >= percentiles["80%"])] <- percentiles["80%"]
            removeProba[which(removeProba  <= percentiles["20%"])] <- percentiles["20%"]
            candidateRemove<-removeProba[as.character(record[i,3])]
            
            print(candidateRemove) 
            print(candidateAdd)
            if (runif(1) < min( 1, exp(estimNew[i] - estimNew[i-1] + log(removeStep) - log(addStep) + log(candidateRemove) - log(candidateAdd) )))  {     ### accept "add species"  with prob min{1, P(D|i)*P(removeSpecies)*P(remove Specific Species)/P(D|i-1)*P(addSpecies)*P(add specific species) 

        ######## Accept move #######
              estimNew[i,]<-estimNew[i,]                                  ### if move is accepted, record new logL
              speciestoUse[[1]]<-speciestoUse[[2]]                        ### if move is accepted, tentative species becomes present species.
              abundUsedSpecies[[1]]<-abundUsedSpecies[[2]]                        ###       -->>-->>--   , abundances of new set of species are kept  
              cat('Present species become:', speciestoUse[[1]], '\n')

            }  else {                                                    ######### Reject Move ########
              estimNew[i,]<-estimNew[i-1,]                               ### if move is rejected, record previous logL
              speciestoUse[[1]]<-speciestoUse[[1]]                       ### if move is rejected, keep present species as it is.
              abundUsedSpecies[[1]]<-abundUsedSpecies[[1]]
              cat('Present species remain:', speciestoUse[[1]], '\n')    ### temporary
            }

          }    ### close  if(stepAdd==TRUE)

    
          else {   ########################################  type of move was to remove species   ##############
            removeProba<-(1/abundUsedSpecies[[1]])/sum(1/abundUsedSpecies[[1]])   ###### sampling weight inversely proportional to size. Do this for present Species
            candidateRemove<-removeProba[as.character(record[i,3])]
            candidateAdd<-allSpecies[allSpecies[,1]==(as.character(record[i,3])),2]
            
      
            if (runif(1) < min( 1, exp(estimNew[i] - estimNew[i-1] + log(addStep) - log(removeStep) + log(candidateAdd) - log(candidateRemove) ))) {############# Accept "remove species"        
              estimNew[i,]<-estimNew[i,]
              speciestoUse[[1]]<-speciestoUse[[2]]                       ### tentative species becomes present species.
              abundUsedSpecies[[1]]<-abundUsedSpecies[[2]]
              cat('Present species become:', speciestoUse[[1]], '\n')    
        
            }  else {                                                    ######## Reject Move ###########
              estimNew[i,]<-estimNew[i-1,]
              speciestoUse[[1]]<-speciestoUse[[1]]
              abundUsedSpecies[[1]]<-abundUsedSpecies[[1]]
              cat('Present species remain:', speciestoUse[[1]], '\n')          
            }

          }   ###close else (i.e (stepAdd==FALSE))
    
#################################edw tha kanw attempt to direct swap between slaves, without going through master
       
          if( i%%exchangeInterval == 0 ){   ### every nth (2nd for now) iteration
            message("\n\nTime to attempt an exchange")      
            oddFlag<-oddFlag+1
            swap<-0
            estim.current<-estimNew[i,]/temper[ind]    ######### need untempered logL
            
      
########################################### CREATE prime tags for object to send around
            allowedLength<-175
            Nsubobjects<-round(length(abundUsedSpecies[[1]])/allowedLength)+1
            object.ids <- list.integers[ (noChains.internal+1):(noChains.internal + 4 + Nsubobjects) ]    ### 4 objcts for logL, swap message, untempered , species PLUS as many as necessary for abundances            
            
            if (ind%%2 == oddFlag%%2) {  ###when oddFlag zero , the following code concerns even-numbered slaves. For oddFlag 1, it concerns odd-numbered slaves.                 
              ind.partner<-ind+1

              if (0<ind.partner && ind.partner<(noChains.internal+1)){           
                estim.partner<-mpi.recv.Robj(ind.partner, tag=object.ids[1]*node.ids[ind.partner])  #### receive the logdensity of above  partner
                message("I received the untempered: ", estim.partner)
                swaps.attempted<-swaps.attempted+1          
                lalpha<-(estim.partner - estim.current)*(temper[ind] - temper[ind.partner] )
                message("This is the acceptance probability: ", min(1,exp(lalpha)))
          
                if (runif(1)< min(1, exp(lalpha))) {    ############# exp((chain2 - chain1)*(T1 - T2))
                  message("I exchanged the values")               
                  swap<-1         
                }   ## end of if runif(1)<lalpha M-H step

                else {message("I didn't exchange the values")}
                
                mpi.send.Robj(obj=swap,dest=ind.partner,tag=object.ids[2]*node.ids[ind])
                message("I send message swap: ", swap)
              }     ### end of   if (0<ind.partner && ind.partner<(noChains.internal+1)){          

              if(swap==1){
                swaps.accepted<-swaps.accepted+1

                mpi.send.Robj(obj=estim.current,dest=ind.partner,tag=object.ids[3]*node.ids[ind])
                message("I just sent the untempered: ", estim.current)
                species.swap<-speciestoUse[[1]]
                mpi.send.Robj(species.swap,dest=ind.partner,tag=object.ids[4]*node.ids[ind])
                speciestoUse[[1]]<-mpi.recv.Robj(ind.partner,tag=object.ids[4]*node.ids[ind.partner])
          
                abund.swap<-abundUsedSpecies[[1]]
                mpi.send.Robj(abund.swap,dest=ind.partner,tag=object.ids[5]*node.ids[ind])
                abundUsedSpecies[[1]]<-mpi.recv.Robj(ind.partner,tag=object.ids[5]*node.ids[ind.partner])
          

                message("I received the following abundances from: ", ind.partner, " here ", ind)
                message("The new logdensity will be the one of the neighbor but tempered with my temperature: ", estim.partner*temper[ind])
                estimNew[i,]<-estim.partner * temper[ind]
              }     ### end of if(swap==1)         
        
            } else {  ##### ###when oddFlag zero , the following code concerns odd-numbered slaves. For oddFlag 1, it concerns even-numbered slaves. I.e say what partners should do
              ind.partner<-ind-1;           
              if(0<ind.partner && ind.partner<(noChains.internal+1)){
                mpi.send.Robj(obj=estim.current,dest=ind.partner,tag=object.ids[1]*node.ids[ind])
                message("I just sent the untempered: ", estim.current)
                swap<-mpi.recv.Robj(ind.partner,tag=object.ids[2]*node.ids[ind.partner])
                message("I received the swap message: ", swap)
                swaps.attempted<-swaps.attempted+1
        }  ####end of  if(0<ind.partner && ind.partner<(noChains.internal+1)){
          
                  
              if(swap==1){
                swaps.accepted<-swaps.accepted+1
                estim.partner<-mpi.recv.Robj(ind.partner, tag=object.ids[3]*node.ids[ind.partner] )  #### receive the logdensity of above  partner
                message("I received the untempered: ", estim.partner)

                species.swap<-speciestoUse[[1]]        
                
                speciestoUse[[1]]<-mpi.recv.Robj(ind.partner,tag=object.ids[4]*node.ids[ind.partner])           
                mpi.send.Robj(species.swap,dest=ind.partner,tag=object.ids[4]*node.ids[ind])
          
                abund.swap<-abundUsedSpecies[[1]]
                abundUsedSpecies[[1]]<-mpi.recv.Robj(ind.partner,tag=object.ids[5]*node.ids[ind.partner])
                mpi.send.Robj(abund.swap,dest=ind.partner,tag=object.ids[5]*node.ids[ind])
                    
                message("I received the following species from: ", ind.partner, " here ", ind)
                message("The new logdensity will be the one of the neighbor but tempered with my temperature: ", estim.partner*temper[ind])
                estimNew[i,]<-estim.partner * temper[ind]          
              }
            }
          }

############ Use indicator variable for species presence per iteration
          record[i, speciestoUse[[1]]]<-ifelse(speciestoUse[[1]] %in% colnames(record), 1, 0)
          record[i,(5+lenSp)]<-estimNew[i,]             ###last colum is logL
    
        }    ###end of for loop for internal iterations
  
        resultSC<-list("estimNew"=estimNew, "record"=record, "usedSp"=speciestoUse[[1]], "abundUsedSp" = abundUsedSpecies[[1]], "swaps.attempt"=swaps.attempted, "swaps.accept"=swaps.accepted, "readSupport"=readSupport.internal, "lpenalty"=lpenalty)
        return(resultSC)  
   
      }    ##end of singleChain function


##########################################----------------------- MAIN ---------------------------
###send function to slaves
      mpi.bcast.Robj2slave(singleChain)

### Now call the single chain function
      message('Start parallel tempering')
      result<-mpi.remote.exec(singleChain(TotalIter, exchangeInterval))
      EndTime<-Sys.time()
      duration<-EndTime-StartTime
      message('PT finished in ', duration, ". ", nrow(ordered.species), " species were explored in ", ExternIter, "x", noChains.internal, " iterations.")
      
      step3<-list("result"=result, "duration"=duration)
      
      if (!is.null(outDir)) {
        step3.name <- paste(outDir, "/step3.RData", sep = "")
        save(step3, file=step3.name)
        rm(list= ls()[!ls() %in% c("step3")])
        gc()    
      } else {
        rm(list= ls()[!ls() %in% c("step3")])
        gc()
      }

  
      mpi.close.Rslaves(dellog=FALSE)
   #   mpi.quit()

      return(step3)
    }
    parallel.temper.wrapped()
  }

}

#' @rdname parallel.temper
#' @title parallel.temper.explicit
#' @description  parallel.temper.explicit is the same function as parallel.temper but with a more involved syntax. 
#' @param pij.sparse.mat sparse matrix of generative probabilities, see value of ?reduce.space. 
#' @param read.weights   see ?reduce.space.
#' @param ordered.species see ?reduce.space.
#' @param gen.prob.unknown  see ?reduce.space.
#' @param outDir see ?reduce.space.
#' @keywords parallel.temper.explicit
#' @export parallel.temper.explicit
#' @import Rmpi Matrix
#' @importFrom gtools rdirichlet

parallel.temper.explicit<-function(readSupport=10, noChains=12, pij.sparse.mat, read.weights, ordered.species, gen.prob.unknown, outDir, seed = 1){

  set.seed(seed);
#print(warnings())
StartTime<-Sys.time()

sieve <- function(n) {
  n <- as.integer(n)
  if(n > 1e6) stop("n too large")
  primes <- rep(TRUE, n)
  primes[1] <- FALSE
  last.prime <- 2L
  for(i in last.prime:floor(sqrt(n)))
    {
      primes[seq.int(2L*last.prime, n, last.prime)] <- FALSE
      last.prime <- last.prime + min(which(primes[(last.prime+1):n]))
    }
          which(primes)
}

list.integers <- sieve(1000)

node.ids <- list.integers[ 1:noChains ]
mpi.spawn.Rslaves(nslaves = noChains)  #number of slaves to spawn, should be equal to individual chains

# In case R exits unexpectedly, have it automatically clean up
# resources taken up by Rmpi (slaves, memory, etc...)
.Last <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0){
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize", PACKAGE='metaMix')
  }
}



pij.sparse.mat<-cBind(pij.sparse.mat, "unknown"=gen.prob.unknown)

#rm(step2)
gc()

allSpecies<-ordered.species[,c("taxonID", "samplingWeight")]
lenSp<-nrow(allSpecies)

lpenalty<-computePenalty(readSupport=readSupport, readWeights=read.weights, pUnknown=gen.prob.unknown)  ### penalty for accepting a new species 

### EMiter
EMiter<-10

### PT parameters
exchangeInterval<-1                             ###leave chains run in parallel for that many iterations before attempting exchange
ExternIter<-5*(nrow(ordered.species))                              ### make chains communicate 1 times
TotalIter<-exchangeInterval * ExternIter

##Tempering Vector --Power Decay
temper<-vector()
K<-0.001
a<-3/2
for (n in 2:noChains){
  temper[1]<-1
  temper[n]<-(temper[n-1]-K)^a
}
for (i in 1:noChains) {names(temper)[i]<- paste("slave",i, sep="")}  ###names

### flag for adding /removing species
stepAdd<-vector(mode = "logical")

#### broadcast necessary objects/functions to slaves
mpi.bcast.Robj2slave(exchangeInterval)
mpi.bcast.Robj2slave(ExternIter)
mpi.bcast.Robj2slave(TotalIter)
mpi.bcast.Robj2slave(list.integers)
mpi.bcast.Robj2slave(node.ids)
mpi.bcast.Robj2slave(ordered.species)
mpi.bcast.Robj2slave(read.weights)
mpi.bcast.Robj2slave(pij.sparse.mat)
mpi.bcast.Robj2slave(lpenalty)
mpi.bcast.Robj2slave(EM)
mpi.bcast.Robj2slave(allSpecies)
mpi.bcast.Robj2slave(lenSp)
mpi.bcast.Robj2slave(EMiter)
mpi.bcast.Robj2slave(noChains)
mpi.bcast.Robj2slave(temper)
mpi.bcast.Robj2slave(gen.prob.unknown)
##########################---------------------------------------SINGLE CHAIN FUNCTION -----------------------------------------------------------------------------------------------#######################
singleChain <- function(TotalIter, exchangeInterval){

  if (file.exists('~/.Rprofile')) source('~/.Rprofile')
  print(.libPaths())
  
       
  ind <- mpi.comm.rank()   # Each slave gets its own copy of ind and chain based on mpi process rank
  print(ind)

  estimNew<-matrix(0, nrow=TotalIter, ncol=1)   ### Create matrix that will hold the estimator of the log-likelihood
  estimNew[1,]<- sum(read.weights[,"weight"])*log(gen.prob.unknown)*temper[ind]  
  print(estimNew[1,])
  
  presentSpecies<-"unknown"   ## create object that receives the species deemed as present, from previous (swapInterv*j) iteration.
  abundUsedSp<-c("unknown"=1)

  ### Create 2 lists. One that holds species names and one with their abundances. Each list has 2 elements. Element1: present species  and Element2: tentative species being explored in current iteration.
  speciestoUse<-list("presentSp"=presentSpecies, "tentativeSp"=NULL)
  abundUsedSpecies<-list("presentSp"= abundUsedSp, "tentativeSp"=NULL)
  
  message("\nThis is the temperature for this slave")
  print(temper[ind])
  
### create matrix to hold info on which species was added or removed, which species were present and logL , through iterations
  record<-matrix(0, ncol=(5+lenSp), nrow=(TotalIter))
  record<-as.data.frame(record)
  colnames(record)<-c("Iter", "Move", "Candidate Species", allSpecies[,1],  "unknown", "logL")
  record[,1]=1:(TotalIter)
  record[1, presentSpecies]<-1
  record[1,(5+lenSp)]<-estimNew[1,]             ###last colum is logL
  
  swaps.attempted<-0
  swaps.accepted<-0
  oddFlag<-0

################ ----------------------------------- Begin MCMC (within single chain)  ------------------------------------------------------------------------------------------############  
  for (i in  2:TotalIter)  {    
    cat('\n',i,'\n')                                 ### temporary - debugging purposes

#### Create 3 functions to use repeatedly in adding/removing species steps.
######################### 1a. For add step:  species  I sample my candidate species from.
    potentialSpecies<-function() {                    
      toChooseFrom <-  allSpecies[,1][!(allSpecies[,1] %in% speciestoUse[[1]])]   ### Species remaining after omitting "present species"
      toChooseFrom<-as.character(toChooseFrom)
      potentialSp<- allSpecies[allSpecies$taxonID %in% toChooseFrom,]
      resultPS<-list("potentialSp"=potentialSp, "toChooseFrom"=toChooseFrom)
      return(resultPS)
    }

    ######################## 1b. For add step, sample candidate organism to add and create object for tentative species
    moveAdd<-function() {
      randSp <- sample(as.character(potentialSp[,1]),1, prob=potentialSp[,2])   ### choose random species, weights need not sum to one. Here weights are based on intial read counts.         
      cat('Add species', randSp, '\n')                     ### temporary - debugging purposes
      speciestoUse[[2]]<-c(randSp, speciestoUse[[1]])      ### tentative present species, to use in Gibbs below.
      record[i, 2]<- "Add"                     ### record info on species we are adding
      record[i,3]<- randSp

      result<-list("record"=record, "speciestoUse"=speciestoUse)      
      return(result)
    }

    ######################## 2. For remove step, sample candidate organism to remove and create object for tentative species
    moveRemove<-function() {
      toremoveFrom <- speciestoUse[[1]][!(speciestoUse[[1]]%in%"unknown")]         ### species set from which to remove (i.e present ones bar unknown)
      toremoveFrom<-as.character(toremoveFrom)          
      removeProb<-(1/abundUsedSpecies[[1]][toremoveFrom])/sum(1/abundUsedSpecies[[1]][toremoveFrom])       ###### sampling weight inversely proportional to assigned read counts.
     ###Flatten removal probabilities
      percentiles<-quantile(removeProb,  probs=c(0.2, 0.8))
      removeProb[which(removeProb  >= percentiles["80%"])] <- percentiles["80%"]
      removeProb[which(removeProb  <= percentiles["20%"])] <- percentiles["20%"]

    ### sample candidate species to remove
      randSp<-sample(toremoveFrom, 1, prob=removeProb )   ### weights need not sum to one
      cat('Remove species', randSp, '\n')                 ### temporary - debugging purposes
      speciestoUse[[2]] <- speciestoUse[[1]][!(speciestoUse[[1]] %in% randSp)]    ### tentative present species, to use in Gibbs below.
      record[i, 2]<- "Remove"                     ### record info on species we are adding
      record[i,3]<- randSp
      result<-list("record"=record, "speciestoUse"=speciestoUse)      
      return(result)
    }

    addStep <-0.5      ####proability of doing add step
    
              ######################################################## Add species ########################################      
    if (runif(1)<= addStep) {                             
      resultPS<-potentialSpecies()
      potentialSp<-as.data.frame(resultPS$potentialSp)
      toChooseFrom<-as.character(resultPS$toChooseFrom)
      
      if (length(toChooseFrom) != 0L){  
        result<-moveAdd()
        speciestoUse<-result$speciestoUse
        record<-result$record
                
      } else {                                 ################## ### If pool of species to add empty, go to remove step  
        message("No more species to choose from, all are kept as present.")
        result<-moveRemove()
        speciestoUse<-result$speciestoUse
        record<-result$record        
        addStep<-0             ### so do remove step with prob=1
      }  

    }  ###close if runif(1)<=addStep

      
    else {    ################################################### Remove species #################################################
      toremoveFrom <- speciestoUse[[1]][!(speciestoUse[[1]]%in%"unknown")]         ### species set from which to remove (i.e present ones bar unknown)
      toremoveFrom<-as.character(toremoveFrom)

         #### First check that more than 1 species are present, so don't risk of remaining only with X bin and wasting iteration.
      if (length(toremoveFrom) > 1) {
        result<-moveRemove()
        speciestoUse<-result$speciestoUse
        record<-result$record        
        
      }  else if  (length(speciestoUse[[1]])==1) {                         #### speciestoUse[[1]]==1 when only unknown bin is there / We can only add
        message("We cannot remove, we only have unknown bin")
        potentialSp<-allSpecies
        result<-moveAdd()
        speciestoUse<-result$speciestoUse
        record<-result$record        
        addStep <- 1        ### so do add step with prob=1
        
      } else {                                                ### we are adding species instead. Removing would leave us just with X-bin again.
        message("\nCannot remove a species, so instead we add\n")   ### temporary
        resultPS<-potentialSpecies()
        potentialSp<-as.data.frame(resultPS$potentialSp)
        toChooseFrom<-as.character(resultPS$toChooseFrom)          
        result<-moveAdd()
        speciestoUse<-result$speciestoUse
        record<-result$record        
        addStep <- 1        
      }        
    } ###close else (remove species)


    noSpecies<-length(speciestoUse[[2]])
    print(noSpecies)
    hyperP<-rep(1, noSpecies)
    startW<-rdirichlet(1, hyperP)
   
    output100Tent<-EM(pij=pij.sparse.mat, iter=EMiter, species=speciestoUse[[2]], abund=startW, readWeights = read.weights)  ### EM function
        
    estimator <- output100Tent$logL[EMiter,2] + (noSpecies * lpenalty) #### penalise likelihood                                  

    message("EM took: ", output100Tent$RunningTime)       
    mean1<-as.numeric(output100Tent$abundances[EMiter,2:(noSpecies+1)])
    names(mean1)<-names(output100Tent$abundances[EMiter,2:(noSpecies+1)])                        
    abundUsedSpecies[[2]]<-mean1

    estimNew[i,]<-estimator*temper[ind]                          ###tempered likelihood
    message('\n', estimNew[i,] , ' ',  estimNew[i-1,])           ### temporary print - which values am I comparing?
    
    ####### flag of adding/removing
    if (record[i,2]=="Add") {
      stepAdd<-TRUE
    } else {stepAdd<-FALSE}

    removeStep <- 1 - addStep    ##remove step


    if (stepAdd) {

      candidateAdd<-allSpecies[allSpecies[,1]==(as.character(record[i,3])),2]
      removeProba<-(1/abundUsedSpecies[[2]])/sum(1/abundUsedSpecies[[2]])       ###### sampling weight inversely proportional to size. Tentative Species

###Flattening the sampling probabilities
      percentiles<-quantile(removeProba,  probs=c(0.2, 0.8))
      removeProba[which(removeProba  >= percentiles["80%"])] <- percentiles["80%"]
      removeProba[which(removeProba  <= percentiles["20%"])] <- percentiles["20%"]
      candidateRemove<-removeProba[as.character(record[i,3])]

      print(candidateRemove) 
      print(candidateAdd)
      if (runif(1) < min( 1, exp(estimNew[i] - estimNew[i-1] + log(removeStep) - log(addStep) + log(candidateRemove) - log(candidateAdd) )))  {     ### accept "add species"  with prob min{1, P(D|i)*P(removeSpecies)*P(remove Specific Species)/P(D|i-1)*P(addSpecies)*P(add specific species) 

        ######## Accept move #######
        estimNew[i,]<-estimNew[i,]                                  ### if move is accepted, record new logL
        speciestoUse[[1]]<-speciestoUse[[2]]                        ### if move is accepted, tentative species becomes present species.
        abundUsedSpecies[[1]]<-abundUsedSpecies[[2]]                        ###       -->>-->>--   , abundances of new set of species are kept  
        cat('Present species become:', speciestoUse[[1]], '\n')

      }  else {                                                    ######### Reject Move ########
        estimNew[i,]<-estimNew[i-1,]                               ### if move is rejected, record previous logL
        speciestoUse[[1]]<-speciestoUse[[1]]                       ### if move is rejected, keep present species as it is.
        abundUsedSpecies[[1]]<-abundUsedSpecies[[1]]
        cat('Present species remain:', speciestoUse[[1]], '\n')    ### temporary
      }

    }    ### close  if(stepAdd==TRUE)

    
    else {   ########################################  type of move was to remove species   ##############
      removeProba<-(1/abundUsedSpecies[[1]])/sum(1/abundUsedSpecies[[1]])   ###### sampling weight inversely proportional to size. Do this for present Species
      candidateRemove<-removeProba[as.character(record[i,3])]
      candidateAdd<-allSpecies[allSpecies[,1]==(as.character(record[i,3])),2]

      
      if (runif(1) < min( 1, exp(estimNew[i] - estimNew[i-1] + log(addStep) - log(removeStep) + log(candidateAdd) - log(candidateRemove) ))) {############# Accept "remove species"        
        estimNew[i,]<-estimNew[i,]
        speciestoUse[[1]]<-speciestoUse[[2]]                       ### tentative species becomes present species.
        abundUsedSpecies[[1]]<-abundUsedSpecies[[2]]
        cat('Present species become:', speciestoUse[[1]], '\n')    
        
      }  else {                                                    ######## Reject Move ###########
        estimNew[i,]<-estimNew[i-1,]
        speciestoUse[[1]]<-speciestoUse[[1]]
        abundUsedSpecies[[1]]<-abundUsedSpecies[[1]]
        cat('Present species remain:', speciestoUse[[1]], '\n')          
      }

    }   ###close else (i.e (stepAdd==FALSE))
    
#################################edw tha kanw attempt to direct swap between slaves, without going through master
       
    if( i%%exchangeInterval == 0 ){   ### every nth (2nd for now) iteration
      message("\n\nTime to attempt an exchange")      
      oddFlag<-oddFlag+1
      swap<-0
      estim.current<-estimNew[i,]/temper[ind]    ######### need untempered logL

      
########################################### CREATE prime tags for object to send around
      allowedLength<-175
      Nsubobjects<-round(length(abundUsedSpecies[[1]])/allowedLength)+1
      object.ids <- list.integers[ (noChains+1):(noChains + 4 + Nsubobjects) ]    ### 4 objcts for logL, swap message, untempered , species PLUS as many as necessary for abundances            

      if (ind%%2 == oddFlag%%2) {  ###when oddFlag zero , the following code concerns even-numbered slaves. For oddFlag 1, it concerns odd-numbered slaves.                 
        ind.partner<-ind+1

        if (0<ind.partner && ind.partner<(noChains+1)){           
          estim.partner<-mpi.recv.Robj(ind.partner, tag=object.ids[1]*node.ids[ind.partner])  #### receive the logdensity of above  partner
          message("I received the untempered: ", estim.partner)
          swaps.attempted<-swaps.attempted+1          
          lalpha<-(estim.partner - estim.current)*(temper[ind] - temper[ind.partner] )
          message("This is the acceptance probability: ", min(1,exp(lalpha)))
          
          if (runif(1)< min(1, exp(lalpha))) {    ############# exp((chain2 - chain1)*(T1 - T2))
            message("I exchanged the values")               
            swap<-1         
          }   ## end of if runif(1)<lalpha M-H step

          else {message("I didn't exchange the values")}
          
          mpi.send.Robj(obj=swap,dest=ind.partner,tag=object.ids[2]*node.ids[ind])
          message("I send message swap: ", swap)
        }     ### end of   if (0<ind.partner && ind.partner<(noChains+1)){          

        if(swap==1){
          swaps.accepted<-swaps.accepted+1

          mpi.send.Robj(obj=estim.current,dest=ind.partner,tag=object.ids[3]*node.ids[ind])
          message("I just sent the untempered: ", estim.current)
          species.swap<-speciestoUse[[1]]
          mpi.send.Robj(species.swap,dest=ind.partner,tag=object.ids[4]*node.ids[ind])
          speciestoUse[[1]]<-mpi.recv.Robj(ind.partner,tag=object.ids[4]*node.ids[ind.partner])
          
          abund.swap<-abundUsedSpecies[[1]]
          mpi.send.Robj(abund.swap,dest=ind.partner,tag=object.ids[5]*node.ids[ind])
          abundUsedSpecies[[1]]<-mpi.recv.Robj(ind.partner,tag=object.ids[5]*node.ids[ind.partner])
          

          message("I received the following abundances from: ", ind.partner, " here ", ind)
          message("The new logdensity will be the one of the neighbor but tempered with my temperature: ", estim.partner*temper[ind])
          estimNew[i,]<-estim.partner * temper[ind]
        }     ### end of if(swap==1)         
        
      } else {  ##### ###when oddFlag zero , the following code concerns odd-numbered slaves. For oddFlag 1, it concerns even-numbered slaves. I.e say what partners should do
        ind.partner<-ind-1;           
        if(0<ind.partner && ind.partner<(noChains+1)){
          mpi.send.Robj(obj=estim.current,dest=ind.partner,tag=object.ids[1]*node.ids[ind])
          message("I just sent the untempered: ", estim.current)
          swap<-mpi.recv.Robj(ind.partner,tag=object.ids[2]*node.ids[ind.partner])
          message("I received the swap message: ", swap)
          swaps.attempted<-swaps.attempted+1
        }  ####end of  if(0<ind.partner && ind.partner<(noChains+1)){
          
                  
        if(swap==1){
          swaps.accepted<-swaps.accepted+1
          estim.partner<-mpi.recv.Robj(ind.partner, tag=object.ids[3]*node.ids[ind.partner] )  #### receive the logdensity of above  partner
          message("I received the untempered: ", estim.partner)

          species.swap<-speciestoUse[[1]]        

          speciestoUse[[1]]<-mpi.recv.Robj(ind.partner,tag=object.ids[4]*node.ids[ind.partner])           
          mpi.send.Robj(species.swap,dest=ind.partner,tag=object.ids[4]*node.ids[ind])
          
          abund.swap<-abundUsedSpecies[[1]]
          abundUsedSpecies[[1]]<-mpi.recv.Robj(ind.partner,tag=object.ids[5]*node.ids[ind.partner])
          mpi.send.Robj(abund.swap,dest=ind.partner,tag=object.ids[5]*node.ids[ind])
                    
          message("I received the following species from: ", ind.partner, " here ", ind)
          message("The new logdensity will be the one of the neighbor but tempered with my temperature: ", estim.partner*temper[ind])
          estimNew[i,]<-estim.partner * temper[ind]          
        }
      }
    }

############ Use indicator variable for species presence per iteration
    record[i, speciestoUse[[1]]]<-ifelse(speciestoUse[[1]] %in% colnames(record), 1, 0)
    record[i,(5+lenSp)]<-estimNew[i,]             ###last colum is logL
    
  }    ###end of for loop for internal iterations
  
  resultSC<-list("estimNew"=estimNew, "record"=record, "usedSp"=speciestoUse[[1]], "abundUsedSp" = abundUsedSpecies[[1]], "swaps.attempt"=swaps.attempted, "swaps.accept"=swaps.accepted, "readSupport"=readSupport, "lpenalty"=lpenalty)
  return(resultSC)  
   
}    ##end of singleChain function


##########################################----------------------- MAIN ---------------------------
###send function to slaves
mpi.bcast.Robj2slave(singleChain)

### Now call the single chain function
message('Start parallel tempering')
result<-mpi.remote.exec(singleChain(TotalIter, exchangeInterval))
EndTime<-Sys.time()
duration<-EndTime-StartTime
message('PT finished in ', duration, ". ", nrow(ordered.species), " species were explored in ", ExternIter, "x", noChains, " iterations.")

step3<-list("result"=result, "duration"=duration)
  
  if (!is.null(outDir)) {
    step3.name <- paste(outDir, "/step3.RData", sep = "")
    save(step3, file=step3.name)
    rm(list= ls()[!ls() %in% c("step3")])
    gc()    
  } else {
    rm(list= ls()[!ls() %in% c("step3")])
    gc()
  }

  
mpi.close.Rslaves(dellog=FALSE)
#mpi.quit()

return(step3)
}

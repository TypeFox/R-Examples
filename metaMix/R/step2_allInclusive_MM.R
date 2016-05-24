################################################ Perform approximating EM ###############################################
#' @name reduce.space
#' @title Reduce the space of potential species by fitting the mixture model with  all potential species as categories
NULL

#' @rdname reduce.space
#' @title reduce.space
#' @description  Having the generative probabilities from step1 (generative.prob() or generative.prob.nucl()), we could proceed directly with the PT MCMC to explore the state space. Typically the number of total potential species is large. Therefore we reduce the size of the state-space, by decreasing the number of species to the low hundreds. We achieve this by fitting a  Mixture Model with as many categories as  all the potential species. Post fitting, we retain only the species categories that are not empty, that is categories that have at least one read assigned to them.
#' @param step1 list. The output from generative.prob() (or generative.prob.nucl(), that is the first step of the pipeline. Alternatively, it can be a character string containing the path name of the ".RData" file where step1 list was saved.
#' @param read.cutoff  numeric vector. This is the used to decide which species to retain for the subsequent MCMC exploration. Default value is 1, i.e keep all species that have at least one read assigned to them. If this number is still in the low thousands as opposed to the low hundreds the user may set this to a higher number, such as 10.
#' @param EMiter  Number of iterations for the EM algorithm. Default value is 500.
#' @param seed Optional argument that sets the random seed (default is 1) to make results reproducible.
#' @return step2: A list with six elements. The first one (ordered.species) is a data.frame containing all the non-empty species categories, as decided by the all inclusive mixture model, ordered by the number of reads assigned to them. The second one (pij.sparse.mat) is a sparse matrix with the generative probability between each read and each species. read.weights, gen.prob.unknown, outDir are all carried forward from the "step1" object. Finally outputEM which records the species abundances throughout the EM iterations (not used in step3 and step4).
#' @keywords reduce.space
#' @export reduce.space
#' @import Matrix data.table
#' @importFrom gtools rdirichlet
#' @examples
#' ## See vignette for more details.
#'
#' \dontrun{
#' # Either load the object created by previous step
#' data(step1)  ## example output of step1, i.e generative.prob()
#' step2 <- reduce.space(step1=step1)
#'
#' # or alternatively point to the location of the step1.RData object
#' step2 <- reduce.space(step1="/pathtoFile/step1.RData")
#' }
######################################################################################################################

reduce.space = function(step1, read.cutoff=1, EMiter=500, seed=1){

  if (is.character(step1)) {     
    load(step1)
  }
  
  should.be.in.the.list <- c("pij.sparse.mat", "ordered.species", "read.weights", "outDir", "gen.prob.unknown")

  if (sum (!( should.be.in.the.list %in% names(step1)) ) > 0) {
    message('Missing the following arguments')
    print(names(step1)[!(should.be.in.the.list %in% names(step1))] )
    stop()    
  }  else {

    reduce.space.wrapped = function(pij.sparse.mat=step1$pij.sparse.mat, ordered.species=step1$ordered.species, read.weights=step1$read.weights, outDir=step1$outDir, gen.prob.unknown=step1$gen.prob.unknown, read.cutoff.internal=read.cutoff, EMiter.internal=EMiter, seed.internal=seed){

      set.seed(seed.internal);
      tentative.species<-colnames(pij.sparse.mat)
      noSpecies<-length(tentative.species)

      hyperP<-rep(1, noSpecies)
      startW<-rdirichlet(1, hyperP)
  
      outputEM<-EM(pij = pij.sparse.mat, iter = EMiter.internal, species = tentative.species, abund = startW, readWeights = read.weights)  ### EM function

      message("EM done")
    
      approxSpecies0<-names(which(round(colMeans(outputEM$abundances[EMiter.internal,])*sum(read.weights[,"weight"]))>0))
      approxSpecies0<-approxSpecies0[-1]

      approxPij<-pij.sparse.mat[, approxSpecies0]
      
      approxSpecies.with.counts<-round(colMeans(outputEM$abundances[EMiter.internal,2:length(colnames(outputEM$abundances))])*sum(read.weights[,"weight"]))
      ordered.approx.species<-cbind(approxSpecies.with.counts, approxSpecies.with.counts/sum(approxSpecies.with.counts))
      colnames(ordered.approx.species)<-c( "countReads", "samplingWeight")
      ordered.approx.species<-  data.frame("taxonID"=rownames(ordered.approx.species), ordered.approx.species, stringsAsFactors=FALSE)
      ordered.species<-ordered.approx.species[order(-ordered.approx.species[,2]) , ]                      #### order them by read count
      ordered.species<-ordered.species[which(ordered.species$countReads>=read.cutoff.internal),]    ###potential species are the ones that have at least one read assigned to them 

      if (!("unknown" %in% ordered.species$taxonID)==T){ordered.species<-rbind(ordered.species, c("unknown", 0, 0))}
      ordered.species$countReads<- as.numeric(ordered.species$countReads)
      ordered.species$samplingWeight<- as.numeric(ordered.species$samplingWeight)
  
      ordered.species<- ordered.species[-which(ordered.species$taxonID=="unknown"),]
      
      approxSpecies<-ordered.species$taxonID
      pij.sparse.mat<-pij.sparse.mat[,approxSpecies]


## ###Flattening the sampling probabilities
      percentiles<-quantile(ordered.species$samplingWeight,  probs=c(0.2, 0.8))
      ordered.species$samplingWeight[which(ordered.species$samplingWeight  >= percentiles["80%"])] <- percentiles["80%"]
      ordered.species$samplingWeight[which(ordered.species$samplingWeight  <= percentiles["20%"])] <- percentiles["20%"]


### remove objects
      step2<-list("outputEM"=outputEM,  "pij.sparse.mat"=pij.sparse.mat,  "ordered.species"=ordered.species, "read.weights"=read.weights, "outDir"=outDir, "gen.prob.unknown"=gen.prob.unknown)

      if (!is.null(outDir)) {
        step2.name <- paste(outDir, "/step2.RData", sep = "")
        save(step2, file=step2.name)
        rm(list= ls()[!ls() %in% c("step2")])
        gc()
      } else {
        rm(list= ls()[!ls() %in% c("step2")])
        gc()
      }

    
      return(step2)
    }
    reduce.space.wrapped()
  }
  
}

#' @rdname reduce.space
#' @title reduce.space.explicit
#' @description  reduce.space.explicit is the same function as reduce.space but with more involved syntax.
#' @param pij.sparse.mat  sparse Matrix of generative probabilities computed by generative.prob() /  generative.prob.nucl().
#' @param ordered.species  data.frame with potential species ordered by numbers of reads matching them. Computed by generative.prob().  
#' @param read.weights  data.frame mapping each read identifier to a weight. For contigs the weight is the number of reads that were used to assemble it. For unassembled reads the weight is equal to one. 
#' @param outDir character vector holding the path to the output directory where the results are written. 
#' @param gen.prob.unknown numeric vector. This is the generative probability for the unknown category. Default value for BLASTx-analysis is 1e-06 while for BLASTn-analysis is 1e-20.
#' @keywords reduce.space.explicit
#' @export reduce.space.explicit
#' @import Matrix data.table
#' @importFrom gtools rdirichlet
##############################################################################################################################################################)

reduce.space.explicit = function(pij.sparse.mat, ordered.species, read.weights, outDir, gen.prob.unknown, read.cutoff=1, EMiter=500, seed=1){

  set.seed(seed);
  tentative.species<-colnames(pij.sparse.mat)
  noSpecies<-length(tentative.species)

  hyperP<-rep(1, noSpecies)
  startW<-rdirichlet(1, hyperP)
  
  outputEM<-EM(pij = pij.sparse.mat, iter = EMiter, species = tentative.species, abund = startW, readWeights = read.weights)  ### EM function

  message("EM done")
  
  approxSpecies0<-names(which(round(colMeans(outputEM$abundances[EMiter,])*sum(read.weights[,"weight"]))>0))
  approxSpecies0<-approxSpecies0[-1]

  approxPij<-pij.sparse.mat[, approxSpecies0]

  approxSpecies.with.counts<-round(colMeans(outputEM$abundances[EMiter,2:length(colnames(outputEM$abundances))])*sum(read.weights[,"weight"]))
  ordered.approx.species<-cbind(approxSpecies.with.counts, approxSpecies.with.counts/sum(approxSpecies.with.counts))
  colnames(ordered.approx.species)<-c( "countReads", "samplingWeight")
  ordered.approx.species<-  data.frame("taxonID"=rownames(ordered.approx.species), ordered.approx.species, stringsAsFactors=FALSE)
  ordered.species<-ordered.approx.species[order(-ordered.approx.species[,2]) , ]                      #### order them by read count
  ordered.species<-ordered.species[which(ordered.species$countReads>=read.cutoff),]    ###potential species are the ones that have at least one read assigned to them 

  if (!("unknown" %in% ordered.species$taxonID)==T){ordered.species<-rbind(ordered.species, c("unknown", 0, 0))}
  ordered.species$countReads<- as.numeric(ordered.species$countReads)
  ordered.species$samplingWeight<- as.numeric(ordered.species$samplingWeight)
  
  ordered.species<- ordered.species[-which(ordered.species$taxonID=="unknown"),]
  
  approxSpecies<-ordered.species$taxonID
  pij.sparse.mat<-pij.sparse.mat[,approxSpecies]


## ###Flattening the sampling probabilities
  percentiles<-quantile(ordered.species$samplingWeight,  probs=c(0.2, 0.8))
  ordered.species$samplingWeight[which(ordered.species$samplingWeight  >= percentiles["80%"])] <- percentiles["80%"]
  ordered.species$samplingWeight[which(ordered.species$samplingWeight  <= percentiles["20%"])] <- percentiles["20%"]


### remove objects
  step2<-list("outputEM"=outputEM,  "pij.sparse.mat"=pij.sparse.mat,  "ordered.species"=ordered.species, "read.weights"=read.weights, "outDir"=outDir, "gen.prob.unknown"=gen.prob.unknown)

  if (!is.null(outDir)) {
    step2.name <- paste(outDir, "/step2.RData", sep = "")
    save(step2, file=step2.name)
    rm(list= ls()[!ls() %in% c("step2")])
    gc()
  } else {
    rm(list= ls()[!ls() %in% c("step2")])
    gc()
  }

    
  return(step2)
}


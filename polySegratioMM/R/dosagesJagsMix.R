`dosagesJagsMix` <-
function(mcmc.mixture, jags.control, seg.ratio, chain=1,
           max.post.prob=TRUE,
           thresholds= c(0.5,0.6,0.7,0.8,0.9,0.95,0.99),
           print=FALSE, print.warning=TRUE, index.sample=20)
{

  ## Purpose: compute and return estimated dosages under specified model
  ##          using posterior probabilities derived from mcmc chains by
  ##          the proportion os samples in each dosage class

  ## Arguments:
  ## sumdata: mcmc summary produced by 'coda' preferably using 'readJagsMix'
  ## marker.names: specify marker names for subset of markers corresponding
  ##               to those with parents heterozygous for markers.
  ##               Note: these may be obtained from 'seg.ratios' function
  ##                     as rownames(seg.ratios(...)$aflp) if no markers
  ##                     removed or
  ##                     names(seg.ratios(...)$nbin) which is safer
  ## max.post.prob: if TRUE then classify dosages using maximum
  ##                posterior probability
  ## thresholds:    vector of cutoff probabilities for dosage class
  ## print:  if TRUE then print a sample of posterior probs for all
  ##                dosages
  ## index.sample:  if a single number then print a random sample of size
  ##                'index.sample', otherwise if a vector then print out those
  ##                markers using 'index.sample' as the index

  ## Values:
  ## p.dosage: matrix of posterior probabilities where columns are dosage
  ##           classes and rows are markers
  ## dosage:   matrix of allocated dosages based on posterior probabilities.
  ##           columns correspond to different 'thresholds' as specified
  ## max.post: maximum dosage posterior probabilties for each marker
  ## max.post.dosage: dosage allocated on basis of 'max.post'


  if (class(mcmc.mixture) != "segratioMCMC")
    stop("'mcmc.mixture' must be of class 'segratioMCMC'")

  if (class(jags.control) != "jagsControl")
    stop("'jags.control' must be of class 'jagsControl'")

  if (class(seg.ratio) != "segRatio") {
    stop("'seg.ratio' must be of class 'segRatio'")
  }

  n.comp <- jags.control$model$n.components
  dosage.names <- names(jags.control$model$E.segRatio$ratio)[1:n.comp]
 
  ind <-  grep("T\\[",coda::varnames(mcmc.mixture$mcmc.list))    # markers
  n.markers <- length(ind)
  marker.names <- names(seg.ratio$r)
  
  p.dosage <- matrix(0,nrow=n.markers,ncol=n.comp)
  colnames(p.dosage) <- dosage.names
  rownames(p.dosage) <- marker.names
  class.dosage <- p.dosage

  mcmc.marker <- mcmc.mixture$mcmc.list[[chain]][,ind]
  
  for (i in 1:n.comp) {
    p.dosage[,i] <- colMeans((mcmc.marker==i))
    class.dosage[,i] <- i
  }
  
  ## print examples of posterior probabilities of dosage
  if (print) {
    if (length(index.sample)==1) {
      cat("A random sample of posterior probabilities\n")
      print(p.dosage[sort(sample(1:dim(p.dosage)[1],index.sample)),])
    } else {
      cat("Posterior probabilities for specified markers\n")
      print(p.dosage[index.sample,])
    }      
  }
  
  ## classify using maximum posterior probability

  if (max.post.prob){
    max.post <- apply(p.dosage, 1, max)
    if (print) {
      cat("\nMaximum posterior probabilities for",dim(p.dosage)[1],"markers\n")
      print(summary(max.post))
      cat("\nProportion of genes classified using maximum posterior",
          "probability\n")
      print(colMeans(p.dosage==max.post))
      cat("Total proportion of markers classified:",
          sum(colMeans(p.dosage==max.post)),"\n")
    }
    no.maxs <- rowSums(p.dosage==max.post)
    class.post.prob <- no.maxs*NA
    if (max(no.maxs)>1) {
      if (print.warning){
        cat("Warning: more than one maximum for some markers\n")
        print(p.dosage[no.maxs>1,])
      }
    }
    for (i in 1:dim(p.dosage)[1]){
      if (no.maxs[i]==1) {
        class.post.prob[i] <-
          c(1:length(dosage.names))[p.dosage[i,]==max.post[i]]
      } 
    }
  }
  
  ## check to see how many markers can be classified assuming a threshold

  dosage <- matrix(NA, nrow=dim(p.dosage)[1], ncol=length(thresholds))
  colnames(dosage) <- thresholds
  rownames(dosage) <- rownames(p.dosage)
  
  for (th in thresholds) {

    if (print) {
      cat("\nProportion of genes classified assuming threshold of",
          th,"\n")
      print(colMeans(p.dosage>th))
      cat("Total proportion of markers classified:",
          sum(colMeans(p.dosage>th)),"\n")
    }
    for (i in 1:n.markers ){
      if (sum(p.dosage[i,]>th))
        dosage[i,th==thresholds] <- c(1:n.comp)[p.dosage[i,]>th]
    }
  }

  res <- list(p.dosage=p.dosage, dosage=dosage, thresholds=thresholds,
              chain=chain, index.sample=index.sample)

  if (max.post.prob) {
    res$max.post=max.post
    res$max.post.dosage=class.post.prob
  }
  res$call <- match.call()
  oldClass(res) <- "dosagesMCMC"
  
  return(res)

}


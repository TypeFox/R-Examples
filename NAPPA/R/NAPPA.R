
NAPPA <- function(data ,
                    tissueType = c("tumour","cells") , 
                    NReferenceSamples = sampleNumber ,
                    sampleNumber = ncol(data)-3 , 
                    scaleFOV = T ,
                    background.method = c("poisson","subtract","poisson.global","subtract.global","subtract.max","subtract.globalmax",
                                          "subtract.mean2sd","subtract.globalmean2sd","none") , 
                    nposcontrols = 4 , 
                    poscontrol.method = c("average","weighted.average","geometric.mean","average.prebc") , 
                    hk.method = c("shrunken.correct","shrunken.subtract","subtract","correct") , 
                    betas = NULL ,
                    hknormfactor.mean = NULL ,
                    sigmoidparameters = NULL ,
                    addconstant = 10 , 
                    imputezeroes.method = c("min","min.retro","none") ,
                    raise.low.counts = 2 ,
                    output=NULL) {
  
  background.method = match.arg(background.method)
  poscontrol.method = match.arg(poscontrol.method)
  hk.method = match.arg(hk.method)
  imputezeroes.method = match.arg(imputezeroes.method)
  tissueType = match.arg(tissueType)
  if(is.null(addconstant) | is.na(addconstant)) addconstant <- 0
  
  
  ### EXTRACT SETS OF DATA
  meta.all <- data[,1:3]
  colnames(meta.all) <- c("Class","Name","Accession")
  meta.all$Class <- gsub(" ","",meta.all$Class)
  meta.counts <- meta.all[meta.all$Class %in% c("Negative","Positive","Endogenous","Housekeeping"),]
  genes <- meta.counts$Name[meta.counts$Class=="Endogenous"]
  if(!is.null(betas) && (length(betas)!=length(genes) || any(names(betas)!=genes))) stop("Parameter 'betas' must have names the same as the genenames in the input data")  
  SampleID <- data[gsub(" ","",meta.all$Class)=="SampleID" , -(1:3)]
  FOVTarget <- as.numeric( data[meta.all$Class=="FOVCount" , 4] )
  FOV <-  as.numeric( data[meta.all$Class=="FOVCounted" , -(1:3)] )
  rawcounts.all <- sapply(data[meta.all$Class %in% c("Negative","Positive","Endogenous","Housekeeping") , -(1:3)] , as.numeric)
  nsamples <- ncol(rawcounts.all)
  
  ### STEP 1 : RAISE LOW COUNTS

  countsr.all <- pmax(rawcounts.all , raise.low.counts)
  
  ### STEP 2 : SCALE BY OBSERVED FIELDS OF VIEW
  
  if(scaleFOV) countsf.all <-  FOVTarget * sweep(countsr.all , 2 , FOV , "/") 
  else countsf.all <- countsr.all
  
  ### STEP 3 : BACKGROUND CORRECTION
  
  subtract.bg <- function(count, b) {
    if(is.na(count)) NA
    else if (count <= 20 | count <= 2 * b) 
      (count - b * (1 - (b^count)/factorial(count)/sum(b^(0:count)/factorial(0:count))))
    else (count - b)
  }
  meanbackgrounds <- apply(countsf.all[meta.counts$Class=="Negative",] , 2 , mean , na.rm=T)
  maxbackgrounds <- apply(countsf.all[meta.counts$Class=="Negative",] , 2 , max , na.rm=T)
  mean2sdbackgrounds <- meanbackgrounds + 2*apply(countsf.all[meta.counts$Class=="Negative",] , 2 , sd , na.rm=T)
  meanglobalbackground <- mean(as.matrix(countsf.all[meta.counts$Class=="Negative",]) , na.rm=T)
  maxglobalbackground <- max(countsf.all[meta.counts$Class=="Negative",] , na.rm=T)
  mean2sdglobalbackground <- meanglobalbackground + 2*sd(as.matrix(countsf.all[meta.counts$Class=="Negative",]) , na.rm=T)
  
  countsb.all <- switch(background.method , 
                        poisson = sapply( 1:nsamples , function(i) sapply(countsf.all[,i] , subtract.bg , b=meanbackgrounds[i])) , 
                        subtract = pmax( sweep(countsf.all , 2 , meanbackgrounds , "-") , 0) ,
                        poisson.global = sapply( 1:nsamples , function(i) sapply(countsf.all[,i] , subtract.bg , b=meanglobalbackground)) ,                
                        subtract.global = pmax( countsf.all - meanglobalbackground , 0) ,
                        subtract.max = pmax( sweep(countsf.all , 2 , maxbackgrounds , "-") , 0) ,  
                        subtract.globalmax = pmax( countsf.all - maxglobalbackground , 0) , 
                        subtract.mean2sd = pmax( sweep(countsf.all , 2 , mean2sdbackgrounds , "-") , 0) ,  
                        subtract.globalmean2sd = pmax( countsf.all - mean2sdglobalbackground , 0) ,
                        none = countsf.all
  )
  
  
  ### STEP 4 : POSITIVE CONTROL NORMALISATION
  
  posnormfactor <- apply( countsb.all[meta.counts$Class=="Positive" , ][1:nposcontrols,] , 2 , sum) / sum(c(128,32,8,2,0.5,0.125)[1:nposcontrols])
  posnormfactor.prebc <- apply( countsf.all[meta.counts$Class=="Positive" , ][1:nposcontrols,] , 2 , sum) / sum(c(128,32,8,2,0.5,0.125)[1:nposcontrols])
  posnormfactor.wa <- apply( sweep( countsb.all[meta.counts$Class=="Positive" , ][1:nposcontrols,] , 1 , c(128,32,8,2,0.5,0.125)[1:nposcontrols] , "/" ) , 2 , mean)
  posnormfactor.gm <- (apply( countsb.all[meta.counts$Class=="Positive" , ][1:nposcontrols,] , 2 , prod) / prod(c(128,32,8,2,0.5,0.125)[1:nposcontrols]))^(1/nposcontrols)
  pcdivide <- function(x,y) {out <- x/y; out[y==0]<-NA ; out}
  countsp.all <- switch(poscontrol.method , 
                        average = sweep(countsb.all , 2 , posnormfactor , pcdivide ) ,
                        average.prebc = sweep(countsb.all , 2 , posnormfactor.prebc , pcdivide ) ,
                        weighted.average = sweep(countsb.all , 2 , posnormfactor.wa , pcdivide ) ,
                        geometric.mean = sweep(countsb.all , 2 , posnormfactor.gm, pcdivide ) 
  )
  
  ### OPTIONAL STEP : REMOVE ZERO COUNTS : WILL ONLY BE RELEVANT IF raise.low.counts==0  
  
  countsz.all <- countsp.all
  countsz.all[countsz.all==0] <- NA
  
  ### STEP 5 : HOUSEKEEPER NORMALISATION
  
  logcountsp.hk <- log2( countsz.all[meta.counts$Class=="Housekeeping",] )
  logcountsp.gene <- log2( countsz.all[meta.counts$Class=="Endogenous",] )
  
  hknormfactor <- apply(logcountsp.hk , 2 , mean , na.rm=T)
  if(is.null(hknormfactor.mean)) hknormfactor.mean <- mean(hknormfactor[1:NReferenceSamples] , na.rm=T)
  gene_means <- apply(logcountsp.gene[ , 1:NReferenceSamples] , 1 , mean , na.rm=T)
  
  if(is.null(sigmoidparameters)) sigmoidparameters <- switch(tissueType , tumour=c(-5.11,1.46) , cells=c(-0.79,1.87))
  if(is.null(betas)) {
    betas <- 1/(1 + exp((sigmoidparameters[1]-gene_means)/sigmoidparameters[2]))
    names(betas) <- genes
  }
  else {
    if(!identical(genes , names(betas))) stop("If providing the shrinkage parameters, beta,, the names of the betas must match the gene names within the data")
  }
  logcountshk.gene <- switch(hk.method , 
                             shrunken.subtract = sapply(1:nsamples , function(i) logcountsp.gene[,i] - betas*hknormfactor[i]) ,
                             shrunken.correct = sapply(1:nsamples , function(i) logcountsp.gene[,i] - betas*(hknormfactor[i]-hknormfactor.mean)) ,
                             subtract = sweep(logcountsp.gene , 2 , hknormfactor , "-") ,
                             correct = sweep(logcountsp.gene , 2 , hknormfactor-hknormfactor.mean , "-")
  )
  
  ### OPTIONAL STEP : IMPUTE ZERO COUNTS : WILL ONLY BE RELEVANT IF raise.low.counts==0  
  
  geneexpression <- switch(imputezeroes.method , 
                           min = t(sapply(1:nrow(logcountshk.gene) , function(i) {out<-as.numeric(logcountshk.gene[i,]) ; if(sum(!is.na(out))>0) out[is.na(out)]<-min(out,na.rm=T) ; out})) ,
                           min.retro = t(sapply(1:nrow(logcountshk.gene) , function(i) {out<-as.numeric(logcountshk.gene[i,]) ; out[is.na(out)]<-min(logcountsp.gene[i,],na.rm=T) ; out})) ,
                           none = logcountshk.gene
  )
  
  
  ### PREPARE RESULTS
  
  geneexpression <- geneexpression + addconstant
  rownames(geneexpression) <- genes
  colnames(geneexpression) <- SampleID
  if(is.null(output)) out <- geneexpression 
  else {
    out <- list()
    out$GeneExpression <- geneexpression
    if("Housekeeping" %in% output || "All" %in% output) { 
      out$Housekeeping <- logcountsp.hk; colnames(out$Housekeeping) <- SampleID ; rownames(out$Housekeeping) <- meta.counts$Name[meta.counts$Class=="Housekeeping"]}
    if("HousekeepingFactor" %in% output | "All" %in% output) { 
      out$HousekeepingFactor <- hknormfactor; out$HousekeepingFactor.Mean <- hknormfactor.mean ; names(out$HousekeepingFactor) <- SampleID }
    if("Betas" %in% output || "All" %in% output) { out$Betas <- betas }
    if("Backgrounds" %in% output || "All" %in% output) { out$Backgrounds <- meanbackgrounds }
    if("PosFactor" %in% output  || "All" %in% output) { 
      out$PosFactor <- switch(poscontrol.method , average = posnormfactor , average.prebc = posnormfactor.prebc ,
                              weighted.average = posnormfactor.wa , geometric.mean = posnormfactor.gm ) }
    if("Description" %in% output || "All" %in% output) { 
      out$Description = list(scaleFOV=scaleFOV , background.method=background.method , poscontrol.method=poscontrol.method , 
                             hk.method=hk.method , imputezeroes.method=imputezeroes.method , addcontstant=addconstant , nposcontrols=nposcontrols , 
                             sigmoidparameters=sigmoidparameters)    }
    if("Steps" %in% output) out$Steps <- list(RAW=rawcounts.all , RAISE=countsr.all , FOV=countsf.all , 
                                              BG=countsb.all , POS=countsp.all , ZEROES=countsz.all , HK=logcountshk.gene ,
                                              META=meta.counts , GE=geneexpression )
  }
  return(out)
}

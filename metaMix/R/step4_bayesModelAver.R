################################################ Perform Bayesian Model Averaging ###############################################
#' @name bayes.model.aver
#' @title Bayesian Model Averaging
NULL
#' @rdname bayes.model.aver
#' @title bayes.model.aver
#' @description Perform Bayesian Model Averaging.  We concentrate on the chain with temperature=1 , i.e the untempered posterior,  to study the distribution over the model choices and perform model averaging.  We consider as present the species that have a posterior probability greater than 0.9. We then fit the mixture model with these species in order to obtain relative abundances and read classification probabilities. A tab seperated file that has a species summary is produced, as well as log-likelihood traceplots and cumulative histogram plots.
#' @param step3 list. The output from parallel.temper(), i.e the third step of the pipeline.  Alternatively, it can be a character string containing the path name of the ".RData" file  where step3 list was saved.         
#' @param step2 list. The output from reduce.space(), i.e the second step of the pipeline.  Alternatively, it can be a character string containing the path name of the ".RData" file  where step2 list was saved.         
#' @param taxon.name.map  The 'names.dmp' taxonomy names file, mapping each taxon identifier to the corresponding scientific name. It can be downloaded from  \url{ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz}
#' @param poster.prob.thr Posterior probability of presence of species threshold for reporting in the species summary. 
#' @keywords bayes.model.aver
#' @export bayes.model.aver
#' @import Matrix ggplot2 data.table
#' @importFrom gtools rdirichlet
#' @useDynLib metaMix
#' @examples
#' ## See vignette for more details
#'
#' \dontrun{
#' # Either load the object created by previous steps
#' data(step2)   ## example output of step2, i.e reduce.space()
#' data(step3)   ## example ouput of step3, i.e  parallel.temper()
#' step4<-bayes.model.aver(step2=step2, step3=step3, taxon.name.map="pathtoFile/taxon.file")
#'
#' # or alternatively point to the location of the step2.RData and step3.RData objects
#' step4<-bayes.model.aver(step2="pathtoFile/step2.RData", step3="pathtoFile/step3.RData",
#'                         taxon.name.map="pathtoFile/taxon.file")
#'
#' }                                                      
######################################################################################################################                                

bayes.model.aver = function(step2, step3,  taxon.name.map=NULL, poster.prob.thr=0.9){

  if (is.character(step2)) {
    load(step2)
  }

  if (is.character(step3)) {
    load(step3)
  }

    

  should.be.in.step3 <- c("duration", "result")

  should.be.in.the.list <- c("pij.sparse.mat",  "read.weights", "outDir", "gen.prob.unknown")

  if (sum (!( should.be.in.step3 %in% names(step3))) > 0) {
    message('Missing the following arguments')
    print(names(step3)[!(should.be.in.step3 %in% names(step3))] )
    stop()
  } else if  (sum (!( should.be.in.the.list %in% names(step2)) ) > 0) {
    message('Missing the following arguments')
    print(names(step2)[!(should.be.in.the.list %in% names(step2))] )
    stop()
  }  else {
    bayes.model.aver.wrapped<-function(result=step3$result, pij.sparse.mat=step2$pij.sparse.mat, read.weights=step2$read.weights, outDir=step2$outDir, gen.prob.unknown=step2$gen.prob.unknown, taxon.name.map.internal=taxon.name.map, poster.prob.thr.internal=poster.prob.thr){

      fast.rmultinom.weight <- function(proba.matrix, z.matrix, seed, weights) {
        return( .Call("C_rmultinom_weight", proba.matrix, z.matrix, weights, PACKAGE='metaMix') )
      }

      ..count..<-NULL  

  
      if (is.null(taxon.name.map.internal)) {
        stop("Please provide the 'names.dmp' file. It can be downloaded from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
      }
       
      nIter<- nrow(result$slave1$record)

      pijSparseUnknown<-cBind(pij.sparse.mat, "unknown"=gen.prob.unknown)


 #     message("Associate taxonIDs with scientific names: reading \"names.dmp\" can take a few minutes")
  
      names.dmp<-fread(input=taxon.name.map.internal, header=F, sep="|", select=c(1,2,4))
      testTaxon<- as.data.table(lapply(names.dmp, function(x) {gsub("\t", "", x)}))
      taxonNames<-testTaxon[which(grepl("scientific", testTaxon[["V4"]])),1:2, with=FALSE]
      setnames(x=taxonNames, old=c("taxonID", "scientName"))
    
      taxonNames<-as.data.frame(taxonNames)
  
      ofInterest<-result$slave1$record[round(nIter/5):nIter,4:(ncol(result$slave1$record)-1)]  ##burn-in 20%
      present.probabilities<- round(apply(ofInterest, MARGIN=2, function(x) sum(x)/length(x)), digits=2)
      poster.prob.all<-present.probabilities[present.probabilities>0]


      if (length(present.probabilities[present.probabilities>=poster.prob.thr.internal])>1) {  ###default                   
        poster.prob<-present.probabilities[present.probabilities>=poster.prob.thr.internal]   ###default 
      } else  if (length(present.probabilities[present.probabilities>=poster.prob.thr.internal])==1) {  
        message("Only the unknown category has posterior probability>=",poster.prob.thr.internal, ". Try post.prob>=", poster.prob.thr.internal-0.1)
        if (length(present.probabilities[present.probabilities>=(poster.prob.thr.internal-0.1)])>1) { 
          poster.prob<-present.probabilities[present.probabilities>=(poster.prob.thr.internal-0.1)]  
        } else { 
          message("Only the unknown category has posterior probability>=",poster.prob.thr.internal-0.1 ,". Try post.prob>=0.5 but be careful with the interpretation of the results")
          poster.prob<-present.probabilities[present.probabilities>=0.5] 
        }
      }


      if (length(present.probabilities[present.probabilities>0.5])==1) {                   
        stop("The method did not find any present organisms in this dataset, defining as present species with posterior probability greater than 0.5. Maybe you used the wrong reference database to annotate your sequences?")
      }
  
      poster.probM<-as.data.frame(poster.prob)
      poster.probM$taxonID<-rownames(poster.probM)
      poster.prob.final<-merge(taxonNames, poster.probM, by.y="taxonID", by.x="taxonID", all.y=T)

      
      namesBF<-rownames(poster.probM)[!(rownames(poster.probM)%in%"unknown")]
      ofInterestBF<-result$slave1$record[round(nIter/5):nIter,c(namesBF, "logL")]  ##burn-in 20%
      BayesFactor<-matrix(0, ncol=2, nrow=(length(colnames(ofInterestBF))-1) )
      BayesFactor<-as.data.frame(BayesFactor)
      colnames(BayesFactor)<-c("taxonID", "log10BF")
            
      EMiter<-10
      pij.sparse.unknown<-cBind(pij.sparse.mat, "unknown"=gen.prob.unknown)


      for ( i in 1:(length(colnames(ofInterestBF))-1) ){
        BayesFactor[i,"taxonID"]<-colnames(ofInterestBF)[i]
        if (!all(ofInterestBF[,colnames(ofInterestBF)[i]]==1)) {
          BayesFactor[i,"log10BF"]<-logmean(ofInterestBF[which(ofInterestBF[,colnames(ofInterestBF)[i]]==1),"logL"])/log(10) -  logmean(ofInterestBF[which(ofInterestBF[,colnames(ofInterestBF)[i]]==0),"logL"])/log(10)
        } else {
          maxLog<- ofInterestBF[which(ofInterestBF[,"logL"]==max(ofInterestBF[which(ofInterestBF[,i]==1),"logL"]))[1],]
          namesSp<-colnames(maxLog[which(maxLog==1)])
          tempSet<- namesSp[!(namesSp %in% BayesFactor[i, "taxonID"])]
          tentSet<- c(tempSet,"unknown")
          noSpecies<-length(tentSet)
          hyperP<-rep(1, noSpecies)
          startW<-rdirichlet(1, hyperP)
          output10Tent<-EM(pij=pij.sparse.unknown, iter=EMiter, species=tentSet, abund=startW, readWeights = read.weights)
          #lpenalty<-(computePenalty(readSupport=result$readSupport, readWeights=read.weights, pUnknown=gen.prob.unknown))/log(10)
          lpenalty<-(result$slave1$lpenalty)/log(10)
          estimator <- (output10Tent$logL[EMiter,2])/log(10) + (noSpecies * lpenalty)
          BayesFactor[i,"log10BF"]<-(maxLog[,"logL"])/log(10) -  estimator
        }
      }
            

      
      
      
      poster.prob.final[which(poster.prob.final[,"taxonID"]=="unknown"),"scientName"]<-"unknown"
  
      poster.prob<-poster.prob.final[order(-poster.prob.final[,"poster.prob"]),]

      finalSpecies<-poster.prob[,"taxonID"]
      noSpecies<-length(finalSpecies)
  
###parameters for gibbs
      hyperP<-rep(1, noSpecies)
      startW<-rdirichlet(1, hyperP)

      BurnIn<-50
      GibbsCycles<-100

#      message("Running final longer chain")
      output100<-Gibbs(pij=pijSparseUnknown, iter=GibbsCycles, species=finalSpecies, abund=startW,  hyperParam=hyperP, fast.rmultinom.weight=fast.rmultinom.weight, readWeights=read.weights)
      finalAssignments<-matrix(output100$assignments[GibbsCycles,], ncol=1, dimnames=list(colnames(output100$assignments[GibbsCycles,])))
      finalAssignmentsDF <- data.frame(taxonID=rownames(finalAssignments), finalAssignments=unlist(finalAssignments))
      finalAssignmentsDF<- finalAssignmentsDF[which(finalAssignmentsDF$taxonID!="Iter"),]


      presentSpecies<-merge(taxonNames, finalAssignmentsDF, by.y="taxonID", by.x="taxonID" ,all.y=T)
      presentSpecies[presentSpecies[,"taxonID"]=="unknown",][,"scientName"]<-"unknown"
      presentSpecies<-presentSpecies[as.numeric(order(presentSpecies[,"finalAssignments"]), decreasing=TRUE),]
  
  
      presentSpecies.allInfo.temp<-merge(presentSpecies, poster.probM, by.y="taxonID", by.x="taxonID", all.y=T)
      presentSpecies.allInfo<- merge(presentSpecies.allInfo.temp, BayesFactor, by.x ="taxonID", all.x=T)

      presentSpecies.allInfo<-presentSpecies.allInfo[order(as.numeric(presentSpecies.allInfo[,"finalAssignments"]), decreasing=TRUE),]
  
      
      summary.name <- paste(outDir, "/presentSpecies_assignedReads.tsv", sep="")
  



###Classification Probability

      noSpecies<-nrow(presentSpecies)
      mean1000<-output100$abundances[GibbsCycles,2:(noSpecies+1)]
      zij<-output100$pijs %*% diag(mean1000)
      sumProd<-rowSums(zij)
      zijFinal <- as.matrix(zij / sumProd)
      colnames(zijFinal)<-colnames(output100$pijs)
      zijFinal<-zijFinal[,presentSpecies[,"taxonID"]]



      assignedReads<-list()
      classProb<-list()
      for (i in presentSpecies[,"taxonID"]){
        assignedReads[[i]]<-rownames(output100$assignedReads[output100$assignedReads[,i]>0,])
        classProb[[i]]<-zijFinal[assignedReads[[i]], i]
      }


      scientNames<-vector()
      for (i in names(classProb)) {
        scN<-presentSpecies[presentSpecies[,"taxonID"]==i, "scientName"]
        scientNames<-append(scientNames, scN)
      }
      names(classProb)<-scientNames


  
      step4<-list("result"=result,  "pij.sparse.mat"=pijSparseUnknown,   "presentSpecies.allInfo"=presentSpecies.allInfo,  "output100"=output100,  "assignedReads"=assignedReads, "classProb"=classProb)


      if (!is.null(outDir)) {
        histograms.name<-paste(outDir, "/histograms_cdf.pdf", sep="")
        pdf(histograms.name)
        for (i in names(classProb)) {

          if ( length(classProb[[i]])>1 ){
            temp<-data.frame(read=names(classProb[[i]]), prob=classProb[[i]], stringsAsFactors=F)

            temp2<-merge(temp, read.weights, by= "read", all.x=T )

            temp3<-temp2[rep(seq_len(nrow(temp2)), temp2$weight),c("read", "prob")]  ### edw einai to provlima sth vignette


            ploti<-ggplot(temp3, aes(x=get('prob'))) + stat_bin(aes(y=..count../sum(..count..), fill = ..count../sum(..count..)), breaks=c(0,0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) ) + stat_ecdf() + labs(list(title = paste(nrow(temp3),"reads assigned to ", i), x = "Classification probability", y = "Percentage of reads")) + guides(fill=guide_legend(title="Percentage of reads")) +  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) ) + theme(plot.title = element_text(size = 12))
          #print(ploti)
          #suppressMessages(print(ploti))

            suppressMessages(suppressWarnings(print(ploti)))
          }
        }
        dev.off()


        traceplot1<-paste(outDir, "/logLikelihood_traceplot_all.pdf" ,sep="")
        pdf(traceplot1)
        plot(result$slave1$record[1:nIter,"logL"], type="l",  xlab="All iterations", ylab="Log-likelihood", main="Parallel Tempering - Coldest Chain", lwd=1.5)
        dev.off()

        traceplot2<-paste(outDir, "/logLikelihood_traceplot_80.pdf", sep="")
        pdf(traceplot2)
        plot(result$slave1$record[(nIter/5):nIter,"logL"], type="l", col="dodgerblue", xlab="Last 80% of iterations", ylab="Log-likelihood", main="Parallel Tempering - Coldest Chain", lwd=1.5)
        dev.off()

        
  #      message("Results in ", summary.name)
        write.table(presentSpecies.allInfo, summary.name, sep="\t")
    
        step4.name <- paste(outDir, "/step4.RData", sep="")
        save(step4, file=step4.name)
        rm(list= ls()[!ls() %in% c("step4")])
        gc()

      } else {
        rm(list= ls()[!ls() %in% c("step4")])
        gc()
      }

    
      return(step4) 
      
    }

    bayes.model.aver.wrapped()

  }

}

#' @rdname bayes.model.aver   
#' @title bayes.model.aver.explicit
#' @description  bayes.model.aver.explicit is the same function as bayes.model.aver with a more involved syntax.
#' @param result The list produced by parallel.temper()  (or paraller.temper.nucl()) . It holds a detailed record for each chain, what moves were proposed, which were accepted and which were rejected as well the log-likelihood through the iterations.
#' @param pij.sparse.mat see ?reduce.space
#' @param read.weights  see ?reduce.space
#' @param gen.prob.unknown  see ?reduce.space
#' @param outDir  see ?reduce.space
#' @keywords bayes.model.aver.explicit
#' @export bayes.model.aver.explicit
#' @import Matrix ggplot2 data.table
#' @importFrom gtools rdirichlet
#' @useDynLib metaMix
##########################-------------------------------------------- MAIN ---------------------------------------------------------------------------------                     

bayes.model.aver.explicit<-function(result, pij.sparse.mat, read.weights, outDir, gen.prob.unknown,  taxon.name.map=NULL, poster.prob.thr=0.9){

  fast.rmultinom.weight <- function(proba.matrix, z.matrix, seed, weights) {
    return( .Call("C_rmultinom_weight", proba.matrix, z.matrix, weights, PACKAGE='metaMix') )
  }

  ..count..<-NULL  

  
  if (is.null(taxon.name.map)) {
    stop("Please provide the 'names.dmp' file. It can be downloaded from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
  }
       
  nIter<- nrow(result$slave1$record)

  pijSparseUnknown<-cBind(pij.sparse.mat, "unknown"=gen.prob.unknown)


  message("Associate taxonIDs with scientific names: reading \"names.dmp\" could take a couple of minutes")
  
  names.dmp<-fread(input=taxon.name.map, header=F, sep="|", select=c(1,2,4))
  testTaxon<- as.data.table(lapply(names.dmp, function(x) {gsub("\t", "", x)}))
  taxonNames<-testTaxon[which(grepl("scientific", testTaxon[["V4"]])),1:2, with=FALSE]
  setnames(x=taxonNames, old=c("taxonID", "scientName"))
    
  taxonNames<-as.data.frame(taxonNames)
  
  ofInterest<-result$slave1$record[round(nIter/5):nIter,4:(ncol(result$slave1$record)-1)]  ##burn-in 20%
  present.probabilities<- round(apply(ofInterest, MARGIN=2, function(x) sum(x)/length(x)), digits=2)
  poster.prob.all<-present.probabilities[present.probabilities>0]


  if (length(present.probabilities[present.probabilities>=poster.prob.thr])>1) {  ###default                   
    poster.prob<-present.probabilities[present.probabilities>=poster.prob.thr]   ###default 
  } else  if (length(present.probabilities[present.probabilities>=poster.prob.thr])==1) {  
    message("Only the unknown category has posterior probability>=",poster.prob.thr, ". Try post.prob>=", poster.prob.thr-0.1)
    if (length(present.probabilities[present.probabilities>=(poster.prob.thr-0.1)])>1) { 
      poster.prob<-present.probabilities[present.probabilities>=(poster.prob.thr-0.1)]  
    } else { 
      message("Only the unknown category has posterior probability>=",poster.prob.thr-0.1 ,". Try post.prob>=0.5 but be careful with the interpretation of the results")
      poster.prob<-present.probabilities[present.probabilities>=0.5] 
    }
  }


    if (length(present.probabilities[present.probabilities>0.5])==1) {                   
      stop("The method did not find any present organisms in this dataset, defining as present species with posterior probability greater than 0.5. Maybe you used the wrong reference database to annotate your sequences?")
    }
  

  poster.probM<-as.data.frame(poster.prob)
  poster.probM$taxonID<-rownames(poster.probM)
  poster.prob.final<-merge(taxonNames, poster.probM, by.y="taxonID", by.x="taxonID", all.y=T)

  poster.prob.final[which(poster.prob.final[,"taxonID"]=="unknown"),"scientName"]<-"unknown"
  
  poster.prob<-poster.prob.final[order(-poster.prob.final[,"poster.prob"]),]

  finalSpecies<-poster.prob[,"taxonID"]
  noSpecies<-length(finalSpecies)
  
###parameters for gibbs
  hyperP<-rep(1, noSpecies)
  startW<-rdirichlet(1, hyperP)

  BurnIn<-50
  GibbsCycles<-100

#  message("Running final longer chain")
  output100<-Gibbs(pij=pijSparseUnknown, iter=GibbsCycles, species=finalSpecies, abund=startW,  hyperParam=hyperP, fast.rmultinom.weight=fast.rmultinom.weight, readWeights=read.weights)
  finalAssignments<-matrix(output100$assignments[GibbsCycles,], ncol=1, dimnames=list(colnames(output100$assignments[GibbsCycles,])))
  finalAssignmentsDF <- data.frame(taxonID=rownames(finalAssignments), finalAssignments=unlist(finalAssignments))
  finalAssignmentsDF<- finalAssignmentsDF[which(finalAssignmentsDF$taxonID!="Iter"),]


  presentSpecies<-merge(taxonNames, finalAssignmentsDF, by.y="taxonID", by.x="taxonID" ,all.y=T)
  presentSpecies[presentSpecies[,"taxonID"]=="unknown",][,"scientName"]<-"unknown"
  presentSpecies<-presentSpecies[as.numeric(order(presentSpecies[,"finalAssignments"]), decreasing=TRUE),]
  
  
  presentSpecies.allInfo<-merge(presentSpecies, poster.probM, by.y="taxonID", by.x="taxonID", all.y=T)
  presentSpecies.allInfo<-presentSpecies.allInfo[order(as.numeric(presentSpecies.allInfo[,"finalAssignments"]), decreasing=TRUE),]
  

  summary.name <- paste(outDir, "/presentSpecies_assignedReads.tsv", sep="")
  



###Classification Probability

  noSpecies<-nrow(presentSpecies)
  mean1000<-output100$abundances[GibbsCycles,2:(noSpecies+1)]
  zij<-output100$pijs %*% diag(mean1000)
  sumProd<-rowSums(zij)
  zijFinal <- as.matrix(zij / sumProd)
  colnames(zijFinal)<-colnames(output100$pijs)
  zijFinal<-zijFinal[,presentSpecies[,"taxonID"]]



  assignedReads<-list()
  classProb<-list()
  for (i in presentSpecies[,"taxonID"]){
    assignedReads[[i]]<-rownames(output100$assignedReads[output100$assignedReads[,i]>0,])
    classProb[[i]]<-zijFinal[assignedReads[[i]], i]
  }


  scientNames<-vector()
  for (i in names(classProb)) {
    scN<-presentSpecies[presentSpecies[,"taxonID"]==i, "scientName"]
    scientNames<-append(scientNames, scN)
  }
  names(classProb)<-scientNames


  
  step4<-list("result"=result,  "pij.sparse.mat"=pijSparseUnknown,   "presentSpecies.allInfo"=presentSpecies.allInfo,  "output100"=output100,  "assignedReads"=assignedReads, "classProb"=classProb)


  if (!is.null(outDir)) {
    histograms.name<-paste(outDir, "/histograms_cdf.pdf", sep="")
    pdf(histograms.name)
    for (i in names(classProb)) {
      temp<-data.frame(read=names(classProb[[i]]), prob=classProb[[i]], stringsAsFactors=F)
      temp2<-merge(temp, read.weights, by= "read", all.x=T )
     
      temp3<-temp2[rep(seq_len(nrow(temp2)), temp2$weight),c("read", "prob")]

#      ploti<-ggplot(temp3, aes(x=get('prob'))) + stat_bin(aes(y=..count../sum(..count..), fill = ..count../sum(..count..))) + stat_ecdf() + labs(list(title = paste(nrow(temp3),"reads assigned to ", i), x = "Classification probability", y = "Percentage of reads")) + guides(fill=guide_legend(title="Percentage of reads"))  + scale_x_continuous(breaks=c(0,0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) ) + theme(plot.title = element_text(size = 12))

      ploti<-ggplot(temp3, aes(x=get('prob'))) + stat_bin(aes(y=..count../sum(..count..), fill = ..count../sum(..count..)), breaks=c(0,0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) + stat_ecdf() + labs(list(title = paste(nrow(temp3),"reads assigned to ", i), x = "Classification probability", y = "Percentage of reads")) + guides(fill=guide_legend(title="Percentage of reads"))  +  scale_x_continuous(breaks=c(0,0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) ) + theme(plot.title = element_text(size = 12))
      #print(ploti)
      #suppressMessages(print(ploti))

      suppressMessages(suppressWarnings(print(ploti)))
    }
    dev.off()


    traceplot1<-paste(outDir, "/logLikelihood_traceplot_all.pdf", sep="")
    pdf(traceplot1)
    plot(result$slave1$record[1:nIter,"logL"], type="l",  xlab="All iterations", ylab="Log-likelihood", main="Parallel Tempering - Coldest Chain", lwd=1.5)
    dev.off()

    traceplot2<-paste(outDir, "/logLikelihood_traceplot_80.pdf", sep="")
    pdf(traceplot2)
    plot(result$slave1$record[(nIter/5):nIter,"logL"], type="l", col="dodgerblue", xlab="Last 80% of iterations", ylab="Log-likelihood", main="Parallel Tempering - Coldest Chain", lwd=1.5)
    dev.off()
    
    #message("Results in ", summary.name)
    write.table(presentSpecies.allInfo, summary.name, sep="\t")
    
    step4.name <- paste(outDir, "/step4.RData", sep="")
    save(step4, file=step4.name)
    rm(list= ls()[!ls() %in% c("step4")])
    gc()

  } else {
    rm(list= ls()[!ls() %in% c("step4")])
    gc()
  }

    
  return(step4) 

}


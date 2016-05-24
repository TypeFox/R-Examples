#SUPER IMPORTANT NOTE: These methods are fairly general, but they assume the last character in the alphabet is a gap character.
setClass("Sequences", representation(alphabet="character", nseqs="numeric"), contains="matrix")
setClass("Descriptors", representation(response="numeric", pvalues="numeric"), contains="data.frame")
setClass("MetricParams", representation(smatrix="matrix", gapOpen="numeric", gapExtension="numeric"))
setClass("MotifModel", representation(mmodel="matrix", bmodel="numeric", width="integer", seqs="Sequences", np="integer"))
setClass("MotifModelSet", representation(motifs="list"))
setClass("OOPS", representation(zmatrix="matrix"), contains="MotifModel")
setClass("SSOOPS", representation(zvector="numeric"), contains="MotifModel")
setClass("ZOOPS", representation(zmatrix="matrix", gamma="numeric", qarray="numeric"), contains="MotifModel")

bs85 <- matrix(c(5, -2, -2, -2, -1, -1, -1,  0, -2, -2, -2, -1, -2, -3, -1,  1,  0, -3, -3, -1, -2, -1, -1, -6, 
-2,  6, -1, -2, -4,  1, -1, -3,  0, -4, -3,  2, -2, -4, -2, -1, -2, -4, -3, -3, -2,  0, -2, -6, 
-2, -1,  7,  1, -4,  0, -1, -1,  0, -4, -4,  0, -3, -4, -3,  0,  0, -5, -3, -4,  4, -1, -2, -6, 
-2, -2,  1,  7, -5, -1,  1, -2, -2, -5, -5, -1, -4, -4, -2, -1, -2, -6, -4, -4,  4,  1, -2, -6, 
-1, -4, -4, -5,  9, -4, -5, -4, -5, -2, -2, -4, -2, -3, -4, -2, -2, -4, -3, -1, -4, -5, -3, -6, 
-1,  1,  0, -1, -4,  6,  2, -3,  1, -4, -3,  1,  0, -4, -2, -1, -1, -3, -2, -3, -1,  4, -1, -6, 
-1, -1, -1,  1, -5,  2,  6, -3, -1, -4, -4,  0, -3, -4, -2, -1, -1, -4, -4, -3,  0,  4, -1, -6, 
 0, -3, -1, -2, -4, -3, -3,  6, -3, -5, -5, -2, -4, -4, -3, -1, -2, -4, -5, -4, -1, -3, -2, -6, 
-2,  0,  0, -2, -5,  1, -1, -3,  8, -4, -3, -1, -3, -2, -3, -1, -2, -3,  2, -4, -1,  0, -2, -6, 
-2, -4, -4, -5, -2, -4, -4, -5, -4,  5,  1, -3,  1, -1, -4, -3, -1, -3, -2,  3, -5, -4, -2, -6, 
-2, -3, -4, -5, -2, -3, -4, -5, -3,  1,  4, -3,  2,  0, -4, -3, -2, -3, -2,  0, -5, -4, -2, -6, 
-1,  2,  0, -1, -4,  1,  0, -2, -1, -3, -3,  6, -2, -4, -2, -1, -1, -5, -3, -3, -1,  1, -1, -6, 
-2, -2, -3, -4, -2,  0, -3, -4, -3,  1,  2, -2,  7, -1, -3, -2, -1, -2, -2,  0, -4, -2, -1, -6, 
-3, -4, -4, -4, -3, -4, -4, -4, -2, -1,  0, -4, -1,  7, -4, -3, -3,  0,  3, -1, -4, -4, -2, -6, 
-1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -4, -2, -3, -4,  8, -1, -2, -5, -4, -3, -3, -2, -2, -6, 
 1, -1,  0, -1, -2, -1, -1, -1, -1, -3, -3, -1, -2, -3, -1,  5,  1, -4, -2, -2,  0, -1, -1, -6, 
 0, -2,  0, -2, -2, -1, -1, -2, -2, -1, -2, -1, -1, -3, -2,  1,  5, -4, -2,  0, -1, -1, -1, -6, 
-3, -4, -5, -6, -4, -3, -4, -4, -3, -3, -3, -5, -2,  0, -5, -4, -4, 11,  2, -3, -5, -4, -3, -6, 
-3, -3, -3, -4, -3, -2, -4, -5,  2, -2, -2, -3, -2,  3, -4, -2, -2,  2,  7, -2, -4, -3, -2, -6, 
-1, -3, -4, -4, -1, -3, -3, -4, -4,  3,  0, -3,  0, -1, -3, -2,  0, -3, -2,  5, -4, -3, -1, -6, 
-2, -2,  4,  4, -4, -1,  0, -1, -1, -5, -5, -1, -4, -4, -3,  0, -1, -5, -4, -4,  4,  0, -2, -6, 
-1,  0, -1,  1, -5,  4,  4, -3,  0, -4, -4,  1, -2, -4, -2, -1, -1, -4, -3, -3,  0,  4, -1, -6, 
-1, -2, -2, -2, -3, -1, -1, -2, -2, -2, -2, -1, -1, -2, -2, -1, -1, -3, -2, -1, -2, -1, -2, -6, 
-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,  1),
    nrow=24, ncol=24)

aabet <- c("A",  "R",  "N",  "D",  "C",  "Q",  "E",  "G",  "H",  "I",  "L",  "K",  "M",  "F",  "P",  "S",  "T",  "W",  "Y",  "V",  "B",  "Z",  "X",  "-")
colnames(bs85) <- aabet
rownames(bs85) <- aabet

default.MetricParams <- new("MetricParams", smatrix=bs85, gapOpen=-10, gapExtension=-0.2)

sdist <- function(s1, s2, params=default.MetricParams) {
  
  distance <- rep(NA, length(s1))

  i1 = "-"
  i2 = "-"
  for(i in 1:length(s1)) {

    #deal with differeing lengths
    if(s1[i] == "-") {
      i1 = "-"
    } else {
      i1 = s1[i]
    }
    if(i > length(s2) || s2[i] == "-") {
      i2 = "-"
    } else {
      i2 = s2[i]
    }

    #gaps
    if(i1 == "-" && i2 != "-") {
      if(i > 1 && s1[i - 1] == "-") {
        distance[i] = params@gapExtension
      }
      else {
        distance[i] = params@gapOpen
      }
    }
    else if(i1 != "-" && i2 == "-") {
      if(i > 1 && s2[i - 1] == "-") {
        distance[i] = params@gapExtension
      }
      else{
        distance[i] = params@gapOpen
      }
    }
    else {
      distance[i] = params@smatrix[i1,i2]
    }
    
  }

  s = sum(distance)
  return(s)

}

sHammingDist <- function(s1, s2, params) {
  distance <- sum(sapply(1:min(length(s1),length(s2)), FUN=function(x) {as.double(s1[x] != s2[x])}))
  distance <- distance + abs(length(s1) - length(s2))
  return(as.double(distance))
}

dist.Sequences <- function(seqs, method="substitution", params=default.MetricParams, ...) {
  seqs <- seqs@.Data
  
  if(method == "substitution" || method=="euclidian") {
    dist <- sdist
  }
  if(method == "hamming") {
    dist <- sHammingDist
  }
  dvec <- vector(mode="numeric", length=nrow(seqs) * (nrow(seqs) - 1) / 2)
  for(i in 1:(nrow(seqs) - 1)){
    for(j in (i + 1):nrow(seqs)) {
      dvec[nrow(seqs) * (i - 1) - i  * (i - 1) / 2 + j - i] <- dist(seqs[i,], seqs[j,], params)
    }
  }
  #Make it a dissimilarity matrix
  dvec <- max(dvec) - dvec

  attr(dvec, "Size") <- nrow(seqs)
  attr(dvec, "Labels") <- rownames(seqs)
  attr(dvec, "Diag") <- FALSE
  attr(dvec, "Upper") <- FALSE
  attr(dvec, "method") <- method
  attr(dvec, "call") <- match.call()
  class(dvec) <- "dist"
    
  return(dvec)
  
}

setGeneric("dist")
setMethod("dist", "Sequences", function(x, method="euclidian", diag=F, upper=F, p=2) dist.Sequences(x, method))
setMethod("[", signature=c("Sequences"), definition=function(x, i, j, ..., drop) {
  if(missing(j)) {
    new("Sequences", x@.Data[i,], alphabet=x@alphabet, nseqs=length(i))
  } else if(missing(i)) {
    new("Sequences", x@.Data[,j], alphabet=x@alphabet, nseqs=0)
  } else {
    new("Sequences", x@.Data[i,j], alphabet=x@alphabet, nseqs=0)
  }
} )


wdist <- function(seqmatrix, sweights=NULL, dist=sdist, params=default.MetricParams) {
  if(is.null(sweights)){
    sweights <- createSWeights(seqmatrix)
  }		      

  dmatrix <- sdist(seqmatrix, dist, params)
  dmatrix <- dmatrix * (sweights %*% t(sweights))
  return(dmatrix)
}

explode <- function(string, max, alphabet) {
  result <- strsplit(string, split="")[[1]]
  while(length(result) < max) {
    result <- c(result,"-")
  }
  result <- sapply(result, FUN=function(x){if(!(x %in% alphabet)) {"B"} else {x}})
  return(result)
}

processError <- function(e, mymessage) {
  cat(paste(mymessage, "\n"))
}

read.sequences <- function(file, header = FALSE, sep = "", quote="\"", dec=".",
                 fill = FALSE, comment.char="", alphabet=aabet) {
  
  #assumes that the first column is the sequence information

  #Try loading the file
  tryCatch(data <- read.table(file, header=header, sep=sep, quote=quote, dec=dec, fill=fill, comment.char=comment.char), error=function(e) {processError(e, "Could not read file")})
  
  #Ok, continue to spread out sequences
  t <- tryCatch(data <- toupper(as.character(data[,1])), error=function(e) {processError(e, "Non-letters in sequence file")})

  if(inherits(t, "try-error")) return
  
  nseqs <- length(unique(data))
  t <- tryCatch(maxlength <- max(sapply(data, FUN=function(x) {length(unlist(strsplit(x, split="")))})), error=function(e) {processError(e, "Failed to split sequence into characters")})

  if(inherits(t, "try-error")) return
  
 t <- tryCatch( seqmatrix <- matrix(sapply(data, FUN=function(x) {explode(x, maxlength, alphabet)}), nrow=length(data), byrow=TRUE), error=function(e) {processError(e, "Failed to process non-cannonical amino acids")})

  if(inherits(t, "try-error")) return

  t <- tryCatch(seqmatrix <- apply(seqmatrix, MARGIN=2, FUN=function(x) {sapply(x, FUN=function(y) {which(alphabet == as.character(y))})}), error=function(e) {processError(e, "Failed to convert to numerical representation")})

  if(inherits(t, "try-error")) return

  rnames <- apply(seqmatrix, MARGIN=1, FUN=function(x) {paste(alphabet[x], collapse="")})
  for(i in 1:(length(rnames) - 1)) {
    if(rnames[i] %in% rnames[(i + 1):length(rnames)]) {
    rnames[i] <- paste(rnames[i], ".", sum(rnames[i] == rnames[(i + 1):length(rnames)]), sep="")
  }
  }

  rownames(seqmatrix) <- rnames


  seqs <- new("Sequences",  seqmatrix, alphabet=alphabet, nseqs=nseqs)

  
  return(seqs)
}

read.fasta <- function(file, header = FALSE, sep = "", quote="\"", dec=".",
                 fill = FALSE, alphabet=aabet) {
  return(read.sequences(file, header, sep, quote, dec, fill, comment.char=">", alphabet))
}


#this will allow the writing of sequences to a file, utilizing R's
#built in write.table method If you pass a motifModel along with the
#sequences, then it will align and output the motifs from the
#sequences.
write.sequences <- function(seqs, motifModel = NULL, file = "", append = FALSE) {

  if(!is.null(motifModel)) {
    
    outSeqs <- matrix(0, ncol=motifModel@width, nrow=nrow(seqs))
    if(class(motifModel) == "SSOOPS") {
      startPos <- which.max(motifModel@zvector)
      
      outSeqs <- seqs[, startPos:(startPos + motifModel@width - 1)]
    } else {
      for(i in 1:nrow(seqs)) {
        startPos <- which.max(motifModel@zmatrix[i, ])
        outSeqs[i,] <- seqs[i, startPos:(startPos + motifModel@width)]
      }
    }
    seqs <- outSeqs
  } 

  
  output <- apply(seqs, MARGIN=1, FUN=function(x) {paste(seqs@alphabet[x], collapse="")})
  
  write.table(output, file=file, append=append, row.names=FALSE, col.names=FALSE)
  
}

write.fasta <- function(seqs, motifModel = NULL, file = "",  eol = "\n") {


  if(!is.null(motifModel)) {
    
    outSeqs <- matrix(0, ncol=motifModel@width, nrow=nrow(seqs))
    if(class(motifModel) == "SSOOPS") {
      startPos <- which.max(motifModel@zvector)
      
      outSeqs <- seqs[, startPos:(startPos + motifModel@width - 1)]
    } else {
      for(i in 1:nrow(seqs)) {
        startPos <- which.max(motifModel@zmatrix[i, ])
        outSeqs[i,] <- seqs[i, startPos:(startPos + motifModel@width)]
      }
    }
    seqs <- outSeqs
  }

  output <- apply(seqs, MARGIN=1, FUN=function(x) {paste(seqs@alphabet[x], collapse="")})

  for(i in 1:length(output)) {
    
    write(paste("> Sequence", i), file=file, append=T)
    write(paste(output[i]), file=file, append=T)
    
  }

}

plot.Sequences <- function(seqs, clusterNumber=3, params=default.MetricParams, distanceMatrix=dist.Sequences(seqs, params=params), clusters=aclust(distanceMatrix, clusterNumber), legendText=c(), main="") {

  #This makes it so that the coloring follows the size of the clusters. Makes plots reproducible since the
  #cluster indices from kmeans are not reproducible.
  idmap <- order(as.integer(lapply(clusters, length)))
    
  colors <- rep("black", nrow(seqs))
  shades <- hcl(h=1:clusterNumber * (360 / clusterNumber), c=90, l=70 )
  for(i in 1:length(clusters)) {
    colors[clusters[[idmap[i]]]] <- shades[i]
  }
  
#  fit <- cmdscale(distanceMatrix, eig=T, k=2)
  fit <- prcomp(distanceMatrix)
  x <- fit$x[,1]
  y <- fit$x[,2]
  plot(x,y,xlab="Component 1", ylab="Component 2", col=colors, main=main)
  if(length(legendText) > 0) {
    legend(max(x), max(y) * 1.25, legendText, col=shades,text.col="black", pch=rep(1, length(legendText)), xpd=TRUE)
  }
}

setGeneric("plot")
setMethod("plot", "Sequences", function(x,y,...) plot.Sequences(x, ...))

createSWeights <- function(seqmatrix, params=default.MetricParams) {

  wmatrix <- matrix(rep(0,nrow(params@smatrix) * ncol(seqmatrix)), nrow=nrow(params@smatrix))
  rownames(wmatrix) <- rownames(params@smatrix)
  tables <- apply(seqmatrix, MARGIN=2, FUN=function(x) table(x))
  for(i in 1:ncol(wmatrix)) {
    for(j in 1:length(tables[[i]])) {
      wmatrix[names(tables[[i]])[j], i] <-  1. / (length(tables[[i]]) * tables[[i]][j])
    }
  }

  sweights <- rep(0,nrow(seqmatrix))
  for(i in 1:nrow(seqmatrix)) {
    for(j in 1:ncol(seqmatrix)) {
      sweights[i] <- sweights[i] + wmatrix[seqmatrix[i,j], j]
    }
  }

  return(sweights)
}

motifModelSet <- function(seqs, motifNumber=NA, type="fixed", width=4, verbose=T, clusterType="kmeans", maxGuess=8, plotMain="")  {

  #can't have a motif wider than the sequences
  if(width > ncol(seqs)) {
    width <- ncol(seqs)
  }

  #check if the number of sequences is accurate
  if(seqs@nseqs == 0) {
    seqs@nseqs <- length(unique(apply(seqs@.Data, MARGIN=1, FUN=function(x) {paste(x,collapse="")})))
  }
  
  if(is.na(motifNumber)) {
    cat("Guessing cluster number, this could take a while...\n")

    #can't have more clusters than unique sequences
    if(maxGuess > seqs@nseqs) {
      maxGuess <- seqs@nseqs
    }
    
    ll <- data.frame(clusterN=1:maxGuess, logLik=rep(0, maxGuess))
    for(i in 1:maxGuess) {
      cat(paste("\r", i, "/", maxGuess))
      ll$logLik[i] <- logLik(motifModelSet(seqs, i, type, width, verbose, clusterType))
    }
    cat("\n")
    plot(ll, xlab="Cluster Number", col="black", main=plotMain)
    lines(ll, col="black", lty=1)
    motifNumber <- ll$clusterN[which.min(ll$logLik)]
  }else if(motifNumber == 1) {
    return(new("MotifModelSet", motifs=list(motifModel(seqs, type, width))))
  }

  if(motifNumber > seqs@nseqs) {
    cat(paste("Warning: must have more unique sequences than motifs. Lowering motif number to", seqs@nseqs))
    motifNumber <- seqs@nseqs
  }

  clusters <- aclust(dist(seqs), clusterNumber=motifNumber, verbose=verbose, type=clusterType)
  seqList <- vector(mode="list", length=motifNumber)
  
  for(i in 1:motifNumber) {
    seqList[[i]] <- new("Sequences", seqs[clusters[[i]],], alphabet=seqs@alphabet)
  }
  
  result <- lapply(seqList, FUN=function(x) { motifModel(x, type, width) })
  mset <- new("MotifModelSet", motifs=result)
  
  return(mset)
}

motifModel <- function(seqs, type="fixed", width=4) {

  if(width > ncol(seqs)) {
    width <- ncol(seqs)
  }

  #check if the number of sequences is accurate
  if(seqs@nseqs == 0) {
    seqs@nseqs <- length(unique(apply(seqs@.Data, MARGIN=1, FUN=function(x) {paste(x,collapse="")})))
  }
  
  if(class(seqs) == "matrix") {

    seqs <- new("Sequences", seqs, alphabet=aabet, nseqs=seqs@nseqs)
    
  }
  
  if(type == "fixed") {
    model <- factory.SSOOPS(seqs, width=width)
  }
  if(type == "variable") {
    model <- factory.OOPS(seqs, width=width)
  }
  if(type == "optional") {
    model <- factory.ZOOPS(seqs, width=width)
  }

  convergence <- 10
  while(convergence > 10**-3) {
    temp <- EMStep(model,seqs)
    convergence <- sum((model@mmodel - temp@mmodel)**2 / (model@mmodel)**2)
    model <- temp
  }
  
  return(model)
}

predict.MotifModel <- function(object, newSequences=object@seqs, ...) {

  logodds <- matrix(object@mmodel, ncol=object@width)
  for(i in 1:nrow(logodds)) {
    logodds[i,] <- log(logodds[i,] / object@bmodel[i])
  }

  result <- rep(NA, nrow(newSequences))
  for(i in 1:length(result)) {
    max <- -10**200
    for(j in 1:(ncol(newSequences) -  object@width)) {
      lsum <- 0
      if((ncol(newSequences) -  object@width) == 0) {
        for(k in 1:ncol(newSequences)) {
          if(newSequences[i,k] < length(newSequences@alphabet)) {
            lsum <- lsum + logodds[newSequences[i,k], k]
          }
          else {
            lsum <- max - 1
            break
          }
        }
      }
      else {
        for(k in j:(j + object@width - 1)) {
          if(newSequences[i,k] < length(newSequences@alphabet)) {
            lsum <- lsum + logodds[newSequences[i,k],k - j + 1]
          }
          else {
            lsum <- max - 1
            break
          }
        }
      }
      if(lsum > max) {
        max <- lsum
      }
    }
    result[i] <- max
  }

  return(result)
}


predict.MotifModelSet <- function(model, newSequences) {

  result <- rep(NA, nrow(newSequences))

  prds <- matrix(rep(NA, length(model@motifs) * nrow(newSequences)), nrow=nrow(newSequences))
  
  for(i in 1:length(model@motifs)) {
    prds[,i] <- predict(model@motifs[[i]], newSequences)
  }
  
  for(i in 1:nrow(newSequences)) {

    result[i] <- max(prds[i,])
    
  }
  
  return(result)
}

classify.MotifModelSet <- function(motifSet, newSequences, threshold = 0) {

  predictions <- lapply(motifSet@motifs, FUN=function(x) {predict(x, newSequences)})
  pmatrix <- matrix(rep(NA, nrow(newSequences) * length(motifSet@motifs)), nrow=nrow(newSequences))

  for(i in 1:length(motifSet@motifs)){
    pmatrix[,i] <- predictions[[i]]
  }
 
  classes <- apply(pmatrix, MARGIN=1, FUN=function(x) {
    if(is.na(threshold) || max(x) > threshold) which.max(x)
    else -1
  })
  
  return(classes)
}


setGeneric("classify", function(model, ...) standardGeneric("classify"))
setMethod("classify", "MotifModelSet", function(model,...) classify.MotifModelSet(model,...))

setGeneric("predict")
setMethod("predict", "MotifModel", function(object,...) predict.MotifModel(object,... ))
setMethod("predict", "MotifModelSet", function(object,...) predict.MotifModelSet(object,... ))

EM.SSOOPS.Linked <- function(fullSeqMatrix, seqids, models) {

  bcounts <- rep(1, length(models[[1]]@bmodel))
  bsum <- 0
  
  for(mindex in 1:length(models)) {

    model <- models[[mindex]]
    seqmatrix <- fullSeqMatrix[seqids[[mindex]],]@.Data
    
    m <- (ncol(seqmatrix) - model@width + 1)
    n <- nrow(seqmatrix)
    ztrial <- rep(1., m)
    mcounts <- matrix(rep(0,length(model@mmodel)), dim(model@mmodel))
    names(bcounts) <- aabet
    rownames(mcounts) <- aabet
    
    msums <- rep(0, ncol(mcounts))
    psum <- 0
    
    for(i in 1:n) {
      for(j in 1:m) {
        for(k in 1:ncol(seqmatrix)) {
          if(k >= j && k < j + model@width) {
            ztrial[j] <- ztrial[j] * model@mmodel[seqmatrix[i,k],k - j + 1]
          }
          else {
            ztrial[j] <- ztrial[j] *  model@bmodel[ seqmatrix[i,k] ]
          }
          
        }
      }
      if(sum(ztrial) > 0.) {
        ztrial <- ztrial / sum(ztrial)
      }
    }

    
    for(i in 1:n) {
      for(j in 1:m) {
        for(k in 1:ncol(seqmatrix)) {
          if(k >= j && k < j + model@width) {
            mcounts[seqmatrix[i,k],k - j + 1 ] <- mcounts[ seqmatrix[i,k], k - j + 1] +  ztrial[j]
            msums[k - j + 1] <- msums[k - j + 1] + ztrial[j]
          } else {
            bcounts[seqmatrix[i,k]] <- bcounts[seqmatrix[i,k]] + 1. - ztrial[j]
            bsum <- bsum + 1. - ztrial[j]
          }
        }
      }
    }
    
    for(i in 1:nrow(model@mmodel)) {
      for(j in 1:ncol(model@mmodel)) {
        models[[mindex]]@mmodel[i, j] <- mcounts[i,j] / msums[j]
      }
    }
    
    models[[mindex]]@zvector <- ztrial
    
  }

  #I don't think this is the correct background model, since each amino acid is represented equally in the peptide library.
#  for(i in 1:length(models)) {
#    for(j in 1:length(models[[i]]@bmodel)) {
#      models[[i]]@bmodel[j] <- bcounts[j] / bsum      
#    }
#  }

  return(models)
  
}


EMStep.SSOOPS <- function(model, seqs) {

  alphabet <- seqs@alphabet

  seqmatrix <- seqs@.Data

  m <- (ncol(seqmatrix) - model@width + 1)
  n <- nrow(seqmatrix)
  
  pseudocount <- length(alphabet)
  if(pseudocount > nrow(seqmatrix)) {
    pseudocount <- nrow(seqmatrix)
  }
  
  ztrial <- rep(1., m)
  bcounts <- rep(pseudocount / length(alphabet),length(model@bmodel))
  mcounts <- matrix(rep(pseudocount / length(alphabet),length(model@mmodel)), dim(model@mmodel))
  names(bcounts) <- alphabet
  rownames(mcounts) <- alphabet

  bsum <- pseudocount
  msums <- rep(pseudocount, ncol(mcounts))
  psum <- 0

  for(i in 1:n) {
    for(j in 1:m) {
      for(k in 1:ncol(seqmatrix)) {
	if(k >= j && k < j + model@width && seqmatrix[i,k] <= nrow(mcounts)) {
          ztrial[j] <- ztrial[j] * model@mmodel[seqmatrix[i,k],k - j + 1]
	}
	else {
          ztrial[j] <- ztrial[j] *  model@bmodel[ seqmatrix[i,k] ]
	}
        
      }
    }
    ztrial <- ztrial / sum(ztrial)
  }
  
  for(i in 1:n) {
     for(j in 1:m) {
       for(k in 1:ncol(seqmatrix)) {
         if(k >= j && k < j + model@width) {
           if(seqmatrix[i,k] <= nrow(mcounts)) {
             mcounts[seqmatrix[i,k],k - j + 1 ] <- mcounts[ seqmatrix[i,k], k - j + 1] +  ztrial[j]
             msums[k - j + 1] <- msums[k - j + 1] + ztrial[j]
           }
          } else {
  	      bcounts[seqmatrix[i,k]] <- bcounts[seqmatrix[i,k]] + 1. - ztrial[j]
              bsum <- bsum + 1. - ztrial[j]
          }
       }
 }
   }
  #I don't think this is the correct background model, since each amino acid is represented equally in the peptide library.
#  for(i in 1:length(model@bmodel)) {
#    model@bmodel[i] <- bcounts[i] / bsum
#  }
  
  for(i in 1:nrow(model@mmodel)) {
    for(j in 1:ncol(model@mmodel)) {
      model@mmodel[i, j] <- mcounts[i,j] / msums[j]
    }
  }
   
   model@zvector <- ztrial

   return(model)
}


EMStep.OOPS <- function(model, seqs) {

  alphabet <- seqs@alphabet
  seqmatrix <- seqs@.Data
  
  m <- (ncol(seqmatrix) - model@width + 1)
  n <- nrow(seqmatrix)

  pseudocount <- length(alphabet)
  if(pseudocount > nrow(seqmatrix)) {
    pseudocount <- nrow(seqmatrix)
  }
  
  ztrial <- matrix(rep(1, m * n), nrow=n, ncol=m)
  bcounts <- rep(pseudocount / length(alphabet),length(model@bmodel))
  mcounts <- matrix(rep(pseudocount / length(alphabet), length(model@mmodel)), dim(model@mmodel))
  names(bcounts) <- alphabet
  rownames(mcounts) <- alphabet

  bsum <- pseudocount
  msums <- rep(pseudocount, ncol(mcounts))


  for(i in 1:n) {
    psum <- 0
    for(j in 1:m) {
      for(k in 1:ncol(seqmatrix)) {
	if(k >= j && k < j + model@width && seqmatrix[i,k] <= nrow(mcounts)) {
            ztrial[i,j] <- ztrial[i,j] * model@mmodel[seqmatrix[i,k],k - j + 1]
	}
	else {
          ztrial[i,j] <- ztrial[i,j] * model@bmodel[ seqmatrix[i,k] ]
	}
      }
      psum <- psum + ztrial[i,j]
    }
    ztrial[i,] <- ztrial[i,] / psum
     for(j in 1:m) {
       for(k in 1:ncol(seqmatrix)) {
         if(k >= j && k < j + model@width) {
           if(seqmatrix[i,k] <= nrow(mcounts)) {
             mcounts[seqmatrix[i,k],k - j + 1 ] <- mcounts[ seqmatrix[i,k], k - j + 1] +  ztrial[i,j]
             msums[k - j + 1] <- msums[k - j + 1] + ztrial[i,j]
           }
          } else {
  	      bcounts[seqmatrix[i,k]] <- bcounts[seqmatrix[i,k]] + 1. - ztrial[i,j]
              bsum <- bsum + 1. - ztrial[i,j]
          }
       }
     }
  }
  #I don't think this is the correct background model, since each amino acid is represented equally in the peptide library.
#  for(i in 1:length(model@bmodel)) {
#    model@bmodel[i] <- bcounts[i] / bsum
# }
  for(i in 1:nrow(model@mmodel)) {
    for(j in 1:ncol(model@mmodel)) {
      model@mmodel[i, j] <- mcounts[i,j] / msums[j]
    }
   }
   
   model@zmatrix <- ztrial

   return(model)
}

EMStep.ZOOPS <- function(model, seqs) {

  alphabet <- seqs@alphabet
  seqmatrix <- seqs@.Data
 

  m <- (ncol(seqmatrix) - model@width + 1)
  n <- nrow(seqmatrix)

  pseudocount <- length(alphabet)
  if(pseudocount > nrow(seqmatrix)) {
    pseudocount <- nrow(seqmatrix)
  }
  
  ztrial <- matrix(rep(1, m * n), nrow=n, ncol=m)
  bcounts <- rep(pseudocount / length(alphabet),length(model@bmodel))
  mcounts <- matrix(rep(pseudocount / length(alphabet),length(model@mmodel)), dim(model@mmodel))
  names(bcounts) <- alphabet
  rownames(mcounts) <- alphabet

  bsum <- pseudocount
  msums <- rep(pseudocount, ncol(mcounts))


  for(i in 1:n) {
    psum = 0
    noseq <- 1
    #calcultate probability of no occurence of the motif
    for(j in 1:ncol(seqmatrix)) {
      noseq <- noseq * model@bmodel[ seqmatrix[i,j] ]
    }
    for(j in 1:m) {
      for(k in 1:ncol(seqmatrix)) {
	if(k >= j && k < j + model@width && seqmatrix[i,k] <= nrow(mcounts)) {
          ztrial[i,j] <- ztrial[i,j] * model@mmodel[seqmatrix[i,k],k - j + 1]
	}
	else {
          ztrial[i,j] <- ztrial[i,j] * model@bmodel[ seqmatrix[i,k] ]
	}
      }
      psum <- psum + ztrial[i,j]
    }

    ztrial[i,] <- ztrial[i,] * model@gamma / (psum * model@gamma + noseq * m * (1 - model@gamma))
    model@qarray[i] <- sum(ztrial[i,])

     for(j in 1:m) {
       for(k in 1:ncol(seqmatrix)) {
         if(k >= j && k < j + model@width) {
           if(seqmatrix[i,k] <= nrow(mcounts)) {
             mcounts[seqmatrix[i,k],k - j + 1 ] <- mcounts[ seqmatrix[i,k], k - j + 1] +  ztrial[i,j]
             msums[k - j + 1] <- msums[k - j + 1] + ztrial[i,j]
           }
          } else {
	    bcounts[seqmatrix[i,k]] <- bcounts[seqmatrix[i,k]] + 1. - ztrial[i,j]
            bsum <- bsum + 1. - ztrial[i,j]
          }
       }
     }
  }
  #I don't think this is the correct background model, since each amino acid is represented equally in the peptide library.
#  for(i in 1:length(model@bmodel)) {
#    model@bmodel[i] <- bcounts[i] / bsum
#  }
  for(i in 1:nrow(model@mmodel)) {
    for(j in 1:ncol(model@mmodel)) {
      model@mmodel[i, j] <- mcounts[i,j] / msums[j]
    }
   }
   
   model@zmatrix <- ztrial

   model@gamma <- sum(model@qarray) / nrow(seqmatrix)

  
   return(model)
}


setGeneric("EMStep", def = function(model, seqs,...) standardGeneric("EMStep"))
setMethod("EMStep", "OOPS", EMStep.OOPS)
setMethod("EMStep", "SSOOPS", EMStep.SSOOPS)
setMethod("EMStep", "ZOOPS", EMStep.ZOOPS)

logLik.OOPS <- function(model) {

  logLik <- 0
  for(i in 1:nrow(model@seqs)) {
    pseq <- 0
    for(j in 1:ncol(model@seqs)) {
      if(j + model@width - 1 <= ncol(model@seqs)) {
        pj <- 0
        for(k in 1:model@width) {
          pj <- pj  + model@zmatrix[i,j] * log(model@mmodel[model@seqs@.Data[i,j + k - 1], k])
        }
        pseq <- pseq + pj
      }#PUT BACKGROUND BACK IN
    }
    logLik <- logLik + pseq
  }
  names(logLik) <- NULL
  return(logLik)
}

logLik.SSOOPS <- function(model) {

  logLik <- 0
  for(i in 1:nrow(model@seqs)) {
    pseq <- 0
    for(j in 1:ncol(model@seqs)) {
      if(j + model@width - 1 <= ncol(model@seqs)) {
        pj <- 0
        for(k in 1:ncol(model@seqs)) {
          if(k >= j && k < j + model@width) {
            pj <- pj  + model@zvector[j] * log(model@mmodel[model@seqs@.Data[i,k], k -j + 1])
          }
          else {
            pj <- pj + model@zvector[j] * log(model@bmodel[model@seqs@.Data[i,k]])
          }
        }
        pseq <- pseq + pj
      }
    }
    logLik <- logLik + pseq
  }
  names(logLik) <- NULL
  
  return(logLik)
}

logLik.ZOOPS <- function(model) {

  logLik <- 0
  for(i in 1:nrow(model@seqs)) {

    pseq <- 0
    #liklihood of sequence occuring
    for(j in 1:ncol(model@seqs)) {
      if(j + model@width - 1 <= ncol(model@seqs)) {
        pj <- 0
        for(k in 1:ncol(model@seqs)){
          if(k >= j && k < j + model@width) {
            pj <- pj  + model@zmatrix[i,j] * log(model@mmodel[model@seqs@.Data[i,k], k - j + 1])
          }
          else {
            pj <- pj + model@zmatrix[i,j] * log(model@bmodel[model@seqs@.Data[i,k]])
          }
        }
        pseq <- pseq + pj
      }
    }
    logLik <- logLik + pseq
    #liklihood of peptide being background
    for(j in 1:ncol(model@seqs)) {
      logLik <- logLik + (1 - model@qarray[i]) * log(model@bmodel[model@seqs@.Data[i,k]])
    }
  }
  names(logLik) <- NULL
  return(logLik)
}

BIC.MotifModel <- function(object) {

  motif <- object
  logLik <- logLik(motif)
  k <- k + motif@np
  n <- n + motif@seqs
  
  return(-2 * logLik + k * log(n))
}

BIC.MotifModelSet <- function(object) {

  mset <- object
  logLik <- 0
  k <- 0 #Number of parameters
  n <- 0 #sample size
  for(i in 1:length(mset@motifs)) {
    logLik <- logLik + logLik(mset@motifs[[i]])
    k <- k + mset@motifs[[i]]@np
    n <- n + nrow(mset@motifs[[i]]@seqs)
  }
  
  return(-2 * logLik + k * log(n))
}

logLik.MotifModelSet <- function(object) {


  mset <- object
  logLik <- 0
  for(i in 1:length(mset@motifs)) {
    logLik <- logLik + logLik(mset@motifs[[i]])
  }
  
  return(logLik)
}

setGeneric("logLik")
setMethod("logLik", "OOPS", function(object, ...) logLik.OOPS(object,...))
setMethod("logLik", "ZOOPS", function(object, ...) logLik.ZOOPS(object, ...))
setMethod("logLik", "SSOOPS", function(object, ...) logLik.SSOOPS(object, ...))
setMethod("logLik", "MotifModelSet", function(object, ...) logLik.MotifModelSet(object, ...))
setGeneric("BIC")
setMethod("BIC", "MotifModelSet", function(object) BIC.MotifModelSet(object))
setMethod("BIC", "MotifModel", function(object) BIC.MotifModel(object))


factory.OOPS <- function(seqs, width=ncol(seqs)) {
  alphabet <- seqs@alphabet
  
  model <- new("OOPS",
               zmatrix=matrix(rep(c(1,rep(0,ncol(seqs) - width + 1)), nrow(seqs)), nrow=nrow(seqs)),
               mmodel=matrix(rep(1.0/length(alphabet),length(alphabet) * width), ncol=width),
               bmodel=rep(1.0/length(alphabet),length(alphabet)), width=as.integer(width),
               seqs=seqs)

  model@np <- length(model@zmatrix) + length(model@mmodel)
  
  names(model@bmodel) <- alphabet
  rownames(model@mmodel) <- alphabet
  return(model)
}

factory.SSOOPS <- function(seqs, width=ncol(seqs)) {
  alphabet <- seqs@alphabet

  model <- new("SSOOPS",
               zvector=rep(1./(ncol(seqs) - width + 1.), ncol(seqs) - width + 1),
               mmodel=matrix(rep(1.0/length(alphabet),length(alphabet) * width), ncol=width),
               bmodel=rep(1.0/length(alphabet),length(alphabet)),
               width=as.integer(width),
               seqs=seqs)
  
  model@np <- length(model@zvector) + length(model@mmodel)
  
  names(model@bmodel) <- alphabet
  rownames(model@mmodel) <- alphabet
  return(model)
}

factory.ZOOPS <- function(seqs, width=ncol(seqs)) {
  alphabet <- seqs@alphabet

  model <- new("ZOOPS",
               zmatrix=matrix(rep(c(1,rep(0,ncol(seqs) - width + 1)), nrow(seqs)), nrow=nrow(seqs), byrow=T),
               mmodel=matrix(rep(1.0/length(alphabet),length(alphabet) * width), ncol=width),
               bmodel=rep(1.0/length(alphabet),length(alphabet)),
               width=as.integer(width),
               gamma=0.5,
               qarray=rep(0.5,nrow(seqs)),
               seqs=seqs)
  model@np <- length(model@zmatrix) + length(model@mmodel) + length(model@qarray)
  names(model@bmodel) <- alphabet
  rownames(model@mmodel) <- alphabet
  return(model)
}


print.Sequences <- function(seqs) {
  
}

print.MotifModel <- function(x,...) {

  model <- x
  
  #find where the motif begins
  ncolz <- 0
  if(class(model) == "SSOOPS") {
    zsum <- model@zvector
    ncolz <- length(model@zvector)
  } else {
    zsum <- apply(model@zmatrix, MARGIN=2, FUN=sum)
    ncolz <- ncol(model@zmatrix)
  }
  start <- which.max(zsum)

  mnum <- 4
  cut <- 0.15

  motifMode <- matrix(rep("", model@width * mnum), nrow=mnum)
  for(i in 1:ncol(model@mmodel)) {
    csort <- order(model@mmodel[,i], decreasing=T)
    motifMode[1,i] <- rownames(model@mmodel)[csort[1]]
    for(j in 2:mnum) {
      if(model@mmodel[csort[j],i] > cut) {
        motifMode[j,i] <- rownames(model@mmodel)[csort[j]]
      }
    }
  }
  motifModeStr <- apply(motifMode, MARGIN=2, FUN=function(x) {paste("[",paste(x,collapse=""),"]", sep="")})

  seqMode <- c(rep("-", start - 1), motifModeStr, rep("-", ncolz - start ))
  cat(paste("Number of sequences:", nrow(model@seqs)))
  cat(paste("\nMotif Mode:", paste(seqMode, collapse="")))
  cat("\n\nBackground:\n")
  print(sort(model@bmodel, decreasing=T)[1:5])
  
}

MotifModel.motifString <- function(model)  {

  #Check for trivial case
  if(nrow(model@seqs) == 1) {
    return(model@seqs[1])
  }
  
  #find where the motif begins
  ncolz <- 0
  if(class(model) == "SSOOPS") {
    zsum <- model@zvector
    ncolz <- length(model@zvector)
  } else {
    zsum <- apply(model@zmatrix, MARGIN=2, FUN=sum)
    ncolz <- ncol(model@zmatrix)
  }
  start <- which.max(zsum)

  mnum <- 4
  cut <- 0.15

  motifMode <- matrix(rep("", model@width * mnum), nrow=mnum)
  for(i in 1:ncol(model@mmodel)) {
    csort <- order(model@mmodel[,i], decreasing=T)
    motifMode[1,i] <- rownames(model@mmodel)[csort[1]]
    for(j in 2:mnum) {
      if(model@mmodel[csort[j],i] > cut) {
        motifMode[j,i] <- rownames(model@mmodel)[csort[j]]
      }
    }
  }
  motifModeStr <- apply(motifMode, MARGIN=2, FUN=function(x) {paste("[",paste(x,collapse=""),"]", sep="")})

  seqMode <- c(rep("-", start - 1), motifModeStr, rep("-", ncolz - start ))
  return(paste(seqMode, collapse=""))
   
}

setGeneric("print")
setMethod("print", "MotifModel", function(x,...) print.MotifModel(x))

setGeneric("motifString", def = function(x,...) standardGeneric("motifString"))
setMethod("motifString", "MotifModel", function(x,...) MotifModel.motifString(x))


plot.MotifModelSet <- function(x,...) {

  model <- x
  clusterNumber <- length(model@motifs)

  clusters <- vector("list", clusterNumber)
  legendText <- vector("character", clusterNumber)
  
  seqs.data <- matrix(0,nrow=sum(sapply(1:clusterNumber, function(x) nrow(model@motifs[[x]]@seqs))), ncol=ncol(model@motifs[[1]]@seqs))

  counter <- 1
  for(i in 1:clusterNumber) {
    clusters[[i]] <- counter:(counter - 1 + nrow(model@motifs[[i]]@seqs))
    counter <- counter + nrow(model@motifs[[i]]@seqs)
    seqs.data[clusters[[i]],] <- model@motifs[[i]]@seqs@.Data
    legendText[i] <- motifString(model@motifs[[i]])


  }

  seqs <- new("Sequences", seqs.data,alphabet=model@motifs[[1]]@seqs@alphabet)


  colors <- rep("black", nrow(seqs))
  shades <- hcl(h=1:clusterNumber * (360 / clusterNumber), c=90, l=70 )
  for(i in 1:length(clusters)) {
    colors[clusters[[i]]] <- shades[i]
  }

  fit <- prcomp(dist(seqs))
  x <- fit$x[,1]
  y <- fit$x[,2]
  par(mar=c(5,4,2,5))
  plot(x,y,xlab="Component 1", ylab="Component 2", col=colors)

  legend(max(x) / 2, max(y) * 1.25, legendText, col=shades,text.col="black", pch=rep(1, clusterNumber), xpd=TRUE)


  
}

MotifModel.plotStartingPosition <- function(motifModel) {

  par(fg="dark gray")
  if(class(motifModel) == "SSOOPS") {
    zsum <- motifModel@zvector
  } else {
    zsum <- apply(motifModel@zmatrix, MARGIN=2, FUN=sum)
  }
  barplot(zsum / sum(zsum), names.arg=as.character(1:length(zsum)), main="", col=hcl(h=1:length(zsum) * (360 / length(zsum))), border="black")

}

MotifModel.plotPositions <- function(motifModel) {

  par(mfrow=c(ceiling(motifModel@width / 2),2), mar=c(3,3,2,2), cex=0.7)
  for(i in 1:motifModel@width) {
    barx <- barplot(motifModel@mmodel[,i], main=paste("Position", i), col=hcl(h=1:ncol(motifModel@seqs) * (360 / ncol(motifModel@seqs))), border=NA, xaxt="n", ylim=c(0,max(motifModel@mmodel[,i]) * 1.3))
    #Find the highest three bars and label those
    morder <- rev(order(motifModel@mmodel[,i]))[1:3]
    text(barx[morder], motifModel@mmodel[morder,i], rownames(motifModel@mmodel)[morder],
         pos=3, col="black")
  }
  par(mfrow=c(1,1))
}

MotifModel.plotFits <- function(motifModel) {

  
  predictions <- sort(predict(motifModel), decreasing=T)
  ncol <- 50
  colsGood <- colorRampPalette(c("white", "green", "green"))(50)
  colsBad <-  colorRampPalette(c("red", "red", "white"))(50)
  cols <- c()
  if(sum(predictions > 0) > 0) {
    cols <- c(colsGood[cut(predictions[predictions > 0], breaks=ncol, labels=F)])
  }
  if(sum(predictions < 0) > 0) {
    cols <- c(cols, colsBad[cut(predictions[predictions < 0], breaks=ncol, labels=F)])
  }
  barplot(predictions, names.arg=NULL, main="", col=cols)

}

setGeneric("plotFits", function(motifModel) standardGeneric("plotFits"))
setMethod("plotFits", "MotifModel", MotifModel.plotFits)

setGeneric("plotStartingPosition", function(motifModel) standardGeneric("plotStartingPosition"))
setMethod("plotStartingPosition", "MotifModel", MotifModel.plotStartingPosition)

setGeneric("plotPositions", function(motifModel) standardGeneric("plotPositions"))
setMethod("plotPositions", "MotifModel", MotifModel.plotPositions)

setMethod("plot", "MotifModelSet", function(x, y, ...) plot.MotifModelSet(x,...))

plot.MotifModel <- function(x,...) {


  if((class(x) == "SSOOPS" && x@width < length(x@zvector)) || (class(x) != "SSOOPS" && x@width < nrow(x@zmatrix))){
    #plot motif start
    if(class(x) == "SSOOPS") {
      zsum <- x@zvector
    } else {
      zsum <- apply(x@zmatrix, MARGIN=2, FUN=sum)
    }
    barplot(zsum, names.arg=as.character(1:length(zsum)), main=paste("Motif Starting Position"), col="gray30")
  }
  
  par(mfrow=c(3,1), mar=c(3,3,2,2), ask=TRUE)
  for(i in 1:x@width) {
    barplot(x@mmodel[,i], main=paste("Motif Position", i), col="gray30")
  }
  par(mfrow=c(1,1))
  predictions <- sort(predict(x), decreasing=T)
  ncol <- 50
  colsGood <- colorRampPalette(c("white", "green", "green"))(50)
  colsBad <-  colorRampPalette(c("red", "red", "white"))(50)
  cols <- c()
  if(sum(predictions > 0) > 0) {
    cols <- c(colsGood[cut(predictions[predictions > 0], breaks=ncol, labels=F)])
  }
  if(sum(predictions < 0) > 0) {
    cols <- c(cols, colsBad[cut(predictions[predictions < 0], breaks=ncol, labels=F)])
  }
  barplot(predictions, names.arg=NULL, main="Log Odds (Fit)", col=cols)
}

setMethod("plot", "MotifModel", function(x, y, ...) plot.MotifModel(x,...))

aclust <- function(sDistMatrix, clusterNumber, verbose=T, type="kmeans", knstart=20) {	   
  
  if(class(sDistMatrix) == "Sequences") {
    sDistMatrix <- dist(sDistMatrix)
  }

  if(type != "kmeans" && type != "agglomerative") {
    stop("type must be either \"kmeans\" or \"agglomerative\"")
  }
  
  if(type == "kmeans") {

    km <- kmeans(sDistMatrix, centers=clusterNumber, nstart=knstart)
    result <- vector(mode="list", length=clusterNumber)
    for(i in 1:clusterNumber) {
      result[[i]] <- which(km$cluster == i)
    }

    return(result)
    
  }
  
  
  if(class(sDistMatrix) == "dist") {
    N <- attr(sDistMatrix,"Size")
    temp <- matrix(rep(NA, N ** 2), nrow=N, ncol=N)
    for(i in 1:(nrow(temp) - 1)) {
      for(j in (i + 1):ncol(temp)) {
        temp[i,j] <- sDistMatrix[N * (i - 1) - i  * (i - 1) / 2 + j - i]
        temp[j,i] <- temp[i,j]
      }
    }

    sDistMatrix <- temp
    
  }

  if(clusterNumber == 1) {
    return(list(1:nrow(sDistMatrix)))
  }

  if(verbose) {
    cat("Finding clusters...\n")
  }

  cn <- nrow(sDistMatrix)

  #while cluster number is greater than clusterNumber
  clusters <- as.list(1:nrow(sDistMatrix))

  while(cn > clusterNumber && cn > 1) {

    #Find minimum element in sDistMatrix
    wmin <- findMinDistElement(sDistMatrix)

    #find which cluster that element refers to
    cA <- 0
    cB <- 0
    for(i in 1:length(clusters)) {
      if(wmin[1] %in% clusters[[i]] || wmin[2] %in% clusters[[i]]) {
        if(cA == 0) {
	  cA <- i
	}else{
	  cB <- i
	}
      }
      if(cA != 0 && cB != 0) {
       break
      }
    }

   

   #clobber rows/columns in sDistMatrix with minimum element
   ca1 <- clusters[[cA]][1]
   cb1 <- clusters[[cB]][1]
   for(i in 1:ncol(sDistMatrix)) {

     #is this a within cluster distance?
     if(!(i %in% c(clusters[[cA]], clusters[[cB]]))) {
       cmin <- cA
       cmax <- cB
       if(sDistMatrix[min(ca1,i), max(ca1,i)] < sDistMatrix[min(cb1,i), max(cb1,i)]) {
         cmin <- cB
         cmax <- cA
       }
       for(j in 1:length(clusters[[cmax]])) {
         sDistMatrix[min(clusters[[cmax]][j], i), max(clusters[[cmax]][j], i)] <- sDistMatrix[min(clusters[[cmin]][1], i), max(clusters[[cmin]][1], i)]
       }
     }
  
   }

  
    #set within cluster distance to NA
    for(i in 1:length(clusters[[cA]])) {
      for(j in 1:length(clusters[[cB]])) {
        if(clusters[[cA]][i] > clusters[[cB]][j]) {
          sDistMatrix[clusters[[cB]][j],clusters[[cA]][i]] <- NA
	} else {
	  sDistMatrix[clusters[[cA]][i],clusters[[cB]][j]] <- NA
	}
      }
    }

    #join two ID sets
  
    clusters[[cA]] <- c(clusters[[cA]], clusters[[cB]])
    clusters <- clusters[-cB]

    cn <- length(clusters)
    if(verbose) {
      cat("\r")
      cat(paste(cn, "   "))
    }

  }

  return(clusters)

}


changeClusterFormat <- function(clusterList) {
  seqnumber <- 0
  for(i in 1:length(clusterList)) {
    seqnumber <- seqnumber + length(clusterList[[i]])
  }
  clusters <- rep(0,seqnumber)
  for(i in 1:seqnumber) {
    for(j in 1:length(clusterList)) {
      if(i %in% clusterList[[j]]) {
        clusters[i] = j
        break
      }
    }
  }
  return(clusters)
}

findMinDistElement <- function(sDistMatrix) {

  minimum <- NA			
  wmin <- c(0,0)

  for(i in 1:(nrow(sDistMatrix) - 1)) {
    for(j in (i + 1):ncol(sDistMatrix)) {
      if(is.na(minimum) || (!is.na(sDistMatrix[i,j]) &&  sDistMatrix[i,j] < minimum)) {
        wmin <- c(i,j)
	minimum <- sDistMatrix[i,j]
      }
    }
  }
  return(wmin)

}

simpleDescriptors <- function(seqs, response=numeric(0), include.statistics=FALSE) {

  desc <- descriptors(seqs, response, base.frame=defaultBaseMatrix[,c("count.BasicGroups", "count.AcidicGroups", "count.PolarGroups", "count.NonPolarGroups", "count.AromaticGroups", "count.ChargedGroups",  "ALogP")], do.var=F,
  alags=c(), do.counts=F, do.mean=T, do.position=F,
  include.statistics=include.statistics, accuracy=0.001)


  return(desc)
  
}

descriptors <- function(seqs, response=numeric(0), base.frame=NA, do.var=TRUE, alags=c(1,2,3), do.mean=TRUE, do.counts=FALSE, do.position=TRUE, alphabet=seqs@alphabet, include.statistics=TRUE, accuracy=0.01) {

  
  if(include.statistics) {
    if(ncol(seqs) >= 10) {
      cat("Warning: the sequence space is very large for calculating statistics")
    }
  }
  
  if(length(base.frame) == 1 && is.na(base.frame)) {
    base.frame <- defaultBaseMatrix
  }

  base.matrix <- data.matrix(base.frame)


  #For speed, create matrix vesion of sequences
  seqmatrix <- data.matrix(seqs@.Data)
  
  
  #get number of descriptors
  base.num <- ncol(base.frame)
  multiplier <- 0
  if(!is.null(alags)) {
    multiplier <- length(alags)
  }
  if(do.var)
    multiplier <- multiplier + 1
  if(do.mean)
    multiplier <- multiplier + 1
  if(do.position)
    multiplier <- multiplier + ncol(seqmatrix)


  desc.num <- base.num * multiplier
  if(do.counts)
    desc.num <- desc.num + length(alphabet)

  #Build the names first
  dnames <- rep("",desc.num)
  index <- 1
  for(i in 1:base.num) {
    if(do.mean) {
      dnames[index] <- paste(names(base.frame)[i], 'AVG', sep=".")
      index <- index + 1
    }
    if(do.var) {
      dnames[index] <- paste(names(base.frame)[i], 'VAR', sep=".")
      index <- index + 1
    }
    if(!is.null(alags)) {
      for(j in 1:length(alags)) {
        dnames[index] <- paste(names(base.frame)[i], 'AUTO', alags[j], sep=".")
        index <- index + 1
      }
    }
    if(do.position) {
      for(j in 1:ncol(seqmatrix)) {
        dnames[index] <- paste(names(base.frame)[i], 'P',j,sep=".")
        index <- index + 1
      }
    }
  }

  if(do.counts) {
    count.start.index <- index
    for(i in 1:length(alphabet)) {
      dnames[index] <- paste('NUM', alphabet[i], sep=".")
      index <- index + 1
    }
  }
  

  
  #Now calculate the descriptors
  desc <- empty.matrix(dnames, rownames(seqmatrix))
  finished <- 0
  res <- rep(NA, base.num * multiplier)
  for(i in 1:nrow(seqmatrix)) {
    index <- 1
    for(j in 1:base.num) {
      
      cur.desc <- base.matrix[seqmatrix[i,], j]
      cur.mean <- mean(cur.desc)
      cur.var <- var(cur.desc)
      
      if(do.mean) {
        res[index] <- cur.mean
        index <- index + 1
      }
      if(do.var) {
        res[index] <- cur.var
        index <- index + 1
      }
      if(!is.null(alags)) {
        for(k in 1:length(alags)) {
          res[index] <- autocorrelation(cur.desc, cur.mean, cur.var, lag=alags[k])
          index <- index + 1
        }
      }
      if(do.position){
        for(k in 1:ncol(seqmatrix)) {
          res[index] <- cur.desc[k]
          index <- index + 1
        }
      }
    }
    if(do.counts) {
      for(j in 1:length(alphabet)) {
        res[index] <- sum(seqmatrix[i,] == j)
        index <- index + 1
      }
    }
    desc[i,] <- res
  }

  cat("\n")


  #Remove non-varying descriptors
  desc.var <- apply(desc, MARGIN=2, FUN=var)

  if(sum(which(desc.var == 0)) > 0) {
    if(do.counts) {
      if(sum(which(desc.var[1:count.start.index] == 0)) > 0) {
        desc <- desc[-which(desc.var[1:count.start.index] == 0)]
      }
    }
    else {
      desc <- desc[-which(desc.var == 0)]
    }
  }

  desc <- factory.descriptor(desc)
  desc@response=response


  #now, calculate means for each row
  if(include.statistics) {

    cat("Standard Error and variances are currently unimplemented\nCalculating Means..\n.")

    dseqs <- decoys(seqs, 500)
    ddesc <- matrix(0, nrow=nrow(dseqs), ncol=ncol(desc))
    for(i in 1:nrow(dseqs)) {
      index <- 1
      for(j in 1:base.num) {
        cat(paste("\r", format((base.num * multiplier * (i- 1) + index) * 100 / (base.num * multiplier * nrow(dseqs)),width=4, digits=3), "%        "))

        cur.desc <- base.matrix[dseqs[i,], j]
        cur.mean <- mean(cur.desc)
        cur.var <- var(cur.desc)

        if(do.mean) {
          ddesc[i, index] <- cur.mean
          index <- index + 1
        }
        if(do.var) {
          ddesc[i, index] <- cur.var
          index <- index + 1
        }
        if(!is.null(alags)) {
          for(k in 1:length(alags)) {
            ddesc[i, index] <- autocorrelation(cur.desc, cur.mean, cur.var, lag=alags[k])
            index <- index + 1
          }
        }
        if(do.position){
          for(k in 1:ncol(dseqs)) {
            ddesc[i, index] <- cur.desc[k]
            index <- index + 1
          }
        }
      }
      if(do.counts) {
        for(k in 1:length(alphabet)) {
          ddesc[i, index] <- sum(dseqs[i,] == k)
          index <- index + 1
        }
      }
    }
    cat("\n")

    #Calculate estimated p-values, the amount of overlap between the two distributions using Mann-Whitney test
    index <- 1
    for(i in 1:ncol(ddesc)) {
      if(desc.var[i] != 0) {
        x <- ddesc[,i]
        y <- desc[,i]
        p.value <- wilcox.test(x,y)$p.value
        desc@pvalues[index] <- p.value
        index <- index + 1
      }
    }    
  }
  
  
  return(desc)
}

#Calculate the autocorrelation of the given function on the data with the given lag
autocorrelation <- function(data,ef=mean(data), v=var(data), lag=1) {

  if(lag >= length(data)) {
    return(0)
  }

  if(v == 0) {
    return(0)
  }

  numer <- 0
  for(i in 1:(length(data) - lag)) {
    numer <- numer + (data[i] - ef) * (data[i + lag] - ef)
  }
  numer = numer / (length(data) - lag)

  return(numer / v)
  
}

                                    
factory.descriptor <- function(data) {

  d <- new("Descriptors", data.frame(data@.Data), row.names=rownames(data), names=colnames(data), response=numeric(0), pvalues=rep(0, ncol(data)))
  names(d@pvalues) <- colnames(data)
  return(d)
  
}

empty.df <- function(cnames, rnames, default=NA) {
  df<-data.frame(matrix(rep(default,length(cnames)*length(rnames)), nrow=length(rnames)))
  colnames(df)<-cnames
  rownames(df) <- rnames
  return(df)
}

empty.matrix <- function(cnames, rnames, default=NA) {
  df<-matrix(rep(default,length(cnames)*length(rnames)), nrow=length(rnames))
  colnames(df)<-cnames
  rownames(df) <- rnames
  return(df)
}


plot.Descriptors <- function(desc,...) {

  #correlation

  if(length(desc@response) > 0) {    

    desc.cor <- cor(desc, desc@response)
    hist(desc.cor, main="Histogram of descriptor/response correlations", col="gray30")


    desc.cor.sorted <- order(abs(desc.cor), decreasing=T)
    par(cex.axis=0.7)
    barplot(desc.cor[desc.cor.sorted[1:10]], names.arg=rownames(desc.cor)[desc.cor.sorted[1:10]], ylab="Correlation", col="gray30")
    par(cex.axis=1)
    
  }

  desc.pc <- prcomp(desc, scale=T, center=T, tol=0.01, retx=T)
  plot(desc.pc, main="Scree Plot", col="gray30")
  plot(desc.pc$x[,1], desc.pc$x[,2], xlab="PC1", ylab="PC1", main="Projection")
  biplot(desc.pc, main="Biplot")
  
}

setMethod("plot", "Descriptors", function(x, y, ...) plot.Descriptors(x, y, ...))

rbind.Sequences <- function(seq1, seq2) {
  seqs <- new("Sequences", rbind(seq1@.Data, seq2@.Data), alphabet=unique(c(seq1@alphabet, seq2@alphabet)))
  return(seqs)
}
setGeneric("rbind")
setMethod("rbind", "Sequences", function(...,deparse.level=1) rbind.Sequences(...))

decoys <- function(seqs, n=nrow(seqs)){
  choices <- seqs@alphabet
  decoys <- t(sapply(1:n, FUN=function(x) {sample(x=choices, size=ncol(seqs), replace=T)}))
  decoyNames <- apply(decoys, MARGIN=1, FUN=function(x) {paste(seqs@alphabet[x], collapse="")})
  rownames(decoys) <- decoyNames
  

  return(new("Sequences", decoys, alphabet=seqs@alphabet))
}

#This function will give the FP and FN values. Classes should be
#positive for being in the motif and negative for not.
residuals.MotifModel <- function(model, seqs=model@seqs, classes=rep(0,nrow(seqs)), threshold = 10**-3) {
  
  predictions <-predict(model, seqs)
  
  result <- list(FP=0, FN=0, accuracy=0)

  for(i in 1:length(predictions)) {
    if(predictions[i] > threshold && classes[i] < 0) result$FP <- result$FP + 1
    if(predictions[i] < threshold && classes[i] >= 0) result$FN <- result$FN + 1
  }

  if(sum(classes < 0 ) == 0) {
    result$accuracy = 1 - result$FN  / length(classes)
  }
  else {
  result$accuracy=sqrt((sum(classes >= 0) - result$FN) * (sum(classes < 0) - result$FP) / (sum(classes >= 0) * sum(classes < 0)))
   }
  
  return(result)
  
}

residuals.MotifModelSet <- function(model, seqs=model@seqs, classes=rep(0, nrow(seqs)), threshold = 0) {
  
  classes.hat <- classify(model, seqs, threshold)
  result <- list(FP=0, FN=0, ACCM=0)

  for(i in 1:length(classes.hat)) {
    if(classes.hat[i] >= 0 && classes[i] < 0) result$FP <- result$FP + 1
    if(classes.hat[i] < 0 && classes[i] >= 0) result$FN <- result$FN + 1
  }


  if(sum(classes < 0 ) == 0) {
    result$accuracy = 1 - result$FN  / length(classes)
  }
  else {
    result$accuracy=sqrt((sum(classes >= 0) - result$FN) * (sum(classes < 0) - result$FP) / (sum(classes >= 0) * sum(classes < 0)))
  }
  
  return(result)  
  
}

setGeneric("residuals")
setMethod("residuals", "MotifModel", function(object, ...) residuals.MotifModel(object, ...))
setMethod("residuals", "MotifModelSet", function(object,...) residuals.MotifModelSet(object,...))


FiveTwoCV.Sequences <- function(seqs, classes, motifNumber=1, motifModelType="fixed", width=3, verbose=F,...) {


  result <- list(FP=0, FN=0, accuracy=0)

  pclasses <- classes >= 0
  pseqs.data <- seqs@.Data[pclasses,]
  nclasses <- classes < 0
  nseqs.data <- seqs@.Data[nclasses,]
  
  
  for(i in 1:5) {

    train.p <- sample(nrow(pseqs.data), size=nrow(pseqs.data) / 2, replace=F)
    train.n <- sample(nrow(nseqs.data), size=nrow(nseqs.data) / 2, replace=F)
    validate.p <- setdiff(1:nrow(pseqs.data), train.p)
    validate.n <- setdiff(1:nrow(nseqs.data), train.n)

    tseqs <- new("Sequences", rbind(pseqs.data[train.p,], nseqs.data[train.n,]), alphabet=seqs@alphabet)
    vseqs <- new("Sequences", rbind(pseqs.data[validate.p,],nseqs.data[validate.n,]), alphabet=seqs@alphabet)
    
    mdl <- motifModelSet(tseqs, motifNumber, motifModelType, width, verbose=verbose)
    temp <- residuals(mdl, vseqs, classes=c(rep(1, length(validate.p)), rep(-1, length(validate.n))))
    result$FP <- temp$FP + result$FP
    result$FN <- temp$FN + result$FN
    result$accuracy <- temp$accuracy + result$accuracy

    mdl <- motifModelSet(vseqs, motifNumber, motifModelType, width, verbose=verbose)
    temp <- residuals(mdl, tseqs, classes=c(rep(1, length(train.p)), rep(-1, length(train.n))))
    result$FP <- temp$FP + result$FP
    result$FN <- temp$FN + result$FN
    result$accuracy <- temp$accuracy + result$accuracy
    
  }

  result$FP <- result$FP / (5 * nrow(nseqs.data))
  result$FN <- result$FN / (5 * nrow(pseqs.data))
  result$accuracy <- result$accuracy / 10
  
  return(result)
  
}

FiveTwoCV.Descriptors <- function(desc, resp=desc@response, modelFxn) {
  return(FiveTwoCV(desc@.data, resp, modelFxn))
}

FiveTwoCV.default <- function(data, resp, modelFxn) {

  
   rows <- 1:nrow(data)
   SMR <- 0
   SVR <- 0
  for(i in 1:5) {

    train <- sample(nrow(data), size=nrow(data)/2, replace=F)
    validate <- setdiff(rows,train)

    mdl <- modelFxn(formula=resp~., data=cbind(resp=resp[train],data[train,]))
    validate.yhat <- predict(object=mdl, newdata=data[validate,])
    SVR <- SVR + sum((resp[validate] - validate.yhat)^2)
    SMR <- SMR + sum((mean(resp[train]) - resp[validate])^2)
    
    mdl <- modelFxn(resp~., data=cbind(resp=resp[validate],data[validate,]))
    train.yhat <- predict(object=mdl, newdata=data[train,])
    SVR <- SVR + sum((resp[train] - train.yhat)^2)
    SMR <- SMR + sum((mean(resp[validate]) - resp[train])^2)
  }

  return(1 - SVR/SMR)
}


setGeneric("FiveTwoCV", def=function(data,...) standardGeneric("FiveTwoCV"), useAsDefault=function(data, ...) FiveTwoCV.default(data,...))
setMethod("FiveTwoCV", signature="Sequences", function(data, classes, ...) FiveTwoCV.Sequences(data, classes, ...))
setMethod("FiveTwoCV", signature="Descriptors", function(data, ...) FiveTwoCV.Descriptors(data, ...))


modelStep <- function(desc, modelFxn, fitness=FiveTwoCV.default, resp=desc@response, cols.start=NA, max=(nrow(desc) - 6) / 3, hypervalidateNum=5, verbose=TRUE){

  if(is.na(cols.start)) {
    cols.start <- colnames(desc)[order(abs(cor(desc, resp)), decreasing=T)[1]]
  }

  if(verbose) {
    cat(paste("Performing stepwise regression starting at", cols.start, "and including up to", max,"descriptors\n"))
  }
  
  desc.hyperval.index <- sample(nrow(desc),size=hypervalidateNum)
  dsub <- desc[-desc.hyperval.index,]

  f.last <- 0
  f.cur <-  0

  h.last <- 0
  h.cur <-  0

  cols <- cols.start

  while(f.cur >= f.last && length(cols) < max && h.cur >= h.last * 0.75) {
    f.last <- f.cur
    h.last <- h.cur

    #attempt to add each column and check the fitness
    cols.toadd <- NULL

    for(i in 1:(ncol(desc))) {
      cols.prop <- colnames(desc[i])

      if(sum(cols == cols.prop)== 0) {
        f.temp <- fitness(dsub[,c(cols,cols.prop)], resp[-desc.hyperval.index], modelFxn)
	if(!(is.na(f.temp) || is.null(f.temp)) && (f.temp > f.cur || f.cur == 0)){
	  f.cur <- f.temp
	  cols.toadd <- cols.prop
        }
      }
    }

    cols <- c(cols,cols.toadd)

    df <- data.frame(resp=resp[-desc.hyperval.index],dsub[,cols])

    colnames(df) <- c("resp", cols)

    mdl <- modelFxn(resp~., data=df)
    yhat <- predict(mdl,desc[desc.hyperval.index,])
    SVR <- sum((resp[desc.hyperval.index] - yhat)^2)
    SMR <- sum((mean(resp[-desc.hyperval.index]) - resp[desc.hyperval.index])^2)    
    h.cur <- SVR/SMR

    if(verbose) {
      cat(paste("Current descriptors:",paste(cols, collapse=" "), "\n"))
      cat(paste("Fitness:", f.cur, "\n"))
    }
  }

    return(cols)

}

ROC <- function(model, seqs, classes, thresholdRange=NA) {

  npos <- sum(classes >= 0)
  nneg <- sum(classes < 0)
  
  result <- data.frame(FP=0:nneg, FN=rep(NA,nneg+1))

  if(is.na(thresholdRange)) {
    preds <- predict(model, seqs)
    thresholdRange <- seq(min(preds), max(preds), (max(preds) - min(preds)) / 100)
  }
  
  print("Building ROC curve...")
  for(i in 1:length(thresholdRange)) {
    cat(paste("\r", format(100 * i/length(thresholdRange), width=4, digits=3), "%"))
    rs <- residuals(model, seqs, classes, threshold=thresholdRange[i])
    fp <- rs$FP
    fn <- rs$FN
    if(is.na(result[result[,1] == fp,2]) ||  result[result[,1] == fp,2] > fn) {
      result[result[,1] == fp, 2] = fn
    }
  }

  result <- result[!is.na(result[,2]),]

  result[,1] <- result[,1] / nneg
  result[,2] <- 1 - result[,2] / npos

  plot(result, xlab="False Positive Rate", ylab="True Positive Rate", type="l", col="red", xlim=c(0,1), ylim=c(0,1))
  lines(seq(0,1,0.1), seq(0,1,0.1), lty=2)
  
  return(result)
}



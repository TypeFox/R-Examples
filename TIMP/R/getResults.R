## all of these functions act on the output of the fitModel function to
## return various results
getC <- function(result, dataset=1, file=""){ 
  C<-getX(result, dataset=dataset, file="",lreturnC=TRUE)
  if(file!="")
    write.table(C, file=paste(file,
                        "_C_dataset_", dataset, ".txt", sep=""),
                row.names = result$currModel@modellist[[dataset]]@x,
                quote=FALSE) 
  C
}
getCLPList <- function(result, getclperr = FALSE, file="") {
  specList <- getSpecList(result$currModel, result$currTheta, getclperr)
  if(file!="")
    for(i in 1:length(specList))
      write.table(specList[[i]], file=paste(file,
                                   "_spec_dataset_", i, ".txt", sep=""),
                  row.names = result$currModel@modellist[[i]]@x2,
                  quote=FALSE) 
  specList
}
getCLP <- function(result, getclperr = FALSE, dataset=1, file=""){ 
  spec <- getSpecList(result$currModel, result$currTheta, getclperr)[[dataset]]
  if(file!="")
    write.table(spec, file=paste(file,
                        "_spec_dataset_", dataset, ".txt", sep=""),
                row.names = result$currModel@modellist[[dataset]]@x2,
                quote=FALSE) 
  spec
}
getData <- function(result, dataset = 1, weighted = FALSE) {
  if(weighted)   
    datamat <- result$currModel@data[[dataset]]@psi.weight
  else
    datamat <- result$currModel@data[[dataset]]@psi.df
  datamat
}
getDAS <- function(result, getclperr = FALSE, dataset=1, file=""){ 
  spec <- getSpecList(result$currModel, result$currTheta, getclperr)[[dataset]]
  if((result$currModel@modellist[[dataset]]@fullk==TRUE)||
  (result$currModel@modellist[[dataset]]@seqmod==TRUE))
  {A<-getX(result, dataset=dataset, file="",lreturnA=TRUE)
  ncomp<-nrow(A)
  specA<-spec[,1:ncomp]%*%t(A)
  spec[,1:ncomp]<-specA}
  if(file!="")
    write.table(spec, file=paste(file,
                        "_DAS_dataset_", dataset, ".txt", sep=""),
                row.names = result$currModel@modellist[[dataset]]@x2,
                quote=FALSE) 
  spec
}
getSVDData <- function(result, numsing = 2, dataset=1) {
  datamat <- getData(result, dataset) 
  doSVD(datamat, numsing, numsing)
}
getResiduals <- function(result, dataset = 1) {
  residlist <- result$currModel@fit@resultlist[[dataset]]@resid
  residmat <- unlist(residlist)
  dim(residmat) <- c(length(residlist[[1]]), length(residlist) )
  residmat
}
getSVDResiduals <- function(result, numsing = 2, dataset = 1) {
  residmat <- getResiduals(result, dataset) 
  doSVD(residmat, numsing, numsing)
}
getTraces <- function(result, dataset=1, file="") {
  fitted <- result$currModel@fit@resultlist[[dataset]]@fitted 
  tracesmat <- unlist(fitted)
  dim(tracesmat) <- c(length(fitted[[1]]), length(fitted) )
  if(file!="")
    write.table(tracesmat, file=paste(file,
                             "fit.txt", sep=""),
                col.names = result$currModel@modellist[[dataset]]@x2,
                row.names = result$currModel@modellist[[dataset]]@x,
                quote=FALSE) 
  tracesmat
}
getdim1 <- function(result, dataset=1) 
  result$currModel@modellist[[dataset]]@x
getdim2 <- function(result, dataset=1) 
  result$currModel@modellist[[dataset]]@x2
parEst <- function(result, param = "", dataset = NA, verbose = TRUE,
                   file = "", stderr=TRUE) {
  currTheta <- result$currTheta
  currModel <- result$currModel
  if(stderr && currModel@optlist[[1]]@sumnls) {
    stderr <- TRUE
    if(class(currModel@fit@nlsres$onls) == "timp.nls.lm") {
    currErr <- getThetaCl(sumnls(result)$coefficients[,2], currModel) 
    } else {
    currErr <- getThetaCl(sumnls(result)$parameters[,2], currModel) 
    }
  }
  else
    stderr <- FALSE
  reslist <- stderrlist <- list()
  if(param == "")
    param <- slotNames(theta()) 
  if(is.na(dataset))
    dataset <- 1:length(currTheta)
  
  for(nm in param) {
    for(j in dataset) {
      if(length( slot(currTheta[[j]], nm)) > 0) {
        if(is.null(reslist[[nm]])) {
          reslist[[nm]] <- stderrlist[[nm]] <- list()
          if(verbose) cat("Parameters:", nm, "\n", file=file)
        }
        reslist[[nm]][[length(reslist[[nm]])+1]] <- slot(currTheta[[j]], nm)
        if(stderr) {
          if(nm %in% currModel@modellist[[j]]@positivepar)
            std <- stdErrTransform(slot(currTheta[[j]], nm),
                                   slot(currErr[[j]], nm))
          else
            std <- slot(currErr[[j]], nm)
          stderrlist[[nm]][[length(stderrlist[[nm]])+1]] <- std
        }
        if(verbose) {
          cat("dataset ", j, ": ", toString(slot(currTheta[[j]], nm)),
              "\n", sep="", file=file)   
          if(stderr)
            cat("    standard errors: ",
                toString(stderrlist[[nm]][[length(stderrlist[[nm]])]]),
                "\n", sep="", file=file)
        }
      }
    }
  }
  invisible(list(reslist=reslist,stderrlist=stderrlist))
}
stdErrTransform <- function(k, err) {
   x <- rep(0,length(err))
   kl <- log(k)
   errl <- log(err)
   for(i in 1:length(k)){
     	      b1 <- abs(exp(kl[i]+errl[i]) - exp(kl[i]))
	      b2 <- abs(exp(kl[i]-errl[i]) - exp(kl[i])) 
	      x[i] <- ifelse(b1>b2,b1,b2)              
        }
   x
}


sumnls <- function(result) {
  result$currModel@fit@nlsres[[2]]
}
onls <- function(result) {
  result$currModel@fit@nlsres[[1]]
}

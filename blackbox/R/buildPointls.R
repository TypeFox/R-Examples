insert.at <- function(a, pos, ...){ ## http://stackoverflow.com/questions/18951248/insert-elements-in-a-vector-in-r
  dots <- list(...)
  stopifnot(length(dots)==length(pos))
  result <- vector("list",2*length(pos)+1)
  result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
  result[c(FALSE,TRUE)] <- dots
  unlist(result)
}

buildPointls <- function(dataFile=blackbox.getOption("dataFile"),
                         respCols=NULL,
                         subsetRows=NULL,
                         ycolname, ## this is set by this function
                         cleanResu=""
) {
  blackbox.options(ycolname=ycolname) ## !! SETS !!
  ParameterNames <- blackbox.getOption("ParameterNames")
  paramnbr <- blackbox.getOption("paramnbr") ## length(ParameterNames)
  if ("Migraine" %in% blackbox.getOption("usedBy")) {
    write(paste("Demographic model:", paste(blackbox.getOption("DemographicModel"),collapse=" ")), file=cleanResu)
    NbrMultilocusLiks <- blackbox.getOption("NbrMultilocusLiks")
    if ( ! is.null(NbrMultilocusLiks)) {stop.redef("'NbrMultilocusLiks' is obsolete")}
    if (is.null(respCols)) respCols <- blackbox.getOption("respCols") ## may still be NULL
    usernames <- sapply(ParameterNames, userunit, format="ASCII")
    write(paste("Canonical parameters:", paste(usernames, collapse=" ")), file=cleanResu)
    write("* N stands for number of gene copies, \n    i.e. 2N = 4 x [the number of diploid individuals] *\n", file=cleanResu)
  }
  #
  pointls <- try(read.table(dataFile, row.names=NULL, header=F))
  if(inherits(pointls,"try-error")) stop.redef("Reading data file failed.")
  if (! is.null(subsetRows)) pointls <- pointls[subsetRows, ]
  #
  pointls <- pointls[ do.call(order, pointls) , ] ## do not change this order later !
  absD <- (apply(abs(diff(t(t(pointls)), lag=1)), 1, max)) ## absD no longer has rownames
  nullabsD <- (absD==0)
  if (("Migraine" %in% blackbox.getOption("usedBy")) && length(which(nullabsD))>0) {
    message.redef("(!) Some likelihood estimates  from independent replicates appear identical. ")
    message.redef("    Although this could occur in normal use, this may well be the result ")
    message.redef("    of appending twice or more the result of the same replicate to the pointls file. ")
    message.redef("    Look in particular for the following coordinates in the pointls file:")
    apply(pointls[which(nullabsD), , drop=FALSE], 1, function(v) {message.redef(v[1:paramnbr])})
    if("automatedCleaning" %innc% blackbox.getOption("miscOptions")) {
      message.redef("    Removing suspicious replicates from pointls according to the 'automatedCleaning' option.")
      pointls <- pointls[-(which(nullabsD)), ]
      #### rownames(pointls) <- as.character(nrow(pointls)) ## should not be necessary; rather keep line numbers of original pointls file
      bla <- strsplit(dataFile,".",fixed=TRUE)[[1]]
      clnst <- paste(insert.at(bla, length(bla)-1, "cln"),collapse=".")
      comment <- readLines(dataFile, n=1)
      write(comment, file=clnst)
      write("## file reconstructed from cleaned pointls", file=clnst, append=T)
      write.table(pointls, file=clnst, append=T, quote=F, row.names=F, col.names=F)
    }
  }
  lenptls <- nrow(pointls)
  if(lenptls==0) {stop.redef("(!) empty pointls file?\n")}
  blackbox.options(lenptls=lenptls)
  nFONRespCols <- ncol(pointls)-paramnbr
  ## 09/2015 new conception:
  if (is.null(respCols)) {
    pointls <- pointls[,c(seq(paramnbr),ncol(pointls))] ## has no useful names
  } else {
    if ( ! is.character(respCols) ) {
      if (max(respCols)> nFONRespCols) {
        stop.redef("(!) max(respCols)> nFONRespCols: check the 'respCols' argument")
      }
      respCols <- colnames(pointls)[paramnbr+respCols]
    }
    yvalues <- apply(pointls[, respCols, drop=FALSE], 1, sum)
    pointls <- cbind(pointls[,seq(paramnbr),drop=FALSE],yvalues)
  }
  colnames(pointls) <- c(ParameterNames,ycolname)
  blackbox.options(respCols=ycolname) ## !! SETS !!
  if ("Migraine" %in% blackbox.getOption("usedBy")) {
    localstring <- paste("\nHighest likelihood in ",dataFile," for", sep="");
  } else localstring <- paste("\nHighest response in ",dataFile," for", sep="");
  message.redef(localstring)
  userNvalues <- pointls[which.min(pointls[,ycolname]), ]
  names(userNvalues) <- sapply(names(userNvalues), formatName, format="ASCII")
  message.redef(userNvalues)
  #
  invisible(pointls)
}

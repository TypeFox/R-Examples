siarhdrs <-
function(siardata) {

if(siardata$SHOULDRUN==FALSE) {
    cat("You must load in some data first (via option 1) in order to use \n")
    cat("this feature of the program. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    return(NULL)
}

if(length(siardata$output)==0) {
    cat("No output found - check that you have run the SIAR model. \n \n")
    return(NULL)
}

cat("Summary information for the output file ... \n")

hdrsummary <- matrix(0,ncol=4,nrow=ncol(siardata$output))
colnames(hdrsummary) <- c("Low 95% hdr","High 95% hdr","mode","mean")
if(length(siardata$targets)>0) {
    sourcenames <- as.character(siardata$sources[,1])
    if(siardata$targets[1,1]%%1!=0) {
      rownames(hdrsummary) <- c(sourcenames,paste("SD",seq(1,siardata$numiso),sep=""))
    } else {
      rownames(hdrsummary) <- paste(rep(paste(c(sourcenames,paste("SD",seq(1,siardata$numiso),sep="")),"G",sep=""),times=siardata$numgroups),sort(rep(seq(1,siardata$numgroups),times=siardata$numsources+siardata$numiso)),sep="")
    }
} else {
    rownames(hdrsummary) <- colnames(siardata$output)
}


for(i in 1:ncol(siardata$output)) {
    if(all(siardata$output[,i]==0)) {
      hdrsummary[i,] <- 0
    } else {
      temp <- hdr(siardata$output[,i],h=bw.nrd0(siardata$output[,i]))
      hdrsummary[i,1] <- max(0,temp$hdr[2,1])
      hdrsummary[i,2] <- temp$hdr[2,2]
      hdrsummary[i,3] <- temp$mode
      hdrsummary[i,4] <- mean(siardata$output[,i])
    }
}

print(hdrsummary)
if(siardata$SIARSOLO==TRUE) cat("Ignore SD columns for siarsolo runs. \n")

if(any(hdrsummary>5)) {
  cat("\n")
  cat("=============== READ THIS =============== \n")
  cat("There may be some problems with this data.\n")
  cat("Some of the standard deviations seem especially large. \n")
  cat("Please check to see whether the target data lie outside \n")
  cat("the convex hull implied by the sources. \n \n")
  cat("SIAR rates the problem with this data set as: \n")
  if(any(hdrsummary>50)) {
    cat("Extremely severe - almost certainly affecting results \n")
  } else if(any(hdrsummary>20)) {
    cat("Severe - possibly severely affecting results \n")
  } else if(any(hdrsummary>5)) {
    cat("Mild - but still may affect results. \n")
  }
  cat("======================================== \n")

}

cat("Press <Enter> to continue...")
readline()
invisible()

cat("\n")
cat("Running convergence diagnostics on output. \n")
cat("Output parameters need to have been loaded in or created. \n \n")

cat("Worst parameters are ... \n")
temp <- geweke.diag(siardata$output)[[1]]
print(sort(c(pnorm(temp[temp<0]),1-pnorm(temp[temp>0])))[1:min(10,ncol(siardata$output))])

if(siardata$SIARSOLO==TRUE) cat("Ignore NAs for siar solo runs. \n \n")

cat("If lots of the p-values are very small, try a longer run of the MCMC. \n")


}

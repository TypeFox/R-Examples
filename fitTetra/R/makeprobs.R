makeprobs <-
function(marker,markername,samplenames,modelname,select,markerlines,rr,resultprobs,dat) {
  nsamp <- length(samplenames)
  probs <- data.frame(marker=rep(marker,nsamp),markername=rep(markername,nsamp), sample=samplenames,
          model=modelname, nsamp=nsamp,
          select=rep(0,nsamp), ratio=rep(NA,nsamp))
          #note that probs has a row for each sample while rr has only the selected non-NA sample ratios
  sel <- rep(select, times=ceiling(nrow(dat)/length(select)))  #extend if select shorter than data, eg if select=TRUE (default)       
  sel <- sel[dat$MarkerName==markername] #note that selection is less stringent than that of rr, so sel may be longer    
  names(sel) <- dat$SampleName[dat$MarkerName==markername] 
  m <- match(probs$sample, names(sel))
  probs$select <- as.numeric(sel[m])
  m <- match(probs$sample, names(rr)) #has NA for samples not in rr
  probs$ratio <- rr[m] #has NA for samples where m=NA
  probs <- cbind(probs, resultprobs[m,1:length(resultprobs)]) #idem
  probs
}

qcplate = function(dat, score.before="scorebefore", score.after="scoreafter", poscont="Control P", negcont="Control N", qc1.val=0.225, qc2.val=2, addcont, welltype="welltype"){
  poscont_BEFORE = mean(unlist(dat[grep(poscont, as.character(dat[[welltype]])), grep(score.before, colnames(dat))]), na.rm=TRUE)
  negcont_BEFORE = mean(unlist(dat[grep(negcont, as.character(dat[[welltype]])), grep(score.before, colnames(dat))]), na.rm=TRUE)

  addcontAFTER = sapply(addcont, function(x){mean(unlist(dat[grep(x, as.character(dat[[welltype]])), grep(score.after, colnames(dat))]), na.rm=TRUE)})
  poscont_AFTER = mean(unlist(dat[grep(poscont, as.character(dat[[welltype]])), grep(score.after, colnames(dat))]), na.rm=TRUE)
  negcont_AFTER = mean(unlist(dat[grep(negcont, as.character(dat[[welltype]])), grep(score.after, colnames(dat))]), na.rm=TRUE)

  qc1 = subset(dat, welltype != "Compound")[,grep(score.before, colnames(dat))]
  if (sum(qc1, na.rm=TRUE)==0)
    qc1 = NA
  else if(sum(qc1>=qc1.val, na.rm=TRUE)==0)
    qc1 = TRUE
  else
    qc1 = FALSE

  qc2 = dat[grep(poscont, as.character(dat[["welltype"]])),][,grep(score.after, colnames(dat))]
  qc2 = colMeans(qc2)

  if (sum(qc2, na.rm=TRUE)==0)
    qc2 = NA
  else if(sum(qc2<=qc2.val, na.rm=TRUE)==0)
    qc2 = TRUE
  else
    qc2 = FALSE
 
  qc3 = all(diff(c(negcont_AFTER, addcontAFTER, poscont_AFTER)) > 0)

  qc = ifelse(qc1 & qc2 & qc3, TRUE, FALSE)

  res = data.frame(passQC1=qc1, passQC2=qc2, passQC3=qc3, passQC=qc)

  return(res)
}

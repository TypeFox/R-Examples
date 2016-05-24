hits = function(dat.raw, dat.norm, score.before="scorebefore", score.after="scoreafter", qc.mainplates, qc1.val=0.225, hit.val=3){
  # selecting compounds that pass QC1
  ind = apply(dat.raw, 1, function(x){ifelse(mean(as.numeric(x[grep(score.before, names(x))]), na.rm=TRUE)>qc1.val, x[["ID"]], NA)})
  ind = ind[!is.na(ind)]
  dataux = dat.norm[!is.element(as.character(dat.norm[["ID"]]), ind),]

  IND2 = apply(dataux, 1, function(x){ifelse(mean(as.numeric(x[grep(score.after, names(x))]), na.rm=TRUE)>mean(as.numeric(x[grep(score.before, names(x))]), na.rm=TRUE), TRUE, FALSE)})
  IND3 = apply(dataux, 1, function(x){ifelse(mean(as.numeric(x[grep(score.after, names(x))]), na.rm=TRUE)>hit.val, TRUE, FALSE)})

  dataux = data.frame(dataux, IND2, IND3)
  dataux = dataux[as.character(dataux[["welltype"]])=="Compound",]

  dataux = dataux[is.element(as.character(dataux[["MainPlate"]]), as.character(qc.mainplates)),]
  dataux = subset(dataux, IND2 & IND3)

  return(dataux)
}

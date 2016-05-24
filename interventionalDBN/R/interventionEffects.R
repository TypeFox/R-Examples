interventionEffects <-
function(d,cellLine,baseline,inhibited) {
  # data is a samples x (4 + nodes) matrix or dataframe
  # column 1 gives the cell line in each sample
  # column 2 gives the inhibitor used in each sample
  # column 3 gives the stimuli used in each sample
  # column 4 gives the time index for each sample
  # cellLine is the cell line to examine
  # baseline is the entry in column 2 for the baseline (uninhibited samples)
  # inhibited is the entry in column 2 for the samples in which the inhibitor was active
  if (length(dim(d))!=2 | dim(d)[2]<5) {stop("d must be a matrix with columns 1 to 4 filled with sample information (cell line, inhibitor, stimulus, timepoint)\n")}
  dm<-data.matrix(d[,5:(dim(d)[2])])
  n.nodes<-dim(dm)[2]
  stimuli<-levels(factor(d[,3]))
  n.stimuli<-length(stimuli)
  n.tps<-max(d[,4])+1
  n.baseline<-rep(0,n.stimuli)
  n.inhibited<-rep(0,n.stimuli)
  n.baseline.used<-rep(0,n.stimuli)
  n.inhibited.used<-rep(0,n.stimuli)
  n.differences<-rep(0,n.stimuli)
  degrees.freedom<-rep(0,n.stimuli)
  names(degrees.freedom)<-stimuli
  t.statistics<-matrix(NA,n.stimuli,n.nodes)
  colnames(t.statistics)<-colnames(dm)
  rownames(t.statistics)<-stimuli
  p.values<-t.statistics
  heatmap.p.values<-t.statistics
  all.stim.t.statistics<-rep(NA,n.nodes)
  names(all.stim.t.statistics)<-colnames(dm)
  all.stim.p.values<-all.stim.t.statistics
  all.stim.heatmap.p.values<-all.stim.t.statistics
  all.stim.degrees.freedom<-NA
  diffs<-array(NA,c(n.stimuli,n.nodes,n.tps))
  for (j in 1:n.stimuli) {
    for (k in 1:n.tps-1) {
      wh1<-which(d[,1]==cellLine & d[,2]==baseline  & d[,3]==stimuli[j] & d[,4]==k)
      wh2<-which(d[,1]==cellLine & d[,2]==inhibited & d[,3]==stimuli[j] & d[,4]==k)
      n.baseline[j]<-n.baseline[j]+length(wh1)
      n.inhibited[j]<-n.inhibited[j]+length(wh2)
      if (length(wh1)>0 & length(wh2)>0) {
        n.baseline.used[j]<-n.baseline.used[j]+length(wh1)
        n.inhibited.used[j]<-n.inhibited.used[j]+length(wh2)
        a1<-dm[wh1[1],]
        a2<-dm[wh2[1],]
        if (length(wh1)>1) {a1<-apply(dm[wh1,],2,mean)}
        if (length(wh2)>1) {a2<-apply(dm[wh2,],2,mean)}
        diffs[j,,k+1]<-(a1-a2)/(1/length(wh1)+1/length(wh2))^(1/2) # This ensures that a positive T-statistic corresponds to a reduction in protein expression (ie inhibition).
      }
    }
    wh<-which(!is.na(diffs[j,1,]))# This assumes that if node 1 is observed, so are the others.
    n.differences[j]<-length(wh)
    degrees.freedom[j]<-length(wh)-1
    t.statistics[j,]<-apply(diffs[j,,wh],1,mean)/apply(diffs[j,,wh],1,sd)*sqrt(n.differences[j])
    p.values[j,]<-(1-pt(abs(t.statistics[j,]),degrees.freedom[j]))*2
    heatmap.p.values[j,]<-sign(t.statistics[j,])*(1-p.values[j,])
  }
  all.stim.degrees.freedom<-sum(n.differences)-1
  for (i in 1:n.nodes) {
    all.stim.t.statistics[i]<-mean(c(diffs[,i,]),na.rm=TRUE)/sd(c(diffs[,i,]),na.rm=TRUE)*sqrt(sum(n.differences))
  }
  all.stim.p.values<-(1-pt(abs(all.stim.t.statistics),all.stim.degrees.freedom))*2
  all.stim.heatmap.p.values<-sign(all.stim.t.statistics)*(1-all.stim.p.values)
  #
  for (j in 1:n.stimuli) {
    cat("  Stimulus",stimuli[j],":", n.baseline.used[j],"/", n.baseline[j]," baseline observations used.\n")
    cat("  Stimulus",stimuli[j],":",n.inhibited.used[j],"/",n.inhibited[j],"inhibited observations used.\n")
  }
  return(list(n.differences=n.differences,t.statistics=t.statistics,degrees.freedom=degrees.freedom,p.values=p.values,heatmap.p.values=heatmap.p.values,
    all.stim.t.statistics=all.stim.t.statistics,all.stim.degrees.freedom=all.stim.degrees.freedom,
    all.stim.p.values=all.stim.p.values,all.stim.heatmap.p.values=all.stim.heatmap.p.values))
}

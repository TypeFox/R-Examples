score2subgroup.all <-
function(score.all,treat,treat.only=TRUE) {
  stopifnot(length(score.all)>0)
  stopifnot(length(score.all)==length(treat))
  if(treat.only) stopifnot(sum(treat==1)>0)
  
  if(treat.only) {
    my.cutoff = sort(unique(score.all[treat==1]),decreasing=TRUE)
  } else {
    my.cutoff = sort(unique(score.all),decreasing=TRUE)
  }
  
  my.subs=matrix(NA,nrow=length(treat),ncol=length(my.cutoff))
  
  for(i in 1:length(my.cutoff)) {
    my.subs[,i]=1*(score.all>=my.cutoff[i])
  }
  
  if(treat.only) {
    my.pct=colMeans(my.subs[treat==1,])*100
  } else {
    my.pct=colMeans(my.subs)*100
  }
  
  return(list(cutoff=my.cutoff,score.all=score.all,subs=my.subs,pct=my.pct))
}

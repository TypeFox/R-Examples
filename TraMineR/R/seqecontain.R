seqecontain<-function(seq, eventList,exclude=FALSE){
  if(is.subseqelist(seq))seq <- seq$subseq
  if(!is.seqelist(seq))stop("seq should be a seqelist. See help on seqecreate.")
  dict<-levels.seqelist(seq)

  elist<-factor(eventList,levels=dict)
  if(exclude)excl=as.integer(c(1))
  else excl=as.integer(c(0))
  return(.Call(TMR_tmrsequencecontainevent, seq, as.integer(elist), excl))
}
print.lensMod <- function(x, rnd.coef=2, rnd.p=3,...) {
  if(class(x) != "lensMod") {stop("x must be of class 'LensMod'")}
  size.in <- data.frame(matrix(x$'Cue Validities', nrow=nrow(x$'Cue Validities'), ncol=ncol(x$'Cue Validities')*2), row.names=rownames(x$'Cue Validities'))
  size.ex <- data.frame(matrix(x$'Cue Utilizations', nrow=nrow(x$'Cue Utilizations'), ncol=ncol(x$'Cue Utilizations')*2), row.names=rownames(x$'Cue Utilizations'))
  size.in[,seq(1,ncol(size.in),by=2)] <- round(x$'Cue Validities', rnd.coef)
  size.ex[,seq(1,ncol(size.ex),by=2)] <- round(x$'Cue Utilizations', rnd.coef)
  size.in[,seq(2,ncol(size.in),by=2)] <- round(x$'Validity p', rnd.p)
  size.ex[,seq(2,ncol(size.ex),by=2)] <- round(x$'Utilization p', rnd.p)
  names.in <- rep(colnames(x$'Cue Validities'), each=2)
  names.in[seq(2,length(names.in),by=2)] <- paste(colnames(x$'Cue Validities'), "p", sep=".")
  names.ex <- rep(colnames(x$'Cue Utilizations'), each=2)
  names.ex[seq(2,length(names.ex),by=2)] <- paste(colnames(x$'Cue Utilizations'), "p", sep=".")
  colnames(size.in) <- names.in
  colnames(size.ex) <- names.ex
  out <- list("Cue Validities"=size.in, "Cue Utilizations"=size.ex)
  out
}
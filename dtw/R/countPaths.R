###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino,gmail.com>               #
#       Istituto di Ingegneria Biomedica (ISIB-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id: countPaths.R 313 2013-06-03 17:43:08Z tonig $
#                                                             #
###############################################################


countPaths <-
function(d,debug=FALSE) {
  
  N<-d$N;
  M<-d$M;
  m<-matrix(NA,nrow=N,ncol=M);
  
  if(d$openBegin) {
    m[1,]<-  1;
  } else {
    m[1,1]<-1;
  }
  
  # Some help functions
  dir<-d$stepPattern;
  npats <- attr(dir,"npat");
  nsteps <- dim(dir)[1];          # number of individual steps (counting all patterns)
  deltas <- .mkDirDeltas(dir);     # Cache the total step for each pattern

  wf <- d$windowFunction;
  
  for (ii in 1:N) {
    for (jj in 1:M) {

      # do nothing if already computed
      if( is.finite(m[ii,jj]) ) { next }

      # set to zero if outside window. Ugly but necessary to recover
      # optional arguments such as window.size
      if( do.call(d$windowFunction,
                  c(list(iw=ii,jw=jj),
                  as.list(d$call)))==FALSE ) {
        m[ii,jj]<-0; 
        next;
      }

      np<-0;
      for (k in 1:npats) {
        ni<-ii-deltas[k,1];
        nj<-jj-deltas[k,2];

        if(ni>=1 && nj>=1 ) {
          np <- np+m[ni,nj];
        }
      }
      m[ii,jj] <- np
    }
  }

  if(debug) {
    return(m)
  } 
  
  if(d$openEnd) {
    return(sum(m[N,]));
  } else {
    return(m[N,M]);
  }
  
}

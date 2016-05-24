returnLocalMaxima <-function(spectrum){
  l=length(spectrum);
  diffs=diff(spectrum);
  diffl=diffs[1:(l-1)];
  diffr=diffs[2:l];
  locMax=unique(c(which((diffl>=0&diffr<0))+1,which((diffl>0&diffr<=0))+1));
  pkMax=spectrum[locMax];
return(list(pkMax=pkMax,locMax=locMax))
}

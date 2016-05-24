getBinInfo=function(binLab,binS) {
  w=regexpr("[0-9]+[\\.]?[0-9]*", binLab,perl=T)  
  LL=as.numeric(substring(binLab, w, w+attr(w, "match.length")-1))
  rest=substring(binLab,w+attr(w, "match.length"))
  w=regexpr("[0-9]+[\\.]?[0-9]*", rest,perl=T)  
  UL=as.numeric(substring(rest, w, w+attr(w, "match.length")-1))
  c(LL=LL,UL=UL,index=which(binLab==binS))
}

plotcont = function(dat, score="score", main, xaxis.marks=seq(0,5,0.025)){
  x = unlist(dat[grep(score, colnames(dat))])
  plot(density(x[!is.na(x)]),  axes=FALSE, main=main)
  axis(1, at=xaxis.marks)
  axis(2)
  box()
}

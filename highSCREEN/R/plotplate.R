plotplate = function(dat, score="score", main){
  x = as.matrix(dat[grep(score, colnames(dat))])
  heatmap.2(x, labRow=dat[["well"]], cexRow=0.35, cexCol=0.6, srtCol=0, density.info="none", main=main)
}

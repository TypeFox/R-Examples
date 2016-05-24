pcf.greedy.kernel.new<-function(dendat,h,leaf=round(dim(dendat)[1]/2),minobs=NULL, 
itermax = 200, type="cpp"){
if (type=="old"){
  eva<-eval.greedy(dendat,leaf)
  pa<-partition(eva)
  fs<-fs.calc.parti(pa,dendat,h)
  pcf<-eva
  pcf$value<-fs
}

if (type=="prune"){
    bt<-densplit(dendat,minobs=minobs)
    treeseq<-prune(bt)
    eva<-eval.pick(treeseq,leaf=treeseq$leafs[1])  
    pa<-partition(eva)
    pcf<-eva
    fs<-fs.calc.parti(pa,dendat,h)
    pcf$value<-fs
}

if (type=="greedy"){
   pa<-densplitter(dendat,minobs=minobs)
   pcf<-list(down=pa$down,high=pa$high,grid=pa$grid,support=pa$support,
   recs=pa$recs)
   fs<-fs.calc.parti(pa,dendat,h)
   pcf$value<-fs
}

#if (type=="cpp"){
#   pcf<-densplitter2(dendat, minobs=minobs, neld = TRUE, itermax = 500)
#   pcf$value <-fsCalcParti(pcf$recs, dendat, h)
#}

return(pcf)
}


## ---- eval=FALSE, tidy=TRUE----------------------------------------------
#  rowSums(cpm(data)>100) >= 2

## ---- eval=FALSE, tidy=TRUE----------------------------------------------
#  data(Idgf2)
#  attach(Idgf2)
#  model.fit <- timeSeq(data, group.label, gene.names, reads, exon.length,  exon.level = TRUE,  offset = T)
#  

## ---- eval=FALSE, tidy=TRUE----------------------------------------------
#  model.fit$NPDE.ratio
#  model.fit$PDE.ratio
#  

## ---- eval=FALSE, tidy=TRUE----------------------------------------------
#  timeSeq.plot(model.fit, 1)

## ---- eval=FALSE, tidy=TRUE----------------------------------------------
#  model.fit <- timeSeq(data, group.label, gene.names, reads, exon.length, exon.level = TRUE, p.values = T,  offset = T, iterations = 100)
#  detach(Idgf2)

## ---- eval=FALSE, tidy=TRUE----------------------------------------------
#  data(object_by_timeSeq)
#  timeSeq.screeplot(model.fit, "lines")

## ----  eval=FALSE, tidy=TRUE---------------------------------------------
#  timeSeq.heatmap(model.fit, 10)


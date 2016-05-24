# psych9-corr   display a correlation matrix as a tableplot


# variables for psych9 Abilities, selected from psych24=Harman74.cor
v9 <-c(1,2,4,6,7,9,10,12,13)
psych9<-Harman74.cor
psych9$cov <-psych9$cov[v9,v9]
corr <- round(100 * psych9$cov)
snames <- abbreviate(rownames(corr),8)
rownames(corr) <- colnames(corr) <- snames

# show variable names on diagonal
text.m <- matrix(NA,9,9)
diag(text.m) <- snames

assign.sets <- matrix(1,9,9)
diag(assign.sets)<-0

specs <-make.specs(shape="circle", label=1, label.size=1, cell.fill="blue", scale.max=100)

tableplot(corr, assign.sets, cell.specs=specs, text.m=text.m, 
	left.space=5, top.space=5, table.label=FALSE,
	v.part=rep(3,3), h.part=rep(3,3))

# use different color & symbols for blocks of 3
assign.sets <- 1+outer(rep(1:3, each=3), rep(1:3, each=3), "!=")
diag(assign.sets)<-0


specs <-make.specs(shape=c("circle","square"), label=1, label.size=1, 
	cell.fill=c("blue","gray"), scale.max=100)
tableplot(corr, assign.sets, cell.specs=specs, text.m=text.m, empty.text.size=0.9,
	left.space=5, top.space=5, table.label=FALSE,
	v.part=rep(3,3), h.part=rep(3,3))

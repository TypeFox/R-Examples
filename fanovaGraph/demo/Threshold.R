#################################################################
# Several methods to support threshold decision
#################################################################

require(fanovaGraph)

### example data set with two interactions
d = 4
x <- matrix(runif(40*d), 40, d)
y <- (x[,1]-0.5) * (x[,2]-0.5) + 0.9*(x[,1]-0.5) * (x[,3]-0.5) 

### kriging prediction model
KM <- km(~1, design = data.frame(x), response = y)

### graph estimation
g <- estimateGraph(kmPredictWrapper, d=d, n.tot = 10000, km.object=KM)

##################################################################
### Threshold decision

### examine full graph

plot(g, plot.i1 = FALSE)

### Compare candidate thresholds on prediction performance
comparison <- thresholdIdentification(g, x, y, n.cand = 2)

### Delta Jump Plot
plotDeltaJumps(g)
plotDeltaJumps(g, mean.clique.size=TRUE)

### see graph changing as delta chages

plotGraphChange(g, fix.layout = TRUE)

### see effect of delta on graph interactively with library tcltk
plotTk(g)

### the same with library manipulate
plotManipulate(g)

### 'unknown' true threshold
g.cut <- threshold(g, delta = 0.2, scale = TRUE)
plot(g.cut)



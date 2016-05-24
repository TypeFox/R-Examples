### R code from vignette source 'streamMOA.Rnw'

###################################################
### code chunk number 1: streamMOA.Rnw:5-6
###################################################
options(width = 75, digits = 3, prompt = 'R> ')


###################################################
### code chunk number 2: streamMOA.Rnw:31-32
###################################################
set.seed(1234) 


###################################################
### code chunk number 3: data_bng
###################################################
library("stream") 
stream <- DSD_Memory(DSD_BarsAndGaussians(noise=0.05), n=5500)
stream
plot(stream)


###################################################
### code chunk number 4: streamMOA.Rnw:57-61
###################################################
sample <- DSC_TwoStage(micro=DSC_Sample(k=100), macro=DSC_Kmeans(k=4)) 
window <- DSC_TwoStage(micro=DSC_Window(horizon=100), macro=DSC_Kmeans(k=4))
dstream <- DSC_DStream(gridsize=.7)
dbstream <- DSC_DBSTREAM(r=.5)


###################################################
### code chunk number 5: streamMOA.Rnw:66-69
###################################################
library("streamMOA")
denstream <- DSC_DenStream(epsilon=.5, mu=1) 
clustream <- DSC_CluStream(m=100, k=4) 


###################################################
### code chunk number 6: streamMOA.Rnw:76-82
###################################################
algorithms <- list(Sample=sample, Window=window, 'D-Stream'=dstream, 
  DBSTREAM=dbstream, DenStream=denstream, CluStream=clustream)
for(a in algorithms) {
  reset_stream(stream) 
  update(a, stream, 5000)
}


###################################################
### code chunk number 7: streamMOA.Rnw:86-87
###################################################
sapply(algorithms, nclusters, type="micro")


###################################################
### code chunk number 8: microclusters
###################################################
op <- par(no.readonly = TRUE)
layout(mat=matrix(1:6, ncol=2))
for(a in algorithms) {
  reset_stream(stream) 
  plot(a, stream, main=description(a), type="micro")
}
par(op)


###################################################
### code chunk number 9: microclusters_assignment
###################################################
op <- par(no.readonly = TRUE)
layout(mat=matrix(1:6, ncol=2))
for(a in algorithms) {
  reset_stream(stream) 
  plot(a, stream, main=description(a), assignment=TRUE, weight=FALSE, type="micro")
}
par(op)


###################################################
### code chunk number 10: streamMOA.Rnw:169-175
###################################################
sapply(algorithms, FUN=function(a) {
  reset_stream(stream, 1001) 
  evaluate(a, stream, 
    measure=c("numMicroClusters", "purity", "SSQ", "silhouette"), 
    n=500, assignmentMethod="auto", type="micro")
})


###################################################
### code chunk number 11: macroclusters
###################################################
op <- par(no.readonly = TRUE)
layout(mat=matrix(1:6, ncol=2))
for(a in algorithms) {
  reset_stream(stream) 
  plot(a, stream, main=description(a), type="both")
}
par(op)


###################################################
### code chunk number 12: streamMOA.Rnw:223-228
###################################################
sapply(algorithms, FUN=function(a) {
  reset_stream(stream, 1001) 
  evaluate(a, stream, measure=c("numMacroClusters","purity", "SSQ", "cRand"), 
    n=500, assign="micro", type="macro")
})


###################################################
### code chunk number 13: streamMOA.Rnw:246-248
###################################################
set.seed(0)
stream <- DSD_Memory(DSD_Benchmark(1), 5000)


###################################################
### code chunk number 14: moa1
###################################################
plot(stream, 250, xlim=c(0,1), ylim=c(0,1))
arrows(.15,.85,.85,.15, col=rgb(.8,.8,.8,.6), lwd=10)
arrows(.15,.15,.85,.85, col=rgb(.8,.8,.8,.6), lwd=10)
reset_stream(stream)


###################################################
### code chunk number 15: streamMOA.Rnw:273-284
###################################################
algorithms <- list(
  Sample = DSC_TwoStage(micro=DSC_Sample(k=100, biased=TRUE), 
    macro=DSC_Kmeans(k=2)), 
  Window = DSC_TwoStage(micro=DSC_Window(horizon=100, lambda=.01), 
    macro=DSC_Kmeans(k=2)),
  
  'D-Stream' = DSC_DStream(gridsize=.05, lambda=.01),
  DBSTREAM = DSC_DBSTREAM(r=.02, lambda=.01),
  DenStream = DSC_DenStream(epsilon=.05, lambda=.01),
  CluStream = DSC_CluStream(m=100, k=2) 
)


###################################################
### code chunk number 16: streamMOA.Rnw:294-305
###################################################
n <- 5000
horizon <- 250
reset_stream(stream)

evaluation <- lapply(algorithms, FUN=function(a) {
  reset_stream(stream)
  evaluate_cluster(a, stream, 
    type="macro", assign="micro",
    measure=c("numMicro","numMacro","SSQ", "cRand"), 
    n=n, horizon=horizon)
})


###################################################
### code chunk number 17: dynamic
###################################################
Position <- evaluation[[1]][,"points"]
cRand <- sapply(evaluation, FUN=function(x) x[,"cRand"])
cRand
matplot(Position, cRand, type="l", lwd=2)
legend("bottomleft", legend=names(evaluation), 
  col=1:6, lty=1:6, bty="n", lwd=2)


###################################################
### code chunk number 18: dynamic_box
###################################################
boxplot(cRand, las=2)


###################################################
### code chunk number 19: dynamic2
###################################################
SSQ <- sapply(evaluation, FUN=function(x) x[,"SSQ"])
SSQ
matplot(Position, SSQ, type="l", lwd=2)
legend("topright", legend=names(evaluation), 
  col=1:6, lty=1:6, bty="n", lwd=2)


###################################################
### code chunk number 20: dynamic_box2
###################################################
boxplot(SSQ, las=2)



### R code from vignette source 'stream.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: stream.Rnw:136-137
###################################################
options(width = 75, digits = 3, prompt = 'R> ', scipen = 3)


###################################################
### code chunk number 2: stream.Rnw:652-655
###################################################
library("stream")
set.seed(1000)
stream <- DSD_Gaussians(k = 3, d = 2)


###################################################
### code chunk number 3: stream.Rnw:664-666
###################################################
dstream <- DSC_DStream(gridsize = .1, Cm = 1.2)
update(dstream, stream, n = 500)


###################################################
### code chunk number 4: initial_example
###################################################
km <- DSC_Kmeans(k = 3)
recluster(km, dstream)
plot(km, stream, type = "both")


###################################################
### code chunk number 5: stream.Rnw:882-887
###################################################
library("stream")
set.seed(1000) 

stream <- DSD_Gaussians(k = 3, d = 3, noise = .05, p = c(.5, .3, .1)) 
stream


###################################################
### code chunk number 6: stream.Rnw:908-910
###################################################
p <- get_points(stream, n = 5) 
p 


###################################################
### code chunk number 7: stream.Rnw:921-923
###################################################
p <- get_points(stream, n = 100, class = TRUE) 
head(p, n = 10)


###################################################
### code chunk number 8: static
###################################################
plot(stream, n = 500)


###################################################
### code chunk number 9: static_pc
###################################################
plot(stream, n = 500, method = "pc")


###################################################
### code chunk number 10: moa1
###################################################
set.seed(1000) 
stream <- DSD_Benchmark(1)
stream


###################################################
### code chunk number 11: stream.Rnw:983-987 (eval = FALSE)
###################################################
## for(i in 1:4) {
##   plot(stream, 250, xlim = c(0, 1), ylim = c(0, 1))
##   tmp <- get_points(stream, n = 1400)
## }


###################################################
### code chunk number 12: moa1
###################################################
plot(stream, 250, xlim = c(0, 1), ylim = c(0, 1))
arrows(.15, .85, .85, .15, col = rgb(.8, .8, .8, .6), lwd = 10)
arrows(.15, .15, .85, .85, col = rgb(.8, .8, .8, .6), lwd = 10)
tmp <- get_points(stream, n = 1400)


###################################################
### code chunk number 13: moa2
###################################################
plot(stream, 250, xlim = c(0, 1), ylim = c(0, 1))
arrows(.15, .85, .85, .15, col = rgb(.8, .8, .8, .6), lwd = 10)
arrows(.15, .15, .85, .85, col = rgb(.8, .8, .8, .6), lwd = 10)
tmp <- get_points(stream, n=1400)


###################################################
### code chunk number 14: moa3
###################################################
plot(stream, 250, xlim = c(0, 1), ylim = c(0, 1))
arrows(.15,.85,.85,.15, col=rgb(.8,.8,.8,.6), lwd=10)
arrows(.15,.15,.85,.85, col=rgb(.8,.8,.8,.6), lwd=10)
tmp <- get_points(stream, n=1400)


###################################################
### code chunk number 15: moa4
###################################################
plot(stream, 250, xlim=c(0,1), ylim=c(0,1))
arrows(.15,.85,.85,.15, col=rgb(.8,.8,.8,.6), lwd=10)
arrows(.15,.15,.85,.85, col=rgb(.8,.8,.8,.6), lwd=10)


###################################################
### code chunk number 16: stream.Rnw:1039-1042 (eval = FALSE)
###################################################
## reset_stream(stream)
## animate_data(stream, n = 10000, horizon = 100, 
##   xlim = c(0, 1), ylim = c(0, 1))


###################################################
### code chunk number 17: stream.Rnw:1048-1051 (eval = FALSE)
###################################################
## library("animation")
## animation::ani.options(interval = .1)
## ani.replay()


###################################################
### code chunk number 18: stream.Rnw:1058-1060 (eval = FALSE)
###################################################
## saveHTML(ani.replay())
## saveGIF(ani.replay())


###################################################
### code chunk number 19: stream.Rnw:1084-1087
###################################################
library("stream") 
set.seed(1000) 
stream <- DSD_Gaussians(k = 3, d = 5) 


###################################################
### code chunk number 20: stream.Rnw:1092-1093 (eval = FALSE)
###################################################
## write_stream(stream, "data.csv", n = 100, sep = ",") 


###################################################
### code chunk number 21: stream.Rnw:1123-1127
###################################################
file <- system.file("examples", "kddcup10000.data.gz", package = "stream")
stream_file <- DSD_ReadCSV(gzfile(file),
  take = c(1, 5, 6, 8:11, 13:20, 23:42), class = 42, k = 7)
stream_file


###################################################
### code chunk number 22: stream.Rnw:1139-1140
###################################################
get_points(stream_file, n = 5)


###################################################
### code chunk number 23: stream.Rnw:1147-1149
###################################################
stream_scaled <- DSD_ScaleStream(stream_file, center = TRUE, scale = TRUE)
get_points(stream_scaled, n = 5)


###################################################
### code chunk number 24: stream.Rnw:1180-1182
###################################################
data("EuStockMarkets", package = "datasets")
head(EuStockMarkets)


###################################################
### code chunk number 25: stream.Rnw:1189-1191
###################################################
replayer <- DSD_Memory(EuStockMarkets, k = NA) 
replayer 


###################################################
### code chunk number 26: stream.Rnw:1197-1199
###################################################
get_points(replayer, n = 5)
replayer


###################################################
### code chunk number 27: stream.Rnw:1207-1208 (eval = FALSE)
###################################################
## get_points(replayer, n = 2000)


###################################################
### code chunk number 28: stream.Rnw:1210-1212
###################################################
err <- try(get_points(replayer, n = 2000))
cat(err)


###################################################
### code chunk number 29: stream.Rnw:1226-1228
###################################################
reset_stream(replayer, pos = 100)
replayer


###################################################
### code chunk number 30: stream.Rnw:1534-1537
###################################################
library("stream") 
set.seed(1000) 
stream <- DSD_Gaussians(k = 3, d = 2, noise = .05)


###################################################
### code chunk number 31: stream.Rnw:1545-1547
###################################################
dstream <- DSC_DStream(gridsize = .1, Cm = 1.2) 
dstream


###################################################
### code chunk number 32: stream.Rnw:1555-1557
###################################################
update(dstream, stream, n = 500) 
dstream


###################################################
### code chunk number 33: stream.Rnw:1566-1567
###################################################
head(get_centers(dstream))


###################################################
### code chunk number 34: cluster
###################################################
plot(dstream, stream)


###################################################
### code chunk number 35: cluster-grid
###################################################
plot(dstream, stream, grid = TRUE)


###################################################
### code chunk number 36: stream.Rnw:1887-1891
###################################################
library("stream") 
stream <- DSD_Gaussians(k = 3, d = 2, noise = .05)
dstream <- DSC_DStream(gridsize = .1) 
update(dstream, stream, n = 2000) 


###################################################
### code chunk number 37: stream.Rnw:1899-1900
###################################################
evaluate(dstream, stream, n = 100)


###################################################
### code chunk number 38: stream.Rnw:1911-1912
###################################################
evaluate(dstream, stream, measure = c("purity", "crand"), n = 500)


###################################################
### code chunk number 39: stream.Rnw:1953-1959
###################################################
set.seed(1000)
stream <- DSD_Benchmark(1)
dstream <- DSC_DStream(gridsize = .05, lambda = .01)
ev <- evaluate_cluster(dstream, stream, 
  measure = c("numMicroClusters", "purity"), n = 5000, horizon = 100)
head(ev)


###################################################
### code chunk number 40: evaluation
###################################################
plot(ev[ , "points"], ev[ , "purity"], type = "l", 
  ylab = "Avg. Purity", xlab = "Points")


###################################################
### code chunk number 41: stream.Rnw:1992-1997 (eval = FALSE)
###################################################
## set.seed(1000)
## stream <- DSD_Benchmark(1)
## dstream <- DSC_DStream(gridsize = .05, lambda = .01)
## r <- animate_cluster(dstream, stream, horizon = 100, n = 5000, 
##      measure = "purity", plot.args = list(xlim = c(0, 1), ylim = c(0, 1)))


###################################################
### code chunk number 42: stream.Rnw:2027-2034
###################################################
library("stream")
set.seed(1000) 
stream <- DSD_Gaussians(k = 3, d = 2, noise = .05)
dstream <- DSC_DStream(gridsize = .05, Cm = 1.5)

update(dstream, stream, n = 1000)
dstream


###################################################
### code chunk number 43: recluster
###################################################
plot(dstream, stream, type = "both")


###################################################
### code chunk number 44: recluster2
###################################################
km <- DSC_Kmeans(k = 3, weighted = TRUE)
recluster(km, dstream)
km
plot(km, stream, type = "both") 


###################################################
### code chunk number 45: stream.Rnw:2095-2096
###################################################
evaluate(km, stream, measure = c("purity", "crand", "SSQ"), n = 1000)


###################################################
### code chunk number 46: stream.Rnw:2101-2103
###################################################
evaluate(km, stream, c(measure = "purity", "crand", "SSQ"), n = 1000, 
  assign = "macro")


###################################################
### code chunk number 47: stream.Rnw:2126-2129
###################################################
points <- get_points(stream, n = 100)
assignment <- get_assignment(dstream, points, type = "macro")
assignment


###################################################
### code chunk number 48: silhouette
###################################################
assignment[is.na(assignment)] <- 0L
library("cluster")
plot(silhouette(assignment, dist = dist(points)))


###################################################
### code chunk number 49: data_bng
###################################################
set.seed(1000) 
library("stream") 
stream <- DSD_Memory(DSD_BarsAndGaussians(noise = .05), n = 1500)
stream
plot(stream)


###################################################
### code chunk number 50: stream.Rnw:2354-2362
###################################################
algorithms <- list(
  'Sample' = DSC_TwoStage(micro = DSC_Sample(k = 100), 
    macro = DSC_Kmeans(k = 4)), 
  'Window' = DSC_TwoStage(micro = DSC_Window(horizon = 100), 
    macro = DSC_Kmeans(k = 4)), 
  'D-Stream' = DSC_DStream(gridsize = .7, Cm = 1.5),
  'DBSTREAM' = DSC_DBSTREAM(r = .45)
)


###################################################
### code chunk number 51: stream.Rnw:2372-2376
###################################################
for(a in algorithms) {
  reset_stream(stream) 
  update(a, stream, n = 1000)
}


###################################################
### code chunk number 52: stream.Rnw:2381-2382
###################################################
sapply(algorithms, nclusters, type = "micro")


###################################################
### code chunk number 53: microclusters
###################################################
op <- par(no.readonly = TRUE)
layout(mat = matrix(1:length(algorithms), ncol = 2))
for(a in algorithms) {
  reset_stream(stream) 
  plot(a, stream, main = description(a), type = "micro")
}
par(op)


###################################################
### code chunk number 54: microclusters_assignment
###################################################
op <- par(no.readonly = TRUE)
layout(mat = matrix(1:length(algorithms), ncol = 2))
for(a in algorithms) {
  reset_stream(stream) 
  plot(a, stream, main = description(a), 
    assignment = TRUE, weight = FALSE, type = "micro")
}
par(op)


###################################################
### code chunk number 55: stream.Rnw:2460-2467
###################################################
sapply(algorithms, FUN=function(a) {
  reset_stream(stream, pos = 1001) 
  evaluate(a, stream, 
    measure = c("numMicroClusters", "purity"),
    type = "micro",
    n = 500)
})


###################################################
### code chunk number 56: macroclusters
###################################################
op <- par(no.readonly = TRUE)
layout(mat=matrix(1:length(algorithms), ncol = 2))
for(a in algorithms) {
  reset_stream(stream) 
  plot(a, stream, main = description(a), type = "both")
}
par(op)


###################################################
### code chunk number 57: stream.Rnw:2511-2517
###################################################
sapply(algorithms, FUN = function(a) {
  reset_stream(stream, pos = 1001) 
  evaluate(a, stream, measure = c("numMacroClusters", "purity", 
      "SSQ", "cRand", "silhouette"), 
    n = 500, assign = "micro", type = "macro")
})


###################################################
### code chunk number 58: stream.Rnw:2539-2541
###################################################
set.seed(0)
stream <- DSD_Memory(DSD_Benchmark(1), n = 5000)


###################################################
### code chunk number 59: stream.Rnw:2551-2559
###################################################
algorithms <- list(
  'Sample' = DSC_TwoStage(micro = DSC_Sample(k = 100, biased = TRUE), 
    macro = DSC_Kmeans(k = 2)), 
  'Window' = DSC_TwoStage(micro = DSC_Window(horizon = 100, lambda = .01), 
    macro = DSC_Kmeans(k = 2)), 
  'D-Stream' = DSC_DStream(gridsize = .1, lambda = .01),
  'DBSTREAM' = DSC_DBSTREAM(r = .05, lambda = .01)
)


###################################################
### code chunk number 60: stream.Rnw:2570-2575
###################################################
evaluation <- lapply(algorithms, FUN = function(a) {
  reset_stream(stream) 
  evaluate_cluster(a, stream, horizon = 100, n = 5000, measure = "crand",
    type = "macro", assign = "micro")
})


###################################################
### code chunk number 61: stream.Rnw:2591-2594
###################################################
Position <- evaluation[[1]][ , "points"]
cRand <- sapply(evaluation, FUN = function(x) x[ , "cRand"])
head(cRand)


###################################################
### code chunk number 62: dynamic
###################################################
matplot(Position, cRand, type = "l", lwd = 2)
legend("bottomleft", legend = names(evaluation), 
  col = 1:6, lty = 1:6, bty = "n", lwd = 2)


###################################################
### code chunk number 63: dynamic_box
###################################################
boxplot(cRand, las = 2, cex.axis = .8)


###################################################
### code chunk number 64: stream.Rnw:2652-2660 (eval = FALSE)
###################################################
## library("stream")
## con <- gzcon(
##   url(paste("http://archive.ics.uci.edu/ml/machine-learning-databases/",
##     "kddcup99-mld/kddcup.data.gz", sep="")))
## 
## stream <- DSD_ReadCSV(con, take=c(1, 5, 6, 8:11, 13:20, 23:42), 
##     class=42, k=7)
## stream2 <- DSD_ScaleStream(stream, n=1000)


###################################################
### code chunk number 65: stream.Rnw:2666-2668 (eval = FALSE)
###################################################
## dstream <- DSC_DStream(gridsize = .5, gaptime = 10000L, lambda = .01)
## update(dstream, stream2, n = 4000000, verbose = TRUE)



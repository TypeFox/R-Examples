## ----setup, include=FALSE, cache=FALSE-----------------------------------
library(knitr)
# set global chunk options
options(formatR.arrow=FALSE)
opts_chunk$set(fig.path='figure/Rplots-',fig.align='center',fig.show='hold',comment=NA,background='white',highlight=FALSE,tidy=TRUE,size="small",continue=" ")
knit_hooks$set(source=function(x,options){
    prp <- c("R> ")
    if(!options$prompt) prp <- ""
    wd <- getOption("width")
    if(!is.null(width <- options$tidy.opts$width))
        options(width = width)
    x <- strwrap(x, width = getOption("width"))
    lenx <- length(x)
    pl <- unlist(sapply(gregexpr("\\(", x), function(el){
        if((length(el) == 1))
            if(unique(el) == -1) 0 else 1
        else length(el)}))
    pr <- unlist(sapply(gregexpr("\\)", x), function(el){
        if((length(el) == 1))
            if(unique(el) == -1) 0 else 1
        else length(el)}))
    wp <- rep(prp, length(x))
    if(length(x) > 1){
        xns <- gsub(" ","",x)
        op <- gregexpr("\\+|-|\\*|\\|=",x)
        ct <- sapply(1:(length(x) - 1), function(i, xns, op)
                     (nchar(x[i]) %in% op[[i]]) | (1 %in% op[[i + 1]]), xns,op)
        for(i in 2:length(x)){
            if((sum(pl[1:(i-1)]) != sum(pr[1:(i-1)])) | ct[i - 1])
                wp[i] <- paste(options$continue, "   ", sep = "")
        }
    }
    options(width = wd)
    paste(c("\\begin{Rinput}",paste(wp, x, sep= ""), "\\end{Rinput}",""), collapse = "\n")
},  output=function(x,options){
     if(all(gregexpr("begin\\{tabu|begin\\{longtab",x)[[1]] > 0)) x
     else paste(c("\\begin{Routput}\n",x, "\\end{Routput}\n"), sep = "")
 })

## ----setup-lib, include=FALSE, cache=FALSE-------------------------------
library(ASMap)

## ----data0, eval = TRUE, echo = TRUE,prompt = TRUE-----------------------
data(mapDHf, package = "ASMap")
data(mapDH, package = "ASMap")
data(mapBCu, package = "ASMap")

## ----mst-df,eval=FALSE,echo=TRUE,prompt=FALSE----------------------------
#  mstmap.data.frame(object, pop.type = "DH", dist.fun = "kosambi",
#        objective.fun = "COUNT", p.value = 1e-06, noMap.dist = 15,
#        noMap.size = 0, miss.thresh = 1, mvest.bc = FALSE, detectBadData = FALSE,
#        as.cross = TRUE, return.imputed = TRUE, trace = FALSE, ...)

## ----datadf,eval=TRUE,echo=TRUE,prompt=TRUE------------------------------
testd <- mstmap(mapDHf, dist.fun = "kosambi", trace = TRUE, as.cross = TRUE)
nmar(testd)
chrlen(testd)

## ----mst-cr,eval = FALSE,echo=TRUE---------------------------------------
#  mstmap.cross(object, chr, id = "Genotype", bychr = TRUE,
#         suffix = "numeric", anchor = FALSE, dist.fun = "kosambi",
#         objective.fun = "COUNT", p.value = 1e-06, noMap.dist = 15,
#         noMap.size = 0, miss.thresh = 1, mvest.bc = FALSE, detectBadData =
#         FALSE, return.imputed = FALSE, trace = FALSE, ...)

## ----data, eval = TRUE, echo = TRUE, prompt = TRUE-----------------------
nmar(mapDH)
pull.map(mapDH)[[4]]

## ----mst1,eval=TRUE,echo=TRUE,prompt=TRUE--------------------------------
mapDHa <- mstmap(mapDH, bychr = FALSE, dist.fun = "kosambi", trace = TRUE)
nmar(mapDHa)
pull.map(mapDHa)[[4]]

## ----mst2,eval=TRUE,echo=TRUE,prompt=TRUE--------------------------------
mapDHb <- mstmap(mapDH, bychr = TRUE, dist.fun = "kosambi", anchor = TRUE, trace = TRUE)
nmar(mapDHb)

## ----mst3,eval=TRUE,echo=TRUE,prompt=TRUE--------------------------------
mapDHc <- mstmap(mapDH, bychr = TRUE, dist.fun = "kosambi", anchor = TRUE, trace = TRUE, p.value = 1e-04)
nmar(mapDHc)

## ----mst4,eval=TRUE,echo=TRUE,prompt=TRUE--------------------------------
mapDHd <- mstmap(mapDH, chr = names(mapDH$geno)[1:3], bychr = FALSE, dist.fun = "kosambi", trace = TRUE, p.value = 1e-04)
nmar(mapDHd)

## ----pp1,eval=FALSE,echo=TRUE--------------------------------------------
#  pullCross(object, chr, type = c("co.located","seg.distortion","missing"),
#                 pars = NULL, replace = FALSE, ...)
#  pushCross(object, chr, type = c("co.located","seg.distortion","missing","unlinked"),
#                 unlinked.chr = NULL, pars = NULL, replace = FALSE, ...)
#  pp.init(seg.thresh = 0.05, seg.ratio = NULL, miss.thresh = 0.1, max.rf =
#               0.25, min.lod = 3)

## ----pp2,eval=TRUE,echo=TRUE,prompt=TRUE---------------------------------
mapDHs <- pullCross(mapDH, type = "co.located")
mapDHs <- pullCross(mapDHs, type = "seg.distortion", pars = list(seg.thresh = 0.02))
mapDHs <- pullCross(mapDHs, type = "missing", pars = list(miss.thresh = 0.03))
names(mapDHs)
names(mapDHs$co.located)

## ----pp3,eval=TRUE,echo=TRUE,prompt=TRUE---------------------------------
mapDHs$seg.distortion$table

## ----pp4,eval=TRUE,echo=TRUE,prompt=TRUE---------------------------------
head(mapDHs$co.located$table)

## ----pp5,eval=TRUE,echo=TRUE,prompt=TRUE---------------------------------
mapDHs <- mstmap(mapDHs, bychr = FALSE, dist.fun = "kosambi", trace = TRUE, anchor = TRUE)
nmar(mapDHs)

## ----pp6,eval=TRUE,echo=TRUE,prompt=TRUE---------------------------------
mapDHs <- pushCross(mapDHs, type = "co.located")
mapDHs <- pushCross(mapDHs, type = "seg.distortion", pars = list(seg.thresh = 0.001))
mapDHs <- pushCross(mapDHs, type = "missing", pars = list(miss.thresh = 0.05))
names(mapDHs)

## ----pp7,eval=TRUE,echo=TRUE,prompt=TRUE---------------------------------
pull.map(mapDHs)[[4]]
pull.map(mapDHs)[[21]]

## ----pp8,eval=TRUE,echo=TRUE,prompt=TRUE---------------------------------
mapDHs <- mstmap(mapDHs, bychr = TRUE, dist.fun = "kosambi", trace = TRUE, anchor = TRUE, p.value = 2)

## ----heat1,eval=FALSE,echo=TRUE------------------------------------------
#  heatMap(x, chr, mark, what = c("both", "lod", "rf"), lmax = 12,
#               rmin = 0, markDiagonal = FALSE, color = rev(rainbow(256, start =
#               0, end = 2/3)), ...)

## ----heat2,echo=TRUE,eval=FALSE,prompt=TRUE------------------------------
#  heatMap(mapDH, lmax = 50)

## ----prof1, eval = FALSE-------------------------------------------------
#  statGen(cross, chr, bychr = TRUE, stat.type = c("xo", "dxo", "miss"), id = "Genotype")
#  profileGen(cross, chr, bychr = TRUE, stat.type = c("xo", "dxo", "miss"), id = "Genotype", xo.lambda = NULL, ...)
#  statMark(cross, chr, stat.type = c("marker", "interval"), map.function = "kosambi")
#  profileMark(cross, chr, stat.type = "marker", use.dist = TRUE, map.function = "kosambi", crit.val = NULL, display.markers = FALSE, mark.line = FALSE, ...)

## ----prof2,fig.width = 15,fig.height = 8,fig.pos = "t",fig.env="figure",fig.scap="NA",fig.cap = "Genotype profiles of missing values, double recombinations and recombinations for \\texttt{mapDH}.",prompt=TRUE----
profileGen(mapDH, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id = "Genotype", xo.lambda = 25, layout = c(1,3), lty = 2)

## ----prof3,fig.width = 15,fig.height = 8,fig.pos = "t",fig.env="figure",fig.scap="NA",fig.cap = "Marker and interval profiles of segregation distortion, double crossovers, estimated recombination fractions and LOD scores for \\texttt{mapDH}.",prompt=TRUE, warning = FALSE----
profileMark(mapDH, stat.type = c("seg.dist", "dxo", "erf", "lod"), id = "Genotype", layout = c(1,4), type = "l")

## ----clones01,eval=FALSE,echo=TRUE,prompt=TRUE---------------------------
#  genClones(object, chr, tol = 0.9, id = "Genotype")
#  fixClones(object, gc, id = "Genotype", consensus = TRUE)

## ----clones02,eval=TRUE,echo=TRUE,prompt=TRUE----------------------------
gc <- genClones(mapDH, tol = 0.9)
gc$cgd

## ----clones03,eval=TRUE,echo=TRUE,prompt=TRUE----------------------------
mapDHg <- fixClones(mapDH, gc$cgd, consensus = TRUE)
levels(mapDHg$pheno[[1]])[grep("_", levels(mapDHg$pheno[[1]]))]

## ----mb1, eval = FALSE---------------------------------------------------
#  breakCross(cross, split = NULL, suffix = "numeric", sep = ".")
#  mergeCross(cross, merge = NULL, gap = 5)

## ----mb2, eval = TRUE,prompt=TRUE----------------------------------------
mapDHb1 <- breakCross(mapDH, split = list("3B" = "3B.m.7","6A" = "6A.m.15"))
nmar(mapDHb1)

## ----mb3, eval = TRUE,prompt=TRUE----------------------------------------
mapDHb2 <- breakCross(mapDH, split = list("3B" = "3B.m.7"), suffix = list("3B" = c("3B1","3B2")))
nmar(mapDHb2)

## ----mb4,eval=TRUE,prompt=TRUE-------------------------------------------
mapDHm <- mergeCross(mapDHb1, merge = list("3B" = c("3B.1","3B.2"),"6A" = c("6A.1","6A.2")))
nmar(mapDHm)

## ----quick1,eval=FALSE---------------------------------------------------
#  quickEst(object, chr, map.function = "kosambi", ...)

## ----quick2,fig.width = 7,fig.height = 5,fig.pos = "t",fig.env="figure",fig.scap="NA",fig.cap = "Comparison of \\texttt{mapDH} using \\texttt{est.map} and \\texttt{quickEst}.",prompt=TRUE----
map1 <- est.map(mapDH, map.function = "kosambi")
map1 <- subset(map1, chr = names(nmar(map1))[6:15])
map2 <- quickEst(mapDH, map.function = "kosambi")
map2 <- subset(map2, chr = names(nmar(map2))[6:15])
plot.map(map1, map2)

## ----sc,eval=TRUE,prompt=TRUE--------------------------------------------
mapDH.s <- pullCross(mapDH, type = "seg.distortion")
mapDH.s <- subsetCross(mapDH.s, ind = 3:218)
dim(mapDH.s$seg.distortion$data)[1]

## ----comb1, eval = FALSE-------------------------------------------------
#  combineMap(..., id = "Genotype", keep.all = TRUE)

## ----comb2,eval=TRUE,prompt=TRUE-----------------------------------------
mapDH1 <- mapDH
names(mapDH1$geno)[5:14] <- paste("L",1:10, sep = "")
mapDH1$geno <- lapply(mapDH1$geno, function(el){
  names(el$map) <- dimnames(el$data)[[2]] <- paste(names(el$map), "A", sep = "")
  el})
mapDHc <- combineMap(mapDH, mapDH1)
nmar(mapDHc)

## ----ex1, eval = FALSE---------------------------------------------------
#  data(mapBCu, package = "ASMap")

## ----ex3,eval=FALSE,echo=TRUE,prompt=TRUE--------------------------------
#  plot.missing(mapBCu)

## ----ex4,echo=TRUE,prompt=TRUE-------------------------------------------
sg <- statGen(mapBCu, bychr = FALSE, stat.type = "miss")
mapBC1 <- subset(mapBCu, ind = sg$miss < 1600)

## ----ex5,eval=TRUE,echo=TRUE,prompt=TRUE---------------------------------
gc <- genClones(mapBC1, tol = 0.95)
gc$cgd

## ----ex6,eval=TRUE,echo=TRUE,prompt=TRUE---------------------------------
cgd <- gc$cgd[-c(1,4,5),]
mapBC2 <- fixClones(mapBC1, cgd, consensus = TRUE)
levels(mapBC2$pheno[[1]])[grep("_", levels(mapBC2$pheno[[1]]))]

## ----ex7,eval=FALSE,echo=TRUE,prompt=TRUE--------------------------------
#  profileMark(mapBC2, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1,4), type = "l", cex = 0.5)

## ----ex8,echo=FALSE,fig.width=17,fig.height=10,warning=FALSE-------------
profileMark(mapBC2, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1,4), type = "l", cex = 0.5)

## ----ex9,eval=TRUE,echo=TRUE,prompt=TRUE---------------------------------
mm <- statMark(mapBC2, stat.type = "marker")$marker$AB
mapBC3 <- drop.markers(mapBC2, c(markernames(mapBC2)[mm > 0.98],markernames(mapBC2)[mm < 0.2]))

## ----ex10,eval=TRUE,echo=TRUE,prompt=TRUE--------------------------------
mapBC3 <- pullCross(mapBC3, type = "missing", pars = list(miss.thresh = 0.1))
mapBC3 <- pullCross(mapBC3, type = "seg.distortion", pars = list(seg.thresh = "bonf"))
mapBC3 <- pullCross(mapBC3, type = "co.located")
names(mapBC3)
sum(ncol(mapBC3$missing$data),ncol(mapBC3$seg.dist$data),ncol(mapBC3$co.located$data))

## ----ex11,eval=TRUE,echo=TRUE,prompt=TRUE,cache=TRUE---------------------
mapBC4 <- mstmap(mapBC3, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-12)
chrlen(mapBC4)

## ----ex12,eval=FALSE,echo=TRUE,prompt=TRUE-------------------------------
#  heatMap(mapBC4, lmax = 70)

## ----ex14,eval=FALSE,echo=TRUE,prompt=TRUE-------------------------------
#  pg <- profileGen(mapBC4, bychr = FALSE, stat.type = c("xo","dxo","miss"), id = "Genotype", xo.lambda = 14, layout = c(1,3), lty = 2, cex = 0.7)

## ----ex15,echo=FALSE,fig.width=17,fig.height=10,warning=FALSE------------
pg <- profileGen(mapBC4, bychr = FALSE, stat.type = c("xo","dxo","miss"), id = "Genotype", xo.lambda = 14, layout = c(1,3), lty = 2, cex = 0.7)

## ----ex16,eval=TRUE,echo=TRUE,cache = TRUE,prompt=TRUE-------------------
mapBC5 <- subsetCross(mapBC4, ind = !pg$xo.lambda)
mapBC6 <- mstmap(mapBC5, bychr = TRUE, dist.fun = "kosambi", trace = TRUE, p.value = 1e-12)
chrlen(mapBC6)

## ----ex17,eval=FALSE,echo=TRUE,prompt=TRUE-------------------------------
#  profileMark(mapBC6, stat.type = c("seg.dist","prop","dxo","recomb"), layout = c(1,5), type = "l")

## ----ex18,echo=FALSE,fig.width=17,fig.height=12,warning=FALSE------------
profileMark(mapBC6, stat.type = c("seg.dist","prop","dxo","recomb"), layout = c(1,5), type = "l")

## ----ex19,eval=TRUE,echo=TRUE,prompt=TRUE--------------------------------
mapBC6 <- pushCross(mapBC6, type = "missing", pars = list(miss.thresh = 0.22, max.rf = 0.3))

## ----ex20,eval=FALSE,echo=TRUE,prompt=TRUE-------------------------------
#  heatMap(mapBC6, chr = c("L.3","L.5","L.8","L.9"), lmax = 70)

## ----ex21,echo=FALSE,fig.width=14,fig.height=8,warning=FALSE,cache=TRUE----
heatMap(mapBC6, chr = c("L.3","L.5","L.8","L.9"), lmax = 70)

## ----ex22,eval=TRUE,echo=TRUE,prompt=TRUE,cache=TRUE---------------------
mapBC6 <- mergeCross(mapBC6, merge = list("L.3" = c("L.3","L.5"), "L.8" = c("L.8","L.9")))
names(mapBC6$geno) <- paste("L.", 1:7, sep = "")
mapBC7 <- mstmap(mapBC6, bychr = TRUE, trace = TRUE, dist.fun = "kosambi", p.value = 2)
chrlen(mapBC7)

## ----ex23,eval=FALSE,echo=TRUE,prompt=TRUE-------------------------------
#  pg1 <- profileGen(mapBC7, bychr = FALSE, stat.type = c("xo","dxo","miss"), id = "Genotype", xo.lambda = 14, layout = c(1,3), lty = 2, cex = 0.7)

## ----ex24,echo=FALSE,fig.width=17,fig.height=10,warning=FALSE------------
pg1 <- profileGen(mapBC7, bychr = FALSE, stat.type = c("xo","dxo","miss"), id = "Genotype", xo.lambda = 14, layout = c(1,3), lty = 2, cex = 0.7)

## ----ex25,eval=TRUE,echo=TRUE,cache = TRUE,prompt=TRUE-------------------
mapBC8 <- subsetCross(mapBC7, ind = !pg1$xo.lambda)
mapBC9 <- mstmap(mapBC8, bychr = TRUE, dist.fun = "kosambi", trace = TRUE, p.value = 2)
chrlen(mapBC9)

## ----ex26,eval=FALSE,echo=TRUE,prompt=TRUE-------------------------------
#  profileMark(mapBC9, stat.type = c("seg.dist","prop","dxo","recomb"), layout = c(1,5), type = "l")

## ----ex27,echo=FALSE,fig.width=17,fig.height=10,warning=FALSE------------
profileMark(mapBC9, stat.type = c("seg.dist","prop"), layout = c(1,3), type = "l")

## ----ex28,eval=TRUE,echo=TRUE,prompt=TRUE--------------------------------
dm <- markernames(mapBC9, "L.2")[statMark(mapBC9, chr = "L.2", stat.type = "marker")$marker$neglog10P > 6]
mapBC10 <- drop.markers(mapBC9, dm)
mapBC11 <- pushCross(mapBC10, type = "seg.distortion", pars = list(seg.ratio = "70:30"))
mapBC12 <- mstmap(mapBC11, bychr = TRUE, trace = TRUE, dist.fun = "kosambi", p.value = 2)
round(chrlen(mapBC12) - chrlen(mapBC9), 5)
nmar(mapBC12) - nmar(mapBC10)

## ----ex29,eval=TRUE,echo=TRUE,prompt=TRUE--------------------------------
mapBC <- pushCross(mapBC12, type = "co.located")
names(mapBC)

## ----add1,eval=TRUE,echo=TRUE,prompt=TRUE--------------------------------
set.seed(123)
add1 <- drop.markers(mapBC, markernames(mapBC)[sample(1:3019, 2700, replace = FALSE)])
mapBCs <- drop.markers(mapBC, markernames(add1))
add3 <- add2 <- add1
add2 <- subset(add2, chr = "L.1")
add3$geno[[1]]$data <- pull.geno(add1)
add3$geno[[1]]$map <- 1:ncol(add3$geno[[1]]$data)
names(add3$geno[[1]]$map) <- markernames(add1)
names(add3$geno)[1] <- "ALL"
add3 <- subset(add3, chr = "ALL")

## ----add2,eval=TRUE,echo=TRUE,prompt=TRUE,cache=TRUE---------------------
add1 <- subset(add1, ind = 2:300)
full1 <- combineMap(mapBCs, add1, keep.all = TRUE)
full1 <- mstmap(full1, bychr = TRUE, trace = TRUE, anchor = TRUE, p.value = 2)

## ----add3,eval=TRUE,echo=TRUE,prompt=TRUE,cache=TRUE---------------------
add2 <- subset(add2, ind = 2:300)
full2 <- combineMap(mapBCs, add2, keep.all = TRUE)
full2 <- mstmap(full2, chr = "L.1", bychr = TRUE, trace = TRUE, anchor = TRUE, p.value = 2)

## ----add4,eval=TRUE,echo=TRUE,prompt=TRUE,cache=TRUE---------------------
add3 <- subset(add3, ind = 2:300)
full3 <- combineMap(mapBCs, add3, keep.all = TRUE)
full3 <- pushCross(full3, type = "unlinked", unlinked.chr = "ALL")
full3 <- mstmap(full3, bychr = TRUE, trace = TRUE, anchor = TRUE, p.value = 2)

## ----dist1,eval=FALSE,echo=TRUE,prompt=TRUE------------------------------
#  plot.missing(mapBC4)

## ----dist3,eval=TRUE,echo=TRUE,prompt=TRUE,cache=TRUE--------------------
mapBC4i <- mstmap(mapBC3, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-12, return.imputed = TRUE)
mapBC4i$geno[[1]]$map[1:14]
mapBC4i$imputed.geno[[1]]$map[1:5]

## ----dist4,echo = TRUE,prompt=TRUE---------------------------------------
len <- apply(mapBC4$geno[[1]]$data[,c(1,5)], 1, function(el)
             length(el[!is.na(el)]))
length(len[len > 1])
bca <- apply(mapBC4i$geno[[1]]$data[,c(1,5)], 1, function(el){
    el <- el[!is.na(el)]
    sum(abs(diff(el)))})
bca[bca > 0]

## ----dist5,echo=TRUE,prompt=TRUE-----------------------------------------
mapBC4i$imputed.geno[[1]]$data[pg$xo.lambda,1:5]

## ----dist6,echo=TRUE,prompt=TRUE-----------------------------------------
mapBC4e <- quickEst(mapBC4)
chrlen(mapBC4)
chrlen(mapBC4e)

## ----mvest1,eval=TRUE,echo=TRUE,prompt=TRUE,cache=TRUE-------------------
mapBC4a <- mstmap(mapBC3, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-12, mvest.bc = TRUE)
nmar(mapBC4)

## ----mvest2,eval=TRUE,echo=TRUE,prompt=TRUE------------------------------
sapply(mapBC4a$geno, function(el) length(unique(round(el$map, 4)))) - sapply(mapBC4$geno, function(el) length(unique(round(el$map, 4))))
chrlen(mapBC4a)

## ----mvest3,eval=TRUE,echo=TRUE,prompt=TRUE------------------------------
mapBC4b <- quickEst(mapBC4a)
chrlen(mapBC4b)

## ----dbd1,eval=TRUE,echo=TRUE,prompt=TRUE--------------------------------
mapBCd <- mapBC
mapBCd$geno <- lapply(mapBCd$geno, function(el){
    ns <- sample(1:ncol(el$data), ncol(el$data)/2, replace = TRUE)
    ns <- cbind(sample(1:nrow(el$data), ncol(el$data)/2, replace = TRUE), ns)
    el$data[ns] <- abs(1 - el$data[ns])
    el$data[el$data == 0] <- 2
    el})
mapBCd <- quickEst(mapBCd)
chrlen(mapBCd)

## ----dbd2,eval=TRUE,echo=TRUE,prompt=TRUE,cache=TRUE---------------------
mapBCda <- mstmap(mapBCd, bychr = TRUE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-12, detectBadData = TRUE)
chrlen(mapBCda)


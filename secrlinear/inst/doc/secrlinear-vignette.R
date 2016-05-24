## ----, eval = FALSE------------------------------------------------------
#  library(secrlinear)   ## also loads secr
#  library(maptools)     ## to read shapefiles
#  library(igraph)       ## to investigate edge lengths
#  options(digits = 4)   ## for more readable output

## ----, eval = TRUE, echo = FALSE, include=FALSE--------------------------
library(secrlinear)
library(maptools)     ## to read shapefiles, later
library(igraph)       ## to investigate edge lengths, later 
options(digits=4)     ## for more readable output

## ----, eval = TRUE-------------------------------------------------------
setwd(system.file("extdata", package = "secrlinear"))
arvicola <- read.capthist("Jun84capt.txt", "glymetrap.txt", covname = "sex")

## ----, eval = TRUE-------------------------------------------------------
setwd(system.file("extdata", package = "secrlinear"))
glymemask <- read.linearmask(file = "glymemap.txt", spacing = 4)

## ----,eval = FALSE, fig.width = 7, fig.height=3--------------------------
#  par(mar = c(1,1,1,1))
#  plot(glymemask)
#  plot(arvicola, add = TRUE, tracks = TRUE)
#  plot(traps(arvicola), add = TRUE)

## ----,eval = TRUE, warning=FALSE-----------------------------------------
## 2-D habitat, Euclidean distance
fit2DEuc <- secr.fit(arvicola, buffer = 200, trace = FALSE)   

## 1-D habitat, Euclidean distance
fit1DEuc <- secr.fit(arvicola, mask = glymemask, trace = FALSE)

## 1-D habitat, river distance
fit1DNet <- secr.fit(arvicola, mask = glymemask, trace = FALSE,
                details = list(userdist = networkdistance))

## ----, eval=TRUE---------------------------------------------------------
predict(fit2DEuc)
predict(fit1DEuc)
predict(fit1DNet)

## ----, eval = TRUE-------------------------------------------------------
setwd(system.file("extdata", package = "secrlinear"))
silverstreammask <- read.linearmask(file = "silverstream.shp", spacing = 50)
par(mar = c(1,1,1,1))
plot(silverstreammask)

## ----, eval = FALSE------------------------------------------------------
#  networklength <- sum(SpatialLinesLengths(attr(silverstreammask, "SLDF"))) / 1000
#  discrepancy <- networklength - masklength(silverstreammask)  ## km

## ----, eval = FALSE------------------------------------------------------
#  setwd(system.file("extdata", package = "secrlinear"))
#  silverstreamSLDF <- readShapeSpatial(fn = "silverstream")
#  silverstreammask <- read.linearmask(data = silverstreamSLDF, spacing = 50)

## ----, eval=TRUE---------------------------------------------------------
x <- seq(0, 4*pi, length = 200)
xy <- data.frame(x = x*100, y = sin(x)*300)
linmask <- read.linearmask(data = xy, spacing = 20)

## ----, eval = FALSE------------------------------------------------------
#  plot(linmask)

## ----, eval = FALSE------------------------------------------------------
#  ## start interactive session and click on two points
#  showpath(silverstreammask, lwd = 3)

## ----, eval = TRUE-------------------------------------------------------
trps <- make.line(linmask, detector = "proximity", n = 40, startbuffer = 0, by = 300,
                  endbuffer = 80, cluster = c(0,40,80), type = 'randomstart')

## ----, eval = FALSE------------------------------------------------------
#  plot(linmask)
#  plot(trps, add = TRUE, detpar = list(pch = 16, cex = 1.5, col='red'))

## ----, eval = FALSE------------------------------------------------------
#  plot(silverstreammask)
#  loc <- locator(30)
#  xy <- snapPointsToLinearMask(data.frame(loc), silverstreammask)
#  tr <- read.traps(data = xy, detector = 'multi')
#  plot(tr, add = TRUE)

## ----, eval = FALSE------------------------------------------------------
#  transects <- read.traps('transectxy.txt', detector = 'transect')
#  capt <- read.table('capt.txt')
#  tempCH <- make.capthist(capt, transects, fmt = 'XY')
#  tempCH <- snip(tempCH, by = 100)   ## for 100-m segments
#  CH <- reduce(tempCH, outputdetector = "count")

## ----eval = TRUE, echo = FALSE-------------------------------------------
setwd(system.file("extdata", package = "secrlinear"))
tr <- read.traps("silverstreamtraps.txt", detector = "multi")

## ----, eval = TRUE-------------------------------------------------------
## simulate population of 2 animals / km
pop <- sim.linearpopn(mask = silverstreammask, D = 2)
## simulate detections using network distances
CH <- sim.capthist(traps = tr, popn = pop, noccasions = 4,
                   detectpar = list(g0 = 0.25, sigma = 500), 
                   userdist = networkdistance)
summary(CH)    ## detector spacing uses Euclidean distances

## ----, eval=FALSE--------------------------------------------------------
#  ## and plot the simulated detections...
#  par(mar = c(1,1,1,1))
#  plot(silverstreammask)
#  plot(CH, add = TRUE, tracks = TRUE, varycol = TRUE, rad = 100, cappar = list(cex = 2))
#  plot(tr, add = TRUE)

## ----, eval = FALSE------------------------------------------------------
#  userd <- networkdistance(tr, silverstreammask)
#  userd[!is.finite(userd)] <- 1e8  ## testing
#  sfit <- secr.fit(CH, mask = silverstreammask, details = list(userdist = userd))
#  predict(sfit)

## ----, eval = TRUE-------------------------------------------------------
region.N(fit2DEuc)
region.N(fit1DNet)

## ----, eval = FALSE------------------------------------------------------
#  par(mfrow = c(1,2), mar = c(1,1,1,1))
#  plot(fit2DEuc$mask)
#  plot(traps(arvicola), add = TRUE)
#  mtext(side = 3,line = -1.5, "fit2DEuc$mask", cex = 1)
#  plot(fit1DNet$mask)
#  plot(traps(arvicola), add = TRUE)
#  mtext(side = 3,line = -1.5,"fit1DNet$mask", cex = 1)

## ----,eval = TRUE--------------------------------------------------------
derived(fit2DEuc)
derived(fit1DNet)

## ----, eval = FALSE------------------------------------------------------
#  ## interactively obtain LineID for central 'spine' by clicking on
#  ## each component line in plot
#  tmp <- getLineID(silverstreammask)
#  ## extract coordinates of 'spine'
#  spine <- subset(silverstreammask, LineID = tmp$LineID)
#  ## obtain network distances to spine and save for later use
#  netd <- networkdistance(spine, silverstreammask)  ## matrix dim = c(nrow(spine), nrow(mask))
#  dfs <- apply(netd, 2, min) / 1000  ## km
#  covariates(silverstreammask)$dist.from.spine <- dfs

## ----, eval = FALSE------------------------------------------------------
#  par(mar=c(1,1,1,4))
#  plot(silverstreammask, covariate = 'dist.from.spine', col = topo.colors(13),
#       cex = 1.5, legend = FALSE)
#  strip.legend('right', legend = seq(0, 6.5, 0.5), col = topo.colors(13),
#               title = 'dist.from.spine km', height = 0.35)
#  plot(spine, add = TRUE, linecol = NA, cex = 0.3)

## ----, eval = TRUE, strip.white = TRUE-----------------------------------
## initially OK (no movement > 1000 m)--
checkmoves(arvicola, mask = glymemask, accept = c(0,1000))
## deliberately break graph of linear mask
attr(glymemask, 'graph')[200:203,201:204] <- NULL
## no longer OK --
out <- checkmoves(arvicola, mask = glymemask, accept = c(0,1000))
## display captures of animals 32 and 35 whose records span break
out$df

## ----, eval = FALSE------------------------------------------------------
#  ## problem shows up where voles recaptured either side of break:
#  showedges(glymemask, col = 'red', lwd = 6)
#  plot(out$CH, add = TRUE, tracks = TRUE, rad=8,cappar=list(cex=1.5))
#  pos <- traps(arvicola)['560.B',]
#  text(pos$x+5, pos$y+80, 'break', srt=90, cex=1.1)

## ----, eval = FALSE------------------------------------------------------
#  plot(glymemask)
#  replot(glymemask)  ## click on corners to zoom in
#  showedges(glymemask, col = 'red', lwd = 2, add=T)
#  glymemask <- addedges(glymemask)

## ----, eval = FALSE------------------------------------------------------
#  par(mfrow = c(1,1), mar = c(1,1,1,5))
#  plot(silverstreammask)
#  centres <- data.frame(locator(4))
#  OK <- networkdistance(centres, silverstreammask) < 1000
#  for (i in 1:nrow(OK)) {
#      m1 <- subset(silverstreammask, OK[i,])
#      plot(m1, add = TRUE, col = 'red', cex = 1.7)
#      ml <- masklength(m1)
#      points(centres, pch = 16, col = 'yellow', cex = 1.4)
#      text (1406000, mean(m1$y), paste(ml, 'km'), cex = 1.2)
#  }
#  

## ----,eval=TRUE----------------------------------------------------------
library(secrdesign)

## create a habitat geometry
x <- seq(0, 4*pi, length = 200)
xy <- data.frame(x = x*100, y = sin(x)*300)
linmask <- read.linearmask(data = xy, spacing = 5)

## define two possible detector layouts
trp1 <- make.line(linmask, detector = "proximity", n = 80, start = 200, by = 30)
trp2 <- make.line(linmask, detector = "proximity", n = 40, start = 200, by = 60)
trplist <- list(spacing30 = trp1, spacing60 = trp2)

## create a scenarios dataframe
scen1 <- make.scenarios(D = c(50,200), trapsindex = 1:2, sigma = 25, g0 = 0.2)

## we specify the mask, rather than construct it 'on the fly',
## we will use a non-Euclidean distance function for both
## simulating detections and fitting the model...
det.arg <- list(userdist = networkdistance)
fit.arg <- list(details = list(userdist = networkdistance))

## run the scenarios and summarise results
sims1 <- run.scenarios(nrepl = 50, trapset = trplist, maskset = linmask, 
    det.args = list(det.arg), fit.args = list(fit.arg),
    scenarios = scen1, seed = 345, fit = FALSE)
summary(sims1)

## ----, eval=FALSE--------------------------------------------------------
#  sims2 <- run.scenarios(nrepl = 5, trapset = trplist, maskset = linmask,
#       det.args = list(det.arg), scenarios = scen1, seed = 345, fit = TRUE)
#  summary(sims2)

## ----, eval = FALSE------------------------------------------------------
#  ## It is efficient to pre-compute a matrix of distances between traps (rows)
#  ## and mask points (columns)
#  distmat <- networkdistance (traps(arvicola), glymemask, glymemask)
#  
#  ## Morning and evening trap checks as a time covariate
#  tcov <- data.frame(ampm = rep(c("am","pm"),3))
#  
#  glymefit1 <- secr.fit(arvicola, mask = glymemask, trace = FALSE,
#                        details = list(userdist = distmat),
#                        model = g0~1, hcov = "sex")
#  glymefit2 <- secr.fit(arvicola, mask = glymemask, trace = FALSE,
#                        details = list(userdist = distmat),
#                        model = g0~ampm, timecov = tcov, hcov = "sex")
#  glymefit3 <- secr.fit(arvicola, mask = glymemask, trace = FALSE,
#                        details = list(userdist = distmat),
#                        model = g0~ampm + h2, timecov = tcov, hcov = "sex")
#  glymefit4 <- secr.fit(arvicola, mask = glymemask, trace = FALSE,
#                        details = list(userdist = distmat),
#                        model = list(sigma~h2, g0~ampm + h2),
#                        timecov = tcov, hcov = "sex")
#  
#  fitlist <- secrlist(glymefit1, glymefit2, glymefit3, glymefit4)
#  ## dropping the detectfn (halfnormal) column to save space...
#  AIC(fitlist)[,-2]
#                                    model npar    logLik     AIC    AICc  dAICc AICcwt
#  secr4 D~1 g0~ampm + h2 sigma~h2 pmix~h2    7 -322.5414 659.083 665.305  0.000      1
#  secr3  D~1 g0~ampm + h2 sigma~1 pmix~h2    6 -347.3402 706.680 711.101 45.796      0
#  secr2       D~1 g0~ampm sigma~1 pmix~h2    5 -353.4807 716.961 719.961 54.656      0
#  secr1          D~1 g0~1 sigma~1 pmix~h2    4 -356.8020 721.604 723.509 58.204      0
#  
#  ## summaries of estimated density and sex ratio under different models
#  options(digits=3)
#  
#  ## model does not affect density estimate
#  collate(fitlist, perm = c(2,3,1,4))[,,1,"D"]
#        estimate SE.estimate  lcl  ucl
#  secr1     26.5        5.27 18.0 39.0
#  secr2     26.4        5.26 18.0 38.9
#  secr3     26.3        5.25 17.9 38.8
#  secr4     27.2        5.45 18.5 40.2
#  
#  ## model does affect the estimate of sex ratio (here proportion female)
#  collate(fitlist, perm=c(2,3,1,4))[,,1,"pmix"]
#        estimate SE.estimate   lcl   ucl
#  secr1    0.615      0.0954 0.421 0.779
#  secr2    0.615      0.0954 0.421 0.779
#  secr3    0.634      0.0938 0.439 0.793
#  secr4    0.669      0.0897 0.477 0.817
#  
#  ## predictions from best model
#  newdata <- expand.grid(ampm = c("am", "pm"), h2 = c("F", "M"))
#  predict(glymefit4, newdata = newdata)
#  
#  $`ampm = am, h2 = F`
#         link estimate SE.estimate    lcl    ucl
#  D       log   27.239      5.4478 18.477 40.158
#  g0    logit    0.218      0.0463  0.141  0.322
#  sigma   log   13.624      1.8764 10.414 17.823
#  pmix  logit    0.669      0.0897  0.477  0.817
#  
#  $`ampm = pm, h2 = F`
#         link estimate SE.estimate     lcl    ucl
#  D       log   27.239      5.4478 18.4768 40.158
#  g0    logit    0.116      0.0293  0.0694  0.186
#  sigma   log   13.624      1.8764 10.4136 17.823
#  pmix  logit    0.669      0.0897  0.4774  0.817
#  
#  $`ampm = am, h2 = M`
#         link estimate SE.estimate     lcl    ucl
#  D       log   27.239      5.4478 18.4768 40.158
#  g0    logit    0.153      0.0392  0.0908  0.246
#  sigma   log   70.958     10.0551 53.8247 93.545
#  pmix  logit    0.331      0.0897  0.1829  0.523
#  
#  $`ampm = pm, h2 = M`
#         link estimate SE.estimate     lcl    ucl
#  D       log  27.2394      5.4478 18.4768 40.158
#  g0    logit   0.0782      0.0201  0.0468  0.128
#  sigma   log  70.9581     10.0551 53.8247 93.545
#  pmix  logit   0.3311      0.0897  0.1829  0.523

## ----, eval = FALSE------------------------------------------------------
#  derived(glymefit4, distribution = 'binomial')
#      estimate SE.estimate   lcl   ucl    CVn     CVa    CVD
#  esa   0.9546          NA    NA    NA     NA      NA     NA
#  D    27.2362       2.867 22.17 33.46 0.1038 0.01748 0.1053


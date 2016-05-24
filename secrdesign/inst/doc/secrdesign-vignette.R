## ---- eval = TRUE, echo=FALSE, message=FALSE---------------------------------------
library(secrdesign)
## setwd('d:/density secr 2.10/secrdesign/vignettes')
load('runsims1.RData')
load('runsims2.RData')
options(width = 85)

## ---- eval = FALSE-----------------------------------------------------------------
#  library(secrdesign)
#  scen1 <- make.scenarios(D = c(5,10), sigma = 25, g0 = 0.2)
#  traps1 <- make.grid()
#  sims1 <- run.scenarios(nrepl = 50, trapset = traps1, scenarios =
#       scen1, seed = 345, fit = TRUE)

## ---- eval = TRUE------------------------------------------------------------------
summary(sims1)$OUTPUT

## ---- eval = TRUE------------------------------------------------------------------
library(secrdesign)
mydetectors <- list(grid6x6 = make.grid(6,6),
                    grid8x9 = make.grid(8,9),
                    grid12x12 = make.grid(12,12))

## ----eval=FALSE--------------------------------------------------------------------
#  make.scenarios (trapsindex = 1, noccasions = 3, nrepeats = 1, D, g0,
#      sigma, lambda0, detectfn = 0, recapfactor = 1, popindex = 1,
#      detindex = 1, fitindex = 1, groups, crosstraps = TRUE)

## ---- eval = TRUE------------------------------------------------------------------
make.scenarios (trapsindex = 1:3, noccasions = 4, D = 5, g0 = 0.2, sigma = c(20,30))

## ---- eval = TRUE------------------------------------------------------------------
make.scenarios (trapsindex = 1:3, noccasions = c(8,4,2), D = 5, g0 = 0.2,
                sigma = c(20,30), crosstraps = FALSE)

## ---- eval = FALSE-----------------------------------------------------------------
#  run.scenarios (nrepl, scenarios, trapset, maskset, xsigma = 4,
#      nx = 32, pop.args, det.args, fit = FALSE, fit.args, extractfn =
#      NULL, multisession = FALSE, ncores = 1, seed = 123)

## ----------------------------------------------------------------------------------
extractfn <- function(x) {
    if (inherits(x, 'capthist')) {
        ## summarised raw data
        counts <- function(CH) {
            ## for single-session CH
            if (nrow(CH)==0) 
                c(n=0, ndet=0, nmov=0, dpa = NA)
            else {
                nmoves <- sum(unlist(sapply(moves(CH), function(y) y>0)))
                ## detectors per animal
                dpa <- if (length(dim(CH)) == 2)
                    mean(apply(abs(CH), 1, function(y) length(unique(y[y>0]))))
                else
                    mean(apply(apply(abs(CH), c(1,3), sum)>0, 1, sum))
                c(n=nrow(CH), ndet=sum(abs(CH)>0), nmov=nmoves, dpa = dpa)
            }
        }
        if (ms(x)) 
            unlist(lapply(x, counts))
        else {
            gp <- covariates(x)$group
            if (is.null(gp)) 
                counts(x)
            else 
                unlist(lapply(split(x,gp,dropnullocc=TRUE), counts))
        }            
    }
    else if (inherits(x,'secr') & (!is.null(x$fit)))
        ## fitted model:
        ## default predictions of 'real' parameters
        predict(x)
    else
        ## null output: dataframe of 0 rows and 0 columns
        data.frame()   
}

## ---- eval = FALSE-----------------------------------------------------------------
#  closedNsim <- run.scenarios (nrepl = 10, scenarios = scen1, trapset = traps1,
#       extractfn = closedN, estimator = c("null", "chao", "chaomod"))

## ---- eval = FALSE-----------------------------------------------------------------
#  sum1 <- function(out) {
#      require(abind)
#      ## collapse replicates to an array, omitting non-numeric column
#      out <- do.call(abind, c(out, along = 3))[,-1,,drop = FALSE]
#      ## convert array from character to numeric
#      mode(out) <- "numeric"
#      ## take the average over replicates (meaningless for some fields)
#      apply(out, 1:2, mean, na.rm = TRUE)
#      }
#  lapply(closedNsim$output, sum1)

## ---- eval=FALSE-------------------------------------------------------------------
#  select.stats(object, parameter = "D", statistics)

## ---- eval=FALSE-------------------------------------------------------------------
#  find.param(object)

## ---- eval = TRUE------------------------------------------------------------------
stats1 <- select.stats(sims1, parameter = "D", statistics = c("estimate",
   "lcl", "ucl", "RB", "RSE", "COV"))
lapply(stats1$output, head, 4)

## ----eval=FALSE--------------------------------------------------------------------
#  x <- validate (x, test, validrange = c(0, Inf), targets = test)

## ----eval=FALSE--------------------------------------------------------------------
#  summary (object, dec = 5, fields = c("n", "mean", "se"), alpha = 0.05,
#      type = c("list", "dataframe", "array"), ...)

## ----echo=FALSE, comment=''--------------------------------------------------------
bullet <-  '*' ## rawToChar(as.raw(149))
tmp <- matrix(bullet, nr=8, ncol=12)
dimnames(tmp) <- list(Statistics = c("estimate","SE.estimate","lcl","ucl","RB","RSE","ERR","COV"),
Fields = c("n","mean","se","sd","min","max","lcl","ucl","rms","median","q025","q975"))
tmp[1:6,9] <- ""
tmp[3:5,7:8] <- ""
tmp[8,7:12] <- ""
print(tmp, qu=F)

## ---- eval = TRUE------------------------------------------------------------------
summary(stats1, c('n', 'mean', 'se', 'median'))

## ---- eval=FALSE-------------------------------------------------------------------
#  par(mfrow = c(2,2))
#  plot(stats1, type = "hist", statistic = "estimate")
#  plot(stats1, type = "CI")

## ---- echo=FALSE, eval=FALSE-------------------------------------------------------
#  png(file='d:/density secr 2.10/secrdesign/vignettes/secrdesign-fig3.png',
#      width = 850, height = 800)
#  par(mfrow = c(2,2), cex=1.2)
#  plot(stats1, type = "hist", statistic = "estimate")
#  plot(stats1, type = "CI")
#  dev.off()

## ---- eval = FALSE-----------------------------------------------------------------
#  sims2 <- run.scenarios(nrepl = 50, trapset = traps1, scenarios = scen1,
#      fit = TRUE, fit.args = list(method = "none"))

## ---- eval = TRUE------------------------------------------------------------------
summary(sims2)

## ---- eval = FALSE-----------------------------------------------------------------
#  scen3 <- make.scenarios(D = c(5,10), sigma = 25, g0 = 0.2)
#  traps3 <- make.grid()
#  raw3 <- run.scenarios(nrepl = 50, trapset = traps3, scenarios =
#      scen3, fit = FALSE, extractfn = identity)
#  summary(raw3)
#  ## fit and summarise models
#  sims3 <- fit.models(raw3, fit.args = list(list(model = g0~1),
#      list(model = g0~T)), fit = TRUE, ncores = 4)
#  summary(sims3)

## ---- eval = FALSE-----------------------------------------------------------------
#  traps4 <- list(grid6x6 = make.grid(6,6),
#                 grid8x9 = make.grid(8,9),
#                 grid12x12 = make.grid(12,12))
#  scen4 <- make.scenarios (trapsindex = 1:3, noccasions = c(8,4,2), D = 5,
#      g0 = 0.2, sigma = c(20,30), crosstraps = FALSE)
#  
#  sims4 <- run.scenarios(nrepl = 500, trapset = traps4, scenarios =
#       scen4, fit = FALSE, ncores = 3)

## ---- eval = TRUE------------------------------------------------------------------
class(sims4)        ## just peeking
find.stats(sims4)   ## just peeking
summary(sims4)

## ---- eval = FALSE-----------------------------------------------------------------
#  par(mfrow=c(4,3))
#  plot(sims4, statistic = "n", breaks = seq(0,80,5))      ## animals
#  plot(sims4, statistic = "nmov", breaks = seq(0,140,5))  ## movements

## ---- eval=FALSE, echo=FALSE-------------------------------------------------------
#  png(file='d:/density secr 2.10/secrdesign/vignettes/secrdesign-fig4.png',
#      width = 850, height = 800)
#  par(mfrow = c(4,3), cex=1)
#  plot(sims4, statistic = "n", breaks = seq(0,80,5))  ## number of animals
#  plot(sims4, statistic = "nmov", breaks = seq(0,140,5))
#  dev.off()

## ----eval = FALSE------------------------------------------------------------------
#  ## set up and run simulations
#  traps5 <- list(grid6x6 = make.grid(6,6),
#                 grid10x10 = make.grid(10,10))
#  scen5 <- make.scenarios (trapsindex = 1:2, noccasions = 5, D = 5,
#      g0 = 0.2, sigma = 25, recapfactor = c(0.5, 1, 2), fitindex = 1:2)
#  sims5 <- run.scenarios(nrepl = 500, trapset = traps5, scenarios =
#      scen5, fit = TRUE, fit.args = list(list(model = g0 ~ 1),
#      list(model = g0 ~ b)), ncores = 6)

## ---- eval = TRUE------------------------------------------------------------------
## select statistics and throw out any replicates with SE > 100
## (there is one -- see reduced n in output for scenario 11)
stats5 <- select.stats(sims5)
stats5 <- validate(stats5, "SE.estimate", c(0,100), "all")
sum5 <- summary(stats5, fields = c("n","mean","se","lcl","ucl", "median"))

## ---- eval = FALSE-----------------------------------------------------------------
#  ## plot
#  plot(c(0.5,6.5), c(-0.2,0.4), type = "n", xlab = "Scenario", ylab = "RB(D-hat)")
#  for (i in 1:12) {
#      xv <- if (i<=6) i else (i-6)+0.05
#      segments (xv, sum5$OUTPUT[[i]]["RB","lcl"], xv, sum5$OUTPUT[[i]]["RB","ucl"])
#      ptcol <- if (i<=6) "white" else "black"
#      points(xv, sum5$OUTPUT[[i]]["RB","mean"], pch = 21, bg = ptcol)
#  }
#  abline(h = 0, col="red")
#  text(c(1.5,3.5,5.5), rep(0.38,3), paste("recapfactor", c(0.5,1,2), sep = " = "))

## ---- eval = FALSE, echo=FALSE-----------------------------------------------------
#  ## plot
#  png(file='d:/density secr 2.10/secrdesign/vignettes/secrdesign-fig5.png',
#      width = 850, height = 500)
#  par(mar = c(4,4,1,1), cex = 1.4)
#  plot(c(0.5,6.5), c(-0.2,0.4), type = "n", xlab = "Scenario", ylab = "RB(D-hat)")
#  for (i in 1:12) {
#      xv <- if (i<=6) i else (i-6)+0.05
#      segments (xv, sum5$OUTPUT[[i]]["RB","lcl"], xv, sum5$OUTPUT[[i]]["RB","ucl"])
#      ptcol <- if (i<=6) "white" else "black"
#      points(xv, sum5$OUTPUT[[i]]["RB","mean"], pch = 21, bg = ptcol)
#  }
#  abline(h = 0, col="red")
#  text(c(1.5,3.5,5.5), rep(0.38,3), paste("recapfactor", c(0.5,1,2), sep = " = "))
#  dev.off()

## ---- eval = TRUE------------------------------------------------------------------
## look at extended output
sum5

## ---- eval = FALSE-----------------------------------------------------------------
#  
#  ## add covariates to builtin secr object possummask
#  ## D1 is homogeneous density
#  ## D2 is artificial SW - NE gradient in density
#  
#  xy <- apply(possummask,1,sum) / 500
#  covariates(possummask)[, "D1"] <- 2
#  covariates(possummask)[, "D2"] <- xy - mean(xy) + 2.5
#  
#  ## Note that this object already had a covariates dataframe
#  ## -- if it didn't we would use
#  ## covariates(possummask) <- data.frame ( D1 = ..., D2 = ...)
#  
#  ## specify scenarios
#  ## anticipate two different sets of arguments for sim.popn
#  ## with popindex = 1:2
#  
#  scen6 <- make.scenarios (g0 = 0.2, sigma = 45, noccasions = 5,
#      popindex = 1:2)
#  
#  ## specify alternate models for distribution of animals
#  
#  poplist <- list(list(model2D = "IHP", D = "D1"),
#                  list(model2D = "IHP", D = "D2"))
#  
#  ## run scenarios and summarise
#  ## we use the trap layout from the builtin secr object possumCH
#  
#  sims6 <- run.scenarios (500, scen6, traps(possumCH), possummask,
#      pop.args = poplist)

## ---- eval = TRUE------------------------------------------------------------------
summary(sims6)

## ----eval = FALSE------------------------------------------------------------------
#  sims6a <- run.scenarios (1, scen6, traps(possumCH), possummask,
#      pop.args = poplist, det.args = list(savepopn = TRUE),
#      extractfn = identity)

## ---- eval = FALSE-----------------------------------------------------------------
#  ## sims6a$output is now a list (one component per scenario) of lists
#  ## (one component per replicate) of simulated capthist objects, each
#  ## with its 'popn' object embedded as an attribute
#  
#  pop1 <- attr(sims6a$output[[1]][[1]], "popn")
#  pop2 <- attr(sims6a$output[[2]][[1]], "popn")
#  par(mfrow = c(1,2), mar=c(1,1,1,6))
#  plot(possummask, covariate = "D1", dots = FALSE, breaks = 0:6)
#  plot(traps(possumCH), detpar = list(col = 'green', pch = 15), add = TRUE)
#  plot(pop1, frame = FALSE, add = TRUE, col = "blue", pch = 16, cex = 0.6)
#  plot(possummask, covariate = 'D2', dots = FALSE, breaks = 0:6)
#  plot(traps(possumCH), detpar = list(col = 'green', pch = 15), add = TRUE)
#  plot(pop2, frame = FALSE, add = TRUE, col = "blue", pch = 16, cex = 0.6)

## ---- eval = FALSE-----------------------------------------------------------------
#  ## click on map to display height; Esc to exit
#  spotHeight(possummask, prefix = "D2")

## ---- eval = FALSE-----------------------------------------------------------------
#  pop1 <- attr(sims6a$output[[1]][[1]], "popn")
#  pop2 <- attr(sims6a$output[[2]][[1]], "popn")
#  png(file='d:/density secr 2.10/secrdesign/vignettes/secrdesign-fig6.png',
#      width=850, height=400)
#  par(mfrow = c(1,2), mar=c(1,1,1,6), cex=1.25)
#  plot(possummask, covariate = "D1", dots = FALSE, breaks = 0:6)
#  plot(traps(possumCH), detpar = list(col = 'green', pch = 15), add = TRUE)
#  plot(pop1, frame = FALSE, add = TRUE, col = "blue", pch = 16, cex = 0.6)
#  plot(possummask, covariate = 'D2', dots = FALSE, breaks = 0:6)
#  plot(traps(possumCH), detpar = list(col = 'green', pch = 15), add = TRUE)
#  plot(pop2, frame = FALSE, add = TRUE, col = "blue", pch = 16, cex = 0.6)
#  dev.off()

## ---- eval = FALSE-----------------------------------------------------------------
#  
#  library(secrlinear)
#  library(secrdesign)
#  
#  ## create a habitat geometry
#  x <- seq(0, 4*pi, length = 200)
#  xy <- data.frame(x = x*100, y = sin(x)*300)
#  linmask <- read.linearmask(data = xy, spacing = 5)
#  
#  ## define two possible detector layouts
#  trp1 <- make.line(linmask, detector = 'proximity', n = 80,
#                    startbuffer = 200, endbuffer = 200, by = 30)
#  trp2 <- make.line(linmask, detector = 'proximity', n = 40,
#                    startbuffer = 200, endbuffer = 200, by = 60)
#  trplist <- list(spacing30 = trp1, spacing60 = trp2)
#  
#  ## create a scenarios dataframe
#  scen7 <- make.scenarios(D = c(50,200), trapsindex = 1:2,
#                          sigma = 25, g0 = 0.2)
#  
#  ## we specify a mask, rather than construct it 'on the fly',
#  ## and must manually add column 'maskindex' to the scenarios
#  scen7$maskindex <- c(1,1)
#  
#  ## we will use a non-Euclidean distance function...
#  det.arg <- list(userdist = networkdistance)
#  
#  ## run the scenarios and summarise results
#  sims7 <- run.scenarios(nrepl = 500, trapset = trplist,
#      maskset = linmask, det.args = list(det.arg),
#      scenarios = scen7, seed = 345, fit = FALSE)

## ---- eval = TRUE------------------------------------------------------------------
summary(sims7)

## ---- eval = TRUE------------------------------------------------------------------
scen8 <- make.scenarios (D = 8, g0 = 0.3, sigma = 30, noccasions = c(4,8), groups = c('F','M'))
male <- scen8$group == 'M'
scen8$D[male] <- 4
scen8$g0[male] <- 0.2
scen8$sigma[male] <- 40
scen8[,1:8]

## ---- eval = TRUE------------------------------------------------------------------
grid <- make.grid(8, 8, spacing = 30)
mask <- make.mask(grid, buffer = 160, type = 'trapbuffer')
## extracts total density and proportion from output for the first group (F)
exfn <- function(x) {    
    if (inherits(x, 'secr') & !is.null(x$fit)) {
        pred <- predict(x)
        pred[[1]][c('D','pmix'),]
    }
    else data.frame()
}

## ---- eval = 2---------------------------------------------------------------------
raw8 <- run.scenarios(20, scen8, trapset = list(grid), fit = FALSE, maskset = list(mask))
summary(raw8)

## ---- eval = FALSE-----------------------------------------------------------------
#  sims8 <- run.scenarios(20, scen8, trapset = list(grid), fit = TRUE, extractfn = exfn,
#                         fit.args = list(model = list(g0~h2, sigma~h2), hcov = 'group'),
#                         maskset = list(mask))

## ---- eval = TRUE------------------------------------------------------------------
summary(select.stats(sims8,'D'))$OUTPUT
summary(select.stats(sims8,'pmix'))$OUTPUT

## ---- eval = FALSE, echo = FALSE---------------------------------------------------
#  save(sims1, sims2,
#       file = 'd:/density secr 2.10/secrdesign/vignettes/runsims1.RData')

## ---- eval = FALSE, echo = FALSE---------------------------------------------------
#  save(raw3, sims3, sims4, sims5, sims6, sims6a, sims7, raw8, sims8,
#       file = 'd:/density secr 2.10/secrdesign/vignettes/runsims2.RData')


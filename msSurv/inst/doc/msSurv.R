### R code from vignette source 'msSurv.Rnw'

###################################################
### code chunk number 1: prelim
###################################################
options(prompt="R> ", continue = "+  ", useFancyQuotes = FALSE)


###################################################
### code chunk number 2: <msSurv
###################################################
library("msSurv")


###################################################
### code chunk number 3: RCdata
###################################################
data("RCdata")


###################################################
### code chunk number 4: ex1dat
###################################################
RCdata[70:76, ]


###################################################
### code chunk number 5: stateinfo
###################################################
Nodes <- c("1", "2", "3", "4", "5")
Edges <- list("1" = list(edges = c("2", "3")),
              "2" = list(edges = c("4", "5")),
              "3" = list(edges = NULL),
              "4" = list(edges = NULL),
              "5" = list(edges = NULL))


###################################################
### code chunk number 6: tree
###################################################
treeobj <- new("graphNEL", nodes = Nodes, edgeL = Edges,
           edgemode = "directed")


###################################################
### code chunk number 7: exrtcens
###################################################
ex1 <- msSurv(RCdata, treeobj, bs = TRUE)


###################################################
### code chunk number 8: printex1
###################################################
print(ex1)


###################################################
### code chunk number 9: msSurv.Rnw:666-667
###################################################
Pst(ex1, s = 1, t = 3.1)


###################################################
### code chunk number 10: addvar (eval = FALSE)
###################################################
## Pst(ex1, s = 1, t = 3.1, covar = TRUE)


###################################################
### code chunk number 11: msSurv.Rnw:685-686
###################################################
SOPt(ex1, t = 0.85)


###################################################
### code chunk number 12: msSurv.Rnw:696-697 (eval = FALSE)
###################################################
## SOPt(ex1, t = 0.85, covar = TRUE)


###################################################
### code chunk number 13: EntryExit
###################################################
EntryExit(ex1, t = 1)


###################################################
### code chunk number 14: EntryExitObjs
###################################################
norm.dist <- EntryExit(ex1, t = 1, var = TRUE)
sub.dist  <- EntryExit(ex1, t = 1, norm = FALSE, var = TRUE)


###################################################
### code chunk number 15: EntryExitOutput
###################################################
names(norm.dist)
names(sub.dist)


###################################################
### code chunk number 16: LTRCata
###################################################
data("LTRCdata")


###################################################
### code chunk number 17: head
###################################################
LTRCdata[489:494, ]


###################################################
### code chunk number 18: LTRCtree
###################################################
Nodes <- c("1", "2", "3")
Edges <- list("1" = list(edges = c("2", "3")),
              "2" = list(edges = c("3")),
              "3" = list(edges = NULL))
LTRCtree <- new("graphNEL", nodes = Nodes, edgeL = Edges,
                edgemode = "directed")


###################################################
### code chunk number 19: exrtcenslt
###################################################
ex2 <- msSurv(LTRCdata, LTRCtree, LT = TRUE)


###################################################
### code chunk number 20: msSurv.Rnw:857-858
###################################################
summary(ex2, digits = 2)


###################################################
### code chunk number 21: sumAll (eval = FALSE)
###################################################
## summary(ex2, all = TRUE)


###################################################
### code chunk number 22: sumSOP (eval = FALSE)
###################################################
## summary(ex2, trans.pr = FALSE, dist = FALSE)


###################################################
### code chunk number 23: configdat
###################################################
if (require("KMsurv")) {
    data("bmt")
    attach(bmt)

    ## NOTES
    ## ta = Time To Acute Graft-Versus-Host Disease
    ## tp = Time To Platelet Recovery (incorrect in help file for 'bmt')
    ## tc = Time To Chronic Graft-Versus-Host Disease
    ## da = Acute GVHD Indicator 1-Developed Acute GVHD 0-Never Developed Acute GVHD)
    ## dp = Platelet Recovery Indicator 1-Platelets Returned To Normal, 0-Platelets Never Returned to Normal
    ## dc = Chronic GVHD Indicator 1-Developed Chronic GVHD 0-Never Developed Chronic GVHD
    ## t2 = Disease Free Survival Time (Time To Relapse, Death Or End Of Study)
    ## t2 used for CENSORING time here

    ## Determining transitions from STATE 1 ###

    ## 1 -> 2 transition:
    ## aGVHD BEFORE platelets return to normal
    s1to2 <- which(ta < tp) #7
    data12 <- data.frame(id = s1to2, stop = ta[s1to2], start.stage = 1, end.stage = 2)

    ## 1 -> 3 transition:
    ## platelets return to normal BEFORE aGVHD (and before cGVHD ... )
    s1to3 <- which(tp < ta & tp < tc) #117
    data13 <- data.frame(id = s1to3, stop = tp[s1to3], start.stage = 1, end.stage = 3)
    ## censored individuals = death or relapse, use t2
    c1 <- which(t2 <= pmin(ta, tp)) ## 13
    data10 <- data.frame(id = c1, stop = t2[c1], start.stage = 1, end.stage = 0)


    ## Determining transitions from STATE 2 ###

    ## 2 -> 4 transition:
    ## Individuals who developed acute GVHD BEFORE platelets returned to normal AND
    ## platelets returned to normal BEFORE chronic GVHD AND not censored
    s2to4 <- which(ta < tp & tp < tc & tp < t2 & dp==1) ## 3 (20, 32, 107) OK
    data24 <- data.frame(id = s2to4, stop = tp[s2to4], start.stage = 2, end.stage = 4)

    ## 2 -> 5 transition:
    ## Individuals who developed acute GVHD before platelets return to normal AND
    ## develop chronic GVHD BEFORE platelets returned to normal AND not censored
    s2to5 <- which(ta < tp & tc < tp & tc < t2 & dc==1)
    data25 <- data.frame(id = s2to5, stop = tc[s2to5], start.stage = 2, end.stage = 5)

    ## 2 -> 0 transition (censored):
    c2 <- which(t2 > ta & t2 <= pmin(tp, tc)) ## 2,  censored in stage 2
    data20 <- data.frame(id = c2, stop = t2[c2], start.stage = 2, end.stage = 0)


    ## Determining transitions from STATE 3 ##

    ## 3 -> 4 transition:
    ## normal platelets BEFORE aGVHD, then aGVHD BEFORE end of study and cGVHD
    s3to4 <- which(tp < ta & ta < tc & ta < t2  & da==1 ) ## 19 OK
    ## 11 of these have dc=1
    data34 <- data.frame(id = s3to4, stop = ta[s3to4], start.stage = 3, end.stage = 4)

    ## 3 -> 5 transition:
    ## normal platelets before aGVHD, and cGVHD BEFORE end of study and aGVHD
    s3to5 <- which(tp < ta & tc < ta & tc < t2 & dc==1) ## 44
    data35 <- data.frame(id = s3to5, stop = tc[s3to5], start.stage = 3, end.stage = 5)

    ## 3 -> 0 transition (censored):
    c3 <- which(t2 > tp & t2 <= pmin(ta, tc))
    data30 <- data.frame(id = c3, stop = t2[c3], start.stage = 3, end.stage = 0)


    ## Determining transitions from STATE 4 ##

    ## 4 -> 5 transition:
    ## TWO ways to reach state 4, via state 2 or state 3
    ## via state 2 to 4 to 5
    s4to5a <- which((ta < tp & tp < tc & tp < t2 & dp==1) & tc < t2 & dc==1) ## 1
    ## via state 3 to 4 to 5
    s4to5b <- which((tp < ta & ta < tc & ta < t2  & da==1) & tc < t2 & dc==1) ## 11
    data45a <- data.frame(id = s4to5a, stop = tc[s4to5a], start.stage = 4, end.stage = 5)
    data45b <- data.frame(id = s4to5b, stop = tc[s4to5b], start.stage = 4, end.stage = 5)

    ## 4 -> 0 transition
    ## via state 2 to 4 to 5
    c4a <- which((ta < tp & tp < tc & tp < t2 & dp==1) & t2 <= tc)
    ## via state 3 to 4 to 5
    c4b <- which((tp < ta & ta < tc & ta < t2  & da==1) & t2 <= tc)
    ## identical to below ...
    ## c4a <-  s2to4[which(!(s2to4%in%s4to5a))] ## 2 (20, 107)
    ## c4b <-  s3to4[which(!(s3to4%in%s4to5b))]
    data40a <- data.frame(id = c4a, stop = t2[c4a], start.stage = 4, end.stage = 0)
    data40b <- data.frame(id = c4b, stop = t2[c4b], start.stage = 4, end.stage = 0)


    ## Combining to make a DATASET ready for msSurv ##
    data.sdc <- rbind(data10, data12, data13, data20, data24, data25,
                      data30, data34, data35, data40a, data40b, data45a, data45b)
    data.sdc <- data.sdc[order(data.sdc$id),  ]

    ## Remove row with STOP time of 0
    ## ID = 124 - basically individual STARTS in stage 3 (platelets normal at t = 0)
    data.sdc <- data.sdc[-which(data.sdc$stop==0),  ]

    detach(bmt)
}



###################################################
### code chunk number 24: depdataframe
###################################################
if(exists("data.sdc")) head(data.sdc)


###################################################
### code chunk number 25: deptree
###################################################
Nodes <- c("1", "2", "3", "4", "5")
Edges <- list("1" = list(edges = c("2", "3")),
              "2" = list(edges = c("4", "5")),
              "3" = list(edges = c("4", "5")),
              "4" = list(edges = c("5")),
              "5" = list(edges = NULL))
deptree <- new("graphNEL", nodes = Nodes, edgeL = Edges,
               edgemode = "directed")


###################################################
### code chunk number 26: depcall
###################################################
if(exists("data.sdc"))
    DepEx <- msSurv(data.sdc, deptree, cens.type = "dep", bs = TRUE)


###################################################
### code chunk number 27: DepExSOP
###################################################
if(exists("DepEx")) plot(DepEx)


###################################################
### code chunk number 28: plotst2 (eval = FALSE)
###################################################
## if(exists("DepEx")) plot(DepEx, state = "2")


###################################################
### code chunk number 29: plotst23 (eval = FALSE)
###################################################
## if(exists("DepEx")) plot(DepEx, state = c("2","3"))


###################################################
### code chunk number 30: DepExTransPlot
###################################################
if(exists("DepEx"))
    plot(DepEx, plot.type = "transprob")


###################################################
### code chunk number 31: DepExTransPlot1 (eval = FALSE)
###################################################
## if(exists("DepEx"))
##     plot(DepEx, plot.type = "transprob", trans = c("1 2", "1 3"))


###################################################
### code chunk number 32: EntryExitDepExData (eval = FALSE)
###################################################
## ## NOTE - For state entry / exit distributions for state j we need to
## ##  define states so that only path between states before and after is
## ##  THROUGH state j, and no recursion to state j
## if(exists("data.sdc")) {
##     data.sdc2 <- data.sdc
##     data.sdc2$end.stage <- with(data.sdc2, ifelse(start.stage == 3 & end.stage == 4,
##                                                   "4-3", end.stage))
##     data.sdc2$end.stage <- with(data.sdc2, ifelse(start.stage == 3 & end.stage == 5,
##                                                   "5-3", end.stage))
##     data.sdc2$end.stage <- with(data.sdc2, ifelse(start.stage == 2 & end.stage == 4,
##                                                   "4-2", end.stage))
##     data.sdc2$end.stage <- with(data.sdc2, ifelse(start.stage == 2 & end.stage == 5,
##                                                   "5-2", end.stage))
##     ## table(data.sdc2$end.stage)
##     idx <- which(data.sdc2$end.stage == "5")
##     data.sdc2$start.stage[idx] <- data.sdc2$end.stage[(idx-1)]
##     data.sdc2$end.stage[idx] <- with(data.sdc2[idx, ], ifelse(start.stage == "4-3", "5-3", "5-2"))
## 
##     ## Drop individuals with start.stage == 4, since has no effect on
##     ##      entry / exit from states 2 and 3
##     idx <- which(data.sdc2$start.stage == "4")
##     data.sdc2 <- data.sdc2[-idx, ]
## 
##     Nodes2 <- c("1", "2", "4-2", "5-2", "3", "4-3", "5-3")
##     ## NOTE - edges give the INDEX values of states
##     Edges2 <- list("1" = list(edges = c(2, 5)),  ## states 2, 3
##                    "2" = list(edges = c(3, 4)),  ## states "4-2", "5-2"
##                    "4-2" = list(edges = c(4)),
##                    "5-2" = list(edges = NULL),
##                    "3" = list(edges = c(6, 7)),
##                    "4-3" = list(edges = c(7)),
##                    "5-3" = list(edges = NULL))
##     deptree2 <- new("graphNEL", nodes = Nodes2, edgeL = Edges2,
##                     edgemode = "directed")
## }


###################################################
### code chunk number 33: EntryExitDepExStateModel (eval = FALSE)
###################################################
## if(exists("deptree2") & require("Rgraphviz"))
##     plot(deptree2)


###################################################
### code chunk number 34: EntryExitDepExFit (eval = FALSE)
###################################################
## if(exists("data.sdc2"))
##     DepEx2 <- msSurv(data.sdc2, deptree2, cens.type = "dep", bs = TRUE)


###################################################
### code chunk number 35: NormEntryDepExPlot (eval = FALSE)
###################################################
## ## Now plot normalized entry distributions for states 2, 3
## if(exists("DepEx2"))
##     plot(DepEx2, plot.type = "entry.norm", states = c("2", "3"),
##          xlim = c(0, 100))


###################################################
### code chunk number 36: NormExitDepExPlot (eval = FALSE)
###################################################
## ## Now plot normalized exit distributions for states 2, 3
## if(exists("DepEx2"))
##     plot(DepEx2, plot.type = "exit.norm", states = c("2", "3"))
## ## Note that while normalized entry distributions will by
## ##  definition -> 1 as t -> infty, the normalized exit distributions
## ##  need not do this (e.g., if some individuals are never observed to
## ##  transition out of that state)


###################################################
### code chunk number 37: StateDependentCensoringIllnessDeath (eval = FALSE)
###################################################
## 
## library("MASS")  ## only needed for mvrnorm function
## 
## ####################################################################
## ## Simulation based on Datta and Satten 2002 Biometrics paper
## ####################################################################
## 
## ## 3-state illness-death model
## Nodes <- c("1", "2", "3")
## Edges <- list("1"=list(edges=c("2", "3")), "2"=list(edges=c("3")),
##            "3"=list(edges=NULL))
## SDCtree <- new("graphNEL", nodes=Nodes, edgeL=Edges, edgemode="directed")
## 
## 
## samp.size <- 1000
## mu <- c(1.81, 1.7, 1.7)
## cov <- matrix(0.91, 3, 3)
## diag(cov) <- rep(0.92, 3)
## cov <- 0.393*cov
## set.seed(101)
## data <- mvrnorm(samp.size, mu, cov)
## ## w1, w2, w3,
## ## if exp(w1) < exp(w2) -> death before illness at time exp(w1)
## ## else illness at time exp(w2) and death at time exp(w1) + exp(w2)
## 
## data <- as.data.frame(data)
## names(data) <- c("w1", "w2", "w3")
## idx <- which(data$w1 < data$w2)
## data13 <- data.frame(id=idx, stop=exp(data$w1[idx]), start.stage=1, end.stage=3)
## 
## ## now folks with illness ...
## idx <- which(data$w1 > data$w2)
## data12 <- data.frame(id = idx, stop = exp(data$w2[idx]), start.stage = 1, end.stage = 2)
## data23 <- data.frame(id = idx, stop = exp(data$w2[idx]) + exp(data$w3[idx]), start.stage = 2, end.stage = 3)
## 
## data.nocens <- rbind(data13, data12, data23)
## data.nocens <- data.nocens[order(data.nocens$id, data.nocens$stop), ]
## 
## ## no censoring estimates
## est.nocens <- msSurv(data.nocens, SDCtree)
## 
## 
## ####################################################################
## ## Now incorporate censoring ...
## ## Weibull distribution:
## ####################################################################
## 
## ####################################################################
## ## State 1: hazard = t^(1/4)/60
## ####################################################################
## ## In R notation:
## ## a = shape = 1+1/4
## ## b = scale = (a*den)^1/a, where den = 60
## a <- 1+1/4
## den <- 60
## 
## c1 <- rweibull(samp.size, shape = a, scale = (a*den)^(1/a))
## idx <- which(c1 < pmin(exp(data$w1), exp(data$w2)))
## length(idx)/samp.size
## ## ~13%
## 
## ####################################################################
## ## State 2: hazard = t/50
## ####################################################################
## 
## a <- 2    ## shape
## den <- 50 ## "adjustment" = denominator of hazard
## ## NOTE - Need to sample from TRUNCATED distribution
## ##        Use inverse-transform method
## ## Just generate uniform restricted to cdf at truncation point to 1,
## ##       then take inverse CDF at that point ...
## idx <- which(exp(data$w2) < pmin(c1, exp(data$w1)))
## lower <- pweibull(exp(data$w2)[idx], shape = a, scale = (a*den)^(1/a))
## u <- runif(length(idx), lower, 1)
## c2 <- qweibull(u, shape = a, scale = (a*den)^(1/a))
## 
## idx2 <- which(c2 < (exp(data$w2)[idx] + exp(data$w3)[idx]))
## length(idx2)/length(idx)
## ## ~58.5%
## 
## 
## ## Now incorporate censoring
## 
## ## 1->0
## idx <- which(c1 < pmin(exp(data$w1), exp(data$w2)))
## ## 130 = 13%
## data10 <- data.frame(id=idx, start = 0, stop=c1[idx], start.stage=1, end.stage=0)
## 
## ## 1->3
## idx <- which(exp(data$w1) < pmin(c1, exp(data$w2)))
## data13 <- data.frame(id=idx, start = 0, stop=exp(data$w1[idx]), start.stage=1, end.stage=3)
## 
## ## 1->2
## idx <- which(exp(data$w2) < pmin(c1, exp(data$w1)))
## data12 <- data.frame(id=idx, start = 0, stop=exp(data$w2[idx]), start.stage=1, end.stage=2)
## 
## ## 2->0 (censored in stage 2)
## idx <- which(exp(data$w2) < pmin(c1, exp(data$w1)))
## sum(c2 > exp(data$w2)[idx])  ## everyone
## idx2 <- which(c2 < (exp(data$w2)[idx] + exp(data$w3)[idx]))
## length(idx2)/length(idx)
## ## ~58.5%
## data20 <- data.frame(id=idx[idx2], start=exp(data$w2[idx][idx2]), stop = c2[idx2], start.stage=2, end.stage=0)
## 
## ## 2->3
## idx2 <- which(c2 > (exp(data$w2)[idx] + exp(data$w3)[idx]))
## data23 <- data.frame(id = idx[idx2], start=exp(data$w2[idx][idx2]),
##                      stop = exp(data$w2[idx][idx2]) + exp(data$w3[idx][idx2]),
##                      start.stage = 2, end.stage = 3)
## 
## 
## data.cens <- rbind(data10, data13, data12, data20, data23)
## data.cens <- data.cens[order(data.cens$id, data.cens$stop), ]
## 
## 
## ####################################################################
## ## Assuming independence
## ####################################################################
## est.ind <- msSurv(data.cens, SDCtree)
## 
## 
## ####################################################################
## ## Assuming dependence
## ####################################################################
## est.dep <- msSurv(data.cens, SDCtree, cens = "dep")
## 


###################################################
### code chunk number 38: StateDependentCensoringIllnessDeathPlot (eval = FALSE)
###################################################
## ####################################################################
## ## PLOT of STATE 2 OCCUPATION PROBABILITY
## ####################################################################
## 
## 
## ## All on same graph
## s2.ind <- ps(est.ind)[,2]
## s2.nocens <- ps(est.nocens)[,2]
## s2.dep <- ps(est.dep)[,2]
## 
## t.ind <- et(est.ind)
## t.nocens <- et(est.nocens)
## t.dep <- et(est.dep)
## 
## plot(t.nocens, s2.nocens, type = "s", xlim = c(0, 50),
##      ylab = "Probability", xlab = "Time",
##      main = "State 2 occupation probability")
## lines(t.ind, s2.ind, col = 2, type = "s")
## lines(t.dep, s2.dep, col = 3, type = "s")
## legend("topright", c("No censoring", "AJ - Independent", "DS - Dependent"),
##        lty = 1, col = 1:3, bty = "n")
## 



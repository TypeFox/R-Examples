## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, prompt=TRUE)

## ------------------------------------------------------------------------
2 + 2

## ------------------------------------------------------------------------
my.variable <- 2 + 2
my.variable * 3

## ----eval=FALSE----------------------------------------------------------
#  install.packages("optmatch")
#  install.packages("RItools")

## ----echo=FALSE,message=FALSE--------------------------------------------
library(optmatch)
library(RItools)

## ----eval=FALSE----------------------------------------------------------
#  library(optmatch)
#  library(RItools)

## ---- echo=FALSE---------------------------------------------------------
data(nuclearplants)

## ---- eval=FALSE---------------------------------------------------------
#  data(nuclearplants)

## ------------------------------------------------------------------------
head(nuclearplants)

## ---- eval=FALSE---------------------------------------------------------
#  help("nuclearplants")

## ----eval=FALSE----------------------------------------------------------
#  nuclearplants$pt
#  table(nuclearplants$pt)
#  with(nuclearplants, table(pt))

## ------------------------------------------------------------------------
nuke.nopt <- subset(nuclearplants, pt == 0)

## ----eval=FALSE----------------------------------------------------------
#  head(nuke.nopt)
#  tail(nuke.nopt)

## ------------------------------------------------------------------------
table(nuke.nopt$pr)

## ------------------------------------------------------------------------
pairmatch(pr ~ cap, data = nuke.nopt)

## ------------------------------------------------------------------------
print(pairmatch(pr ~ cap, data = nuke.nopt), grouped = TRUE)

## ---- results="asis", echo=FALSE, warning=FALSE--------------------------
library(pander)
a <- with(nuke.nopt, data.frame(
                         Plant=row.names(nuke.nopt),
                         Date=round(date-65, 1),
                         Capacity=round(x=(cap-400),digits=-1))[as.logical(pr),])

b <- with(nuke.nopt, data.frame(
                         Plant=row.names(nuke.nopt),
                         Date=round(date-65, 1),
                         Capacity=round(x=(cap-400),digits=-1))[!as.logical(pr),])

rownames(a) <- NULL
rownames(b) <- NULL

c <- cbind(data.frame(rbind(as.matrix(a), matrix(nrow=nrow(b)-nrow(a), ncol=3))), b)
pandoc.table(c, style="multiline", missing="",
             caption='New-site (left columns) versus existing-site (right columns) plants. "date" is `date-65`; "capacity" is `cap-400`.')

## ---- eval=FALSE---------------------------------------------------------
#  summary(pairmatch(pr ~ cap, data = nuke.nopt))

## ------------------------------------------------------------------------
pm <- pairmatch(pr ~ cap, data = nuke.nopt)

## ---- eval=FALSE---------------------------------------------------------
#  summary(lm(cost ~ pr + pm, data = nuke.nopt))

## ------------------------------------------------------------------------
tm <- pairmatch(pr ~ cap, controls = 2, data = nuke.nopt)

## ---- eval=FALSE---------------------------------------------------------
#  pairmatch(pr ~ cap, controls = 3, data=nuke.nopt)

## ---- error=TRUE---------------------------------------------------------
pairmatch(pr ~ cap + cost, caliper=.001, data = nuke.nopt)

## ------------------------------------------------------------------------
summary(pm)

## ---- eval=FALSE---------------------------------------------------------
#  cap.noadj <- lm(cap ~ pr, data = nuke.nopt)
#  summary(cap.noadj)

## ------------------------------------------------------------------------
summary(lm(cap ~ pr, data = nuke.nopt))$coeff["pr",]

## ------------------------------------------------------------------------
summary(lm(cap ~ pr + pm, data = nuke.nopt))$coeff["pr",]

## ------------------------------------------------------------------------

xBalance(pr ~ cap + t2, report="all", data=nuke.nopt)
xBalance(pr ~ cap + t2 + strata(pm),
         data=nuke.nopt,
         report=c("adj.mean.diffs", "std", "z"))

## ------------------------------------------------------------------------
psm <- glm(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n + pt,
           family = binomial, data = nuclearplants)

## ----fig.width=5, fig.height=5-------------------------------------------
boxplot(psm)

## ------------------------------------------------------------------------
ps.pm <- pairmatch(psm, data = nuclearplants)
summary(ps.pm)

## ------------------------------------------------------------------------
psm.dist <- match_on(psm, data=nuclearplants)

## ------------------------------------------------------------------------
caliper(psm.dist, 2)

## ------------------------------------------------------------------------
ps.pm2 <- pairmatch(psm.dist, data = nuclearplants)
ps.pm3 <- pairmatch(psm.dist + caliper(psm.dist, 2), data = nuclearplants)
all.equal(ps.pm, ps.pm2, check.attributes=FALSE)
all.equal(ps.pm, ps.pm3, check.attributes=FALSE)
summary(ps.pm3)

## ------------------------------------------------------------------------
mhd1 <- match_on(pr ~ date + cap + scores(psm), data=nuclearplants)
mhpc.pm <- pairmatch(mhd1, caliper=1, data=nuclearplants)
summary(mhpc.pm) # oops
mhpc.pm <- pairmatch(mhd1, caliper=2, data=nuclearplants)
summary(mhpc.pm) # better!

## ---- eval=FALSE---------------------------------------------------------
#  library(RItools)
#  xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n, data = nuclearplants)
#  xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n + pt +
#               strata(ps.pm2) -1, # the `-1` just focuses the output a little
#           data = nuclearplants)

## ---- fig.width=5, fig.height=5------------------------------------------
myb <- xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n + 
                strata(ps.pm2),
                data = nuclearplants,
                report = c("adj.means", "std.diffs",
                           "z.scores", "chisquare.test"))
plot(myb)
print(myb, digits=1)

## ------------------------------------------------------------------------
summary(ps.pm2, psm)

## ---- eval=FALSE---------------------------------------------------------
#  summary(fullmatch(pr ~ date + cap, data = nuke.nopt))
#  summary(fullmatch(pr ~ date + cap, data = nuke.nopt, min = 1))
#  summary(fullmatch(pr ~ date + cap, data = nuke.nopt, min = 2, max = 3))

## ------------------------------------------------------------------------
pairmatch(pr ~ date + cap + scores(psm), data=nuclearplants)
pairmatch(pr ~ date + cap + scores(psm) + strata(pt), data=nuclearplants)

## ------------------------------------------------------------------------
cap.dist <- match_on(pr ~ cap, data = nuke.nopt)
pm1 <- pairmatch(pr ~ cap, data=nuke.nopt)
pm2 <- pairmatch(cap.dist, data=nuke.nopt)
all.equal(pm1, pm2, check.attributes = FALSE)
summary(pm2)

## ------------------------------------------------------------------------
round(cap.dist[1:3, 1:3], 1)

## ------------------------------------------------------------------------
round(cap.dist + caliper(cap.dist, 2), 1)

## ---- eval=FALSE---------------------------------------------------------
#  pairmatch(cap.dist + caliper(cap.dist, 2), data = nuke.nopt)

## ---- eval=FALSE---------------------------------------------------------
#  library(foreign)
#  ?read.dta
#  ?write.dta

## ---- eval=FALSE---------------------------------------------------------
#  my.plants <- read.csv("nuclearplants.csv", header = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  plant.match <- pairmatch(pr ~ cap, data = my.plants)
#  my.plants.extended <- data.frame(my.plants, matches = plant.match, check.rows=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  write.csv(my.plants.extended, file = "nuclearplants-with-matches.csv")

## ------------------------------------------------------------------------
data(tli, package="xtable")
head(tli)

## ---- eval=FALSE---------------------------------------------------------
#  download.file("http://www-stat.wharton.upenn.edu/~rosenbap/DOSdata.RData",
#                destfile="./DOSdata.RData")
#  load("./DOSdata.RData")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("arm", dep=T) # if not already installed
#  data(lalonde, package="arm")
#  help("lalonde", package="arm")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("Hmisc", dep=T) # if not already installed
#  Hmisc:::getHdata(rhc, what = "all")


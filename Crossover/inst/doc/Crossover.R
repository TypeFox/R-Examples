## ----OptionsAndLibraries, include=FALSE, message=FALSE------------------------
if (exists("opts_chunk")) {
  opts_chunk$set(concordance=TRUE)
  opts_chunk$set(tidy.opts=list(keep.blank.line=FALSE, width.cutoff=80))
  #opts_chunk$set(size="footnotesize")
  #opts_chunk$set(size="tiny")
  opts_chunk$set(size="scriptsize") 
  opts_chunk$set(cache=TRUE)
  opts_chunk$set(autodep=TRUE)
}

# See http://yihui.name/knitr/hooks for the following code.
CrossoverNamespace <- function(before, options, envir) {
  if (before) {
    ## code to be run before a chunk
    #attach(loadNamespace("Crossover"), name="namespace:Crossover", pos=3)
    attach(loadNamespace("Crossover"), name="namespace:Crossover", pos=3, warn.conflicts=FALSE)
  } else {
    ## code to be run after a chunk
    detach("namespace:Crossover")
  }
}

if (exists("opts_chunk")) {
  knit_hooks$set(withNameSpace = CrossoverNamespace)
}

library(Crossover, quietly=TRUE)
options(width=80)
options(digits=4)
startGUI <- function(...) {invisible(NULL)}
#options(prompt="> ", continue="+ ")
library(MASS)
library(multcomp)
library(ggplot2)
library(Matrix)
# knitr has to be loaded for 'set_parent' and CRAN checks.
library(knitr)
bibCall <- TRUE


## ----williams3t, echo=TRUE----------------------------------------------------
getDesign("williams3t")

## ----general-carryover--------------------------------------------------------

design <- getDesign("williams3t")
general.carryover(design, model=1)


## ----models, echo=FALSE-------------------------------------------------------

cat(paste(1:9, ": \"", Crossover:::models, "\"", sep=""), sep="\n")


## ----pidgeon1, echo=TRUE------------------------------------------------------

getDesign("pidgeon1")


## ----pidgenVar, echo=FALSE----------------------------------------------------
design <- getDesign("pidgeon1")
design.efficiency(design)$var.trt.pair.adj


## ----echo=FALSE---------------------------------------------------------------
v <- length(levels(as.factor(design)))
im <- matrix(0, v, v)
vn <- sapply(1:v, function(x) {sum(design==x)})
for (i in 1:v) {
  for (j in 1:v) {
    if (i!=j) {
      im[i, j] <- 1/vn[i]+1/vn[j]
    }
  }
}
im


## ----pidgenEff, echo=FALSE----------------------------------------------------

design.efficiency(design)$eff.trt.pair.adj


## ----set-parent-models, echo=FALSE, cache=FALSE, include=FALSE----------------
set_parent('../Crossover.Rnw')
library(multcomp)
library(Crossover)

## ----StandardAdditive, echo=TRUE----------------------------------------------
# Design:
design <- rbind(c(3,2,1),
                c(2,1,3),
                c(1,2,3),
                c(3,2,1))
design
v <- 3 # number of treatments
# Link matrix:
H <- Crossover:::linkMatrix(model="Standard additive model", v)
H
# Row-Column-Design: (cf. John et al. 2004, Table II and page 2649f.)
rcDesign <- Crossover:::rcd(design, v=v, model=1)
rcDesign
# Design Matrix of Row-Column Design:
Xr <- Crossover:::rcdMatrix(rcDesign, v, model=1)
Xr
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----FullInteractions, echo=TRUE----------------------------------------------

H <- Crossover:::linkMatrix(model="Full set of interactions", v)
H
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----SelfAdjacency, echo=TRUE-------------------------------------------------

H <- Crossover:::linkMatrix(model="Self-adjacency model", v)
H
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----echo=TRUE, eval=TRUE-----------------------------------------------------
# Link matrix:
H <- Crossover:::linkMatrix(model="Placebo model", v, placebos=1)
H
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----NoIntoSelf, echo=TRUE----------------------------------------------------

H <- Crossover:::linkMatrix(model="No carry-over into self model", v)
H
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----TreatmentDecay, echo=TRUE------------------------------------------------

H <- Crossover:::linkMatrix(model="Treatment decay model", v)
H
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----Proportionality, echo=TRUE-----------------------------------------------

H <- Crossover:::linkMatrix(model="Proportionality model", v)
H
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----SecondOrder, echo=TRUE---------------------------------------------------
# Link matrix:
H <- Crossover:::linkMatrix(model="Second-order carry-over effects", v)
H
# Row-Column-Design: (cf. John et al. 2004, Table II and page 2649f.)
rcDesign <- Crossover:::rcd(design, v=v, model=8)
rcDesign
# Design Matrix of Row-Column Design:
Xr <- Crossover:::rcdMatrix(rcDesign, v, model=8)
Xr
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----bibtex-models, results='asis', echo=FALSE--------------------------------
if (!exists("bibCall")) {
  # RStudio / bibtex / knitr child document workaround from http://tex.stackexchange.com/questions/31373/citations-with-no-bibliography-at-the-end
  cat("\\newsavebox\\mytempbib \\savebox\\mytempbib{\\parbox{\\textwidth}{\\bibliography{../literatur}}}")
}

## ----set-parent-search, echo=FALSE, cache=FALSE, include=FALSE----------------
library(knitr)
set_parent('../Crossover.Rnw')
library(multcomp)
library(Crossover)
library(MASS)
library(xtable)

## ----SearchStrategy, echo=TRUE, eval=TRUE, cache=TRUE, dev='png', dpi=150-----

set.seed(42)
x <- searchCrossOverDesign(s=9, p=5, v=4, model=4)
plot(x)
plot(x, type=2)


## ----attachNameSpace, echo=TRUE, eval=FALSE, include=FALSE--------------------
#   # We will use a lot of internal commands since we will test and evaluate things the normal user will probably not be interested in. Therefore we load and attach the namespace.
#  # attach(loadNamespace("Crossover"), name="namespace:Crossover", pos=3)
#  # When we are finished we call 'detach("namespace:Crossover")'.

## ----TestOfDifferentApproaches, echo=TRUE, eval=TRUE, withNameSpace=TRUE------

attach(loadNamespace("Crossover"), name="namespace:Crossover", pos=3, warn.conflicts=FALSE)

s <- 6
p <- 3
v <- 3
model <- 1
design <- getDesign("williams3t")
  
rcDesign <- rcd(design, v, model)
# JRW, p 2650, first equation on that page, whithout number
Ar <- infMatrix(rcDesign, v, model)
Xr <- rcdMatrix(rcDesign, v, model)
# JRW, p 2650, second equation on that page, number 11
Ar2 <- t(Xr) %*% (diag(s*p)-Crossover:::getPZ(s,p)) %*% Xr
max(abs(Ar-Ar2))

# Testing the Projection of Z: P_Z times Z should equal Z:
max(abs(Crossover:::getPZ(s,p)%*%getZ(s,p)-getZ(s,p)))

fXr <- cbind(Xr, getZ(s,p))
Ar3 <- t(fXr) %*% fXr
max(abs(ginv(Ar3)[1:12,1:12]-ginv(Ar2)))

H <- linkMatrix(model=model,v=v)
fX <- cbind(Xr%*%H, getZ(s,p))
A1 <- t(fX) %*% fX
A2 <- t(H)%*%Ar%*%H

# While A1 and A2 of course differ (max(abs(A1[1:6,1:6]-A2))=2):

ginv(A1)[1:6,1:6]
ginv(A2)
max(abs(ginv(A1)[1:6,1:6]-ginv(A2)))

max(abs(ginv(A1, tol=10^-15)[1:6,1:6]-ginv(A2, tol=10^-15)))

# The variances for the estimable contrasts are the same:

C <- matrix(0,nrow=15,ncol=1)
C[1:2,1] <- c(-1,1)
tdiff1 <- t(C)%*%ginv(A1)%*%C
tdiff2 <- t(C[1:6,])%*%ginv(A2)%*%C[1:6,]
tdiff1 - tdiff2

C <- matrix(0,nrow=6,ncol=1)
C[1:2,1] <- c(-1,1)
tdiff1 <- t(C)%*%ginv(A1)[1:6,1:6]%*%C
tdiff2 <- t(C)%*%ginv(A2)%*%C
tdiff1 - tdiff2


## ----include=FALSE------------------------------------------------------------
data(exampleSearchResults2t)

## ----echo=FALSE, results='asis'-----------------------------------------------

options(xtable.include.rownames=FALSE, xtable.include.colnames=FALSE, xtable.floating=FALSE)

df <- c()

for (i in 1:length(resultL)) { 
  design <- resultL[[i]]
  cat("Design ",i,":\n")
  print(xtable(design, digits=0)) 
  var <- c()
  var2 <- c()
  for (m in models[1:8]) {
    
    #C <- Crossover:::contrMat2(type="Tukey", v=2, model=m, eff.factor=c(1,1,1))
    #C <- contrMat(rep(1, Crossover:::nrOfParameters(model=i, v=v)), "Tukey")
    
    if (Crossover:::estimable_R(design, v=2, model=m)) {
      var <- c(var, general.carryover(design, model=m)$Var.trt.pair[1,2]/4)      
    } else {
      var <- c(var, NA)      
    }
  }
  df <- rbind(df, var)  
}

rownames(df) <- NULL
df <- as.data.frame(df)
colnames(df) <- c("Additive", "Self-adjacency", "Proportional", "Placebo", "No into self", "Decay", "Interaction", "2nd-order carry-over")
df <- cbind(data.frame(Design=1:length(resultL)), df)

options(xtable.include.colnames=TRUE, xtable.NA.string="Not estimable")

cat("\\[\\]")

cat("\\scriptsize")

print(xtable(df, digits=3))

# cat("\\[\\]")
# print(xtable(df2, digits=3))
# max(abs(df-df2))


## ----echo=FALSE, include=FALSE, eval=FALSE------------------------------------
#  
#  TRIES <- 25
#  resultSubL <- list()
#  i <- 6
#  model <- models[i]
#  cat("======= ", model, " =======","\n")
#  for (k in 1:TRIES) {
#    result <- sortDesign(getDesign(searchCrossOverDesign(s=s, p=p, v=v, model=model, v.rep=c(12,12), eff.factor=c(1,0.01), contrast="Tukey"))) #, start.designs=list(designs[[i]]))
#    already.found <- FALSE
#    for (design in resultSubL) {
#      if (all(result==design)) already.found <- TRUE
#    }
#    if (!already.found) {
#      resultSubL <- c(resultSubL, list(result))
#      print(getDesign(result))
#      cat("\nTreatment: ", general.carryover(designs[[i]], model=i)$Var.trt.pair[1,2]/4, "(literature) vs. ", general.carryover(result, model=i)$Var.trt.pair[1,2]/4,"\n")
#      gco <- general.carryover(designs[[i]], model=i)
#      if (!is.null(dim(gco[[2]]))) {
#        cat("(1st) Carryover: ", gco[[2]][1,2]/4, "(literature) vs. ", general.carryover(result, model=i)[[2]][1,2]/4,"\n\n")
#      }
#      cat(paste("designB",i,sep=""), "<-", Crossover:::dputMatrix(getDesign(result)))
#    }
#  }
#  #resultL <- c(resultL, resultSubL)
#  #save(resultL, file="resultL.RData")
#  

## ----bibtex-search, results='asis', echo=FALSE--------------------------------
if (!exists("bibCall")) {
  # RStudio / bibtex / knitr child document workaround from http://tex.stackexchange.com/questions/31373/citations-with-no-bibliography-at-the-end
  cat("\\newsavebox\\mytempbib \\savebox\\mytempbib{\\parbox{\\textwidth}{\\bibliography{../literatur}}}")
}

## ----detachNameSpace, echo=TRUE, eval=FALSE, include=FALSE--------------------
#  detach("namespace:Crossover")

## ----set-parent-appendices, echo=FALSE, cache=FALSE, include=FALSE------------
set_parent('../Crossover.Rnw')
library(multcomp)
library(Crossover)

## ----mixed1, echo=TRUE, include=FALSE-----------------------------------------

tc <- textConnection("subject period treatment outcome
1 1 C 5.15
1 2 B 5.97
2 1 B 3.19
2 2 C 4.74
3 1 A 6.59
3 2 B 6.28
4 1 C 2.26
4 2 A 4.12
5 1 B 5.87
5 2 A 2.99
6 1 A 4.94
6 2 C 3.71
7 1 C 3.81
7 2 A 1.54
8 1 A 6.18
8 2 C 5.56
9 1 A 2.37
9 2 B 5.76
10 1 C 5.15
10 2 B 5.87
11 1 B 3.09
11 2 A 1.44
12 1 B 3.91
12 2 C 4.32
13 1 A 4.32
13 2 B 6.07
14 1 B 4.94
14 2 A 0.62
15 1 C 2.68
15 2 A 5.76
16 1 B 3.60
16 2 C 1.85
17 1 C 4.43
17 2 B 5.15
18 1 A 0.82
18 2 C 0.62")            

example5.2 <- read.table(tc, header = TRUE, as.is=TRUE)
close(tc)
example5.2$treatment <- factor(example5.2$treatment, levels=c("C", "A", "B"))

fit <- lm(outcome~subject+period+treatment, data=example5.2)
summary(glht(fit, linfct=mcp(treatment="Dunnett")))

library(nlme)

fit <- lme(outcome~period+treatment, random=~1| subject, data=example5.2)
summary(fit)
# Compare standard error with book (REML, p.219): check.
summary(glht(fit, linfct=mcp(treatment="Dunnett")))


## ----bibtex-appendices, results='asis', echo=FALSE----------------------------
if (!exists("bibCall")) {
  # RStudio / bibtex / knitr child document workaround from http://tex.stackexchange.com/questions/31373/citations-with-no-bibliography-at-the-end
  cat("\\newsavebox\\mytempbib \\savebox\\mytempbib{\\parbox{\\textwidth}{\\bibliography{../literatur}}}")
}


pkgname <- "gRc"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('gRc')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("add1drop1")
### * add1drop1

flush(stderr()); flush(stdout())

### Name: add1drop1
### Title: Add or drop colour classes to RCOX models
### Aliases: add1.rcox drop1.rcox
### Keywords: htest

### ** Examples

data(math)
gc.sat <- ~me:ve:al:st:an
gc.1   <- ~me+ve+al+st+an

m.sat <- rcox(gm=gc.sat, data=math)
m.1   <- rcox(gm=gc.1,   data=math)

t.sat <- drop1(m.sat)
t.sat$tab
t.sat$cc

t.1   <- add1(m.1)
t.1$tab
t.1$cc



cleanEx()
nameEx("comparecc")
### * comparecc

flush(stderr()); flush(stdout())

### Name: comparecc
### Title: Compare colour classes of an RCOX model
### Aliases: comparecc
### Keywords: htest

### ** Examples


data (math)

gm  = ~al:an:st
vcc = list(~me+st, ~ve+an, ~al)
ecc = list(~me:ve+me:al, ~ve:al+al:st)

m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math)
m1

comparecc(m1, type="vcc")
comparecc(m1, type="ecc")



cleanEx()
nameEx("fit")
### * fit

flush(stderr()); flush(stdout())

### Name: fit
### Title: Fit RCOX models
### Aliases: fit fit.rcox
### Keywords: models

### ** Examples


data(math)
gm  = ~al:an:st
vcc = list(~me+st, ~ve+an, ~al)
ecc = list(~me:ve+me:al, ~ve:al+al:st)

m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, fit=FALSE)

fit(m1, method="matching")
fit(m1, method="scoring")
fit(m1, method="ipm")

## MISSING



cleanEx()
nameEx("fitInternal")
### * fitInternal

flush(stderr()); flush(stdout())

### Name: fitInteral
### Title: Functions used in connection with fitting of RCOX models
### Aliases: scoring scoring.rcox ipm ipm.rcon ipm.rcor matching
###   matching.rcon matching.rcor rconIPM rcorIPM rconScoreMatch
###   rcorScoreMatch rconScoreTheta rcorScoreTheta fitIPSedge fitIPSset
###   fitNR2 fitNR modNewt refitA
### Keywords: internal

### ** Examples

gm  = ~al:an:st
vcc = list(~me+st, ~ve+an, ~al)
ecc = list(~me:ve+me:al, ~ve:al+al:st)
data(math)

m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, fit=FALSE)
f1 <- matching(m1)

## Use f1$K as starting value
scoring(m1, K0=f1$K)
ipm(m1, K0=f1$K)




cleanEx()
nameEx("getSlot")
### * getSlot

flush(stderr()); flush(stdout())

### Name: getSlot
### Title: Accessing RCOX model objects
### Aliases: getSlot getSlot fitInfo dataRep intRep getedges getecc getvcc
### Keywords: utilities

### ** Examples

data(math)
gm  = ~al:an:st
vcc = list(~me+st, ~ve+an, ~al)
ecc = list(~me:ve+me:al, ~ve:al+al:st)

m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math)
getecc(m1)

getSlot(m1,"type")
fitInfo(m1)
fitInfo(m1,"K")




cleanEx()
nameEx("join1split1")
### * join1split1

flush(stderr()); flush(stdout())

### Name: join1split1
### Title: Joining and splitting of colour classes in RCOX models
### Aliases: join1 split1
### Keywords: htest

### ** Examples

data(math)
g1     <- ~me:ve:al+al:st:an
m1     <- rcox(gm=g1, data=math)
join1(m1)

gm  = ~al:an:st
vcc = list(~me+st, ~ve+an)
ecc = list(~me:ve+me:al, ~ve:al+al:st)
m2 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, type="rcon")
split1(m2)



cleanEx()
nameEx("rcox")
### * rcox

flush(stderr()); flush(stdout())

### Name: rcox
### Title: Main function for specifying RCON/RCOR models
### Aliases: rcox
### Keywords: models

### ** Examples


data(math)
gm  = ~al:an:st
vcc = list(~me+st, ~ve+an, ~al)
ecc = list(~me:ve+me:al, ~ve:al+al:st)

m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, method='matching')
m2 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, method='scoring')
m3 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, method='ipm')

m1
m2
m3

summary(m1)
summary(m2)
summary(m3)

coef(m1)
coef(m2)
coef(m3)

vcov(m1)
vcov(m2)
vcov(m3)



cleanEx()
nameEx("tr")
### * tr

flush(stderr()); flush(stdout())

### Name: tr
### Title: Calculate trace of various matrix products
### Aliases: trA trAW trAWB trAWBW trAWBV
### Keywords: utilities

### ** Examples

d <- 5
W <- matrix(rnorm(d*d),nr=d,nc=d); 
V <- W <- W+t(W)

## Turn list into matrix
##
tomat <- function(x){
  ans <- do.call("rbind", x)
  storage.mode(ans)<-"double"
  return(ans)
}

A1 <- tomat(list(c(1,2),c(1,3)))
A2 <- tomat(list(1,3,5))


## Just for checking the calculations
##
symMat <- function(A,d){
  ans <- matrix(0,nr=d,nc=d)
  for (i in 1:length(A)){
    e <- A[[i]]
    if (length(e)==1){
      ans[e,e] <- 1
    } else { 
      ans[e[1],e[2]] <-   ans[e[2],e[1]] <- 1 
    }
  }

  return(ans)
}

trAW(A1, W)
#sum(diag(symMat(A1,d=d) %*% W))

trAW(A2, W)
#sum(diag(symMat(A2,d=d) %*% W))

trAWB(A1, W, A2)
#sum(diag(symMat(A1,d=d) %*% W %*% symMat(A2,d=d)))

trAWBV(A1, W, A2, V)
#sum(diag(symMat(A1,d=d) %*% W %*% symMat(A2,d=d) %*% V))



cleanEx()
nameEx("update")
### * update

flush(stderr()); flush(stdout())

### Name: update.rcox
### Title: Update an RCOX model
### Aliases: update.rcox
### Keywords: models

### ** Examples


data(math)
gm  = ~al:an:st
vcc = list(~me+st, ~ve+an, ~al)
ecc = list(~me:ve+me:al, ~ve:al+al:st)

m1 <- rcox(gm=gm, vcc=vcc, ecc=ecc, data=math, method='matching', trace=0)

update(m1, joinvcc=list(~me+st, ~ve+an))
update(m1, joinecc=list(~al:an, ~an:st))

update(m1, splitvcc=~ve+an)
update(m1, splitecc=~me:ve+me:al)


update(m1, dropecc=list(~me:st+st:an,~al:an,~st:al))
update(m1, addecc=list(~an:me+st:ve))





### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

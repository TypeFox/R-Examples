### R code from vignette source 'gRbase-arrayops.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: gRbase-arrayops.Rnw:31-34
###################################################
require( gRbase )
prettyVersion <- packageDescription("gRbase")$Version
prettyDate <- format(Sys.Date())


###################################################
### code chunk number 2: gRbase-arrayops.Rnw:80-83
###################################################
##library(gRbase)
options("width"=100)
options(useFancyQuotes="UTF-8")


###################################################
### code chunk number 3: gRbase-arrayops.Rnw:114-131
###################################################
## 1-dimensional array
x1 <- 1:8
dim(x1) <- 8
x1
c(is.array(x1), is.matrix(x1))

## 2-dimensional array (matrix)
x2 <- 1:8
dim(x2) <- c(2,4)
x2
c(is.array(x2), is.matrix(x2))

## 3-dimensional array
x3 <- 1:8
dim(x3) <- c(2,2,2)
x3
c(is.array(x3), is.matrix(x3))


###################################################
### code chunk number 4: gRbase-arrayops.Rnw:143-145
###################################################
data( lizard, package="gRbase" )
lizard


###################################################
### code chunk number 5: gRbase-arrayops.Rnw:151-153
###################################################
dim( lizard )
dimnames(lizard)


###################################################
### code chunk number 6: gRbase-arrayops.Rnw:160-170
###################################################
## because lizard is a vector
lizard[1:2]
## because lizard is an array
lizard[,1,1]
## useful if subsetting programatically
margin <- 2:3; level <- c(1,1)
z <- as.list(rep(TRUE,3))
z[margin] <- level
str(z)
do.call("[", c(list(lizard), z))


###################################################
### code chunk number 7: gRbase-arrayops.Rnw:174-177
###################################################
is.array( lizard[1:2] )
is.array( lizard[,1,1] )
is.array( lizard[,1,] )


###################################################
### code chunk number 8: gRbase-arrayops.Rnw:192-193
###################################################
T.DHS <- lizard


###################################################
### code chunk number 9: gRbase-arrayops.Rnw:206-210
###################################################
T.DS <- tabMarg(lizard, ~diam+species); T.DS
## Alternative forms
T.DS <- tabMarg(lizard, c("diam","species"))
T.DS <- tabMarg(lizard, c(1,3))


###################################################
### code chunk number 10: gRbase-arrayops.Rnw:214-216
###################################################
T.DS <- tableMargin(lizard, ~diam+species); T.DS
T.HS <- tableMargin(lizard, ~height+species); T.HS


###################################################
### code chunk number 11: gRbase-arrayops.Rnw:225-230
###################################################
T.SHD <- tabPerm(T.DHS, ~species+height+diam); ftable( T.SHD )
## Alternative forms:
T.SHD <- tabPerm(T.DHS, c("species","height","diam"))
T.SHD <- tabPerm(T.DHS, c(3,2,1))
T.SHD <- aperm(T.DHS, c(3,2,1))


###################################################
### code chunk number 12: gRbase-arrayops.Rnw:240-241
###################################################
tabSlice(lizard, slice=list(species="anoli"))


###################################################
### code chunk number 13: gRbase-arrayops.Rnw:255-258
###################################################
T.HS <- tabMarg(lizard, ~height+species)
T.DHS.mult = tabMult( T.DS, T.HS )
ftable( T.DHS.mult )


###################################################
### code chunk number 14: gRbase-arrayops.Rnw:270-273
###################################################
T.DHS.div  = tabDiv( T.DS, T.HS )
T.DHS.add  = tabAdd( T.DS, T.HS )
T.DHS.subt = tabSubt( T.DS, T.HS )


###################################################
### code chunk number 15: gRbase-arrayops.Rnw:288-290
###################################################
T.HSD2 <- tabExpand(T.DS, T.HS); T.HSD2
names(dimnames(T.HSD2))


###################################################
### code chunk number 16: gRbase-arrayops.Rnw:302-303
###################################################
tabAlign(T.SHD, T.DHS)


###################################################
### code chunk number 17: gRbase-arrayops.Rnw:308-315
###################################################
n.new <- names(dimnames(T.DHS)); n.new
n.old <- names(dimnames(T.SHD)); n.old
if (setequal( n.new, n.old )){
    tabPerm( T.SHD, n.new )
} else {
    numeric(0)
}


###################################################
### code chunk number 18: gRbase-arrayops.Rnw:321-322
###################################################
tabEqual( T.SHD, T.DHS )


###################################################
### code chunk number 19: gRbase-arrayops.Rnw:340-343
###################################################
yn <- c("y","n")
T.AB <- array(c(5,95,1,99), dim=c(2,2), dimnames=list("A"=yn, "B"=yn))
T.AB <- parray(c("A","B"), levels=list(yn, yn), values=c(5,95,1,99))


###################################################
### code chunk number 20: gRbase-arrayops.Rnw:352-356
###################################################
T.AB <- parray(c("A","B"), levels=list(yn, yn), values=c(5,95,1,99),
               normalize="first")
T.AB <- parray(c("A","B"), levels=list(yn, yn), values=c(5,95,1,99),
               normalize="all")


###################################################
### code chunk number 21: gRbase-arrayops.Rnw:377-385
###################################################
r <- parray("rain", levels=list(yn), values=c(.2, .8))
s.r <- parray(c("sprinkler","rain"), levels=list(yn,yn),
              values=c(.01, .99, .4, .6))
w.sr <- parray(c("wet","sprinkler","rain"), list(yn,yn,yn),
               values=c(.99, .01, .8, .2, .9, .1, 0, 1))
r
s.r
ftable(w.sr, col.vars = "wet")


###################################################
### code chunk number 22: gRbase-arrayops.Rnw:393-397
###################################################
joint <- tabMult( tabMult(r, s.r), w.sr )
ftable(joint)
## Alternative
joint <- tabListMult( list( r, s.r, w.sr ) )


###################################################
### code chunk number 23: gRbase-arrayops.Rnw:403-407
###################################################
wr.marg <- tabMarg(joint, ~wet+rain); wr.marg
tabDiv( wr.marg, tabMarg(wr.marg, ~wet))
## Alternative -- and shorter
tabCondProb(wr.marg, cond=~wet)


###################################################
### code chunk number 24: gRbase-arrayops.Rnw:412-415
###################################################
x <- tabSliceMult(wr.marg, slice=list(wet="y")); x
tabMarg(x, ~rain)
tabCondProb( tabMarg(x, ~rain) )


###################################################
### code chunk number 25: gRbase-arrayops.Rnw:449-469
###################################################
myips <- function(indata, glist){
    fit   <- indata
    fit[] <-  1
    ## List of sufficient marginal tables
    md    <- lapply(glist, function(g) tabMarg(indata, g))

    for (i in 1:4){
        for (j in seq_along(glist)){
            mf  <- tabMarg(fit, glist[[j]])
            adj <- tabDiv( md[[j]], mf)
            fit <- tabMult( fit, adj )
        }
    }
    pearson=sum( (fit-indata)^2 / fit)
    pearson
}

glist<-list(c("species","diam"),c("species","height"),c("diam","height"))
str( myips(lizard, glist), max.level=2)
str( loglin(lizard, glist), max.level = 2)


###################################################
### code chunk number 26: gRbase-arrayops.Rnw:483-485
###################################################
dim2222 <- c(2,2,2,2)
dim2323 <- c(2,3,2,3)


###################################################
### code chunk number 27: gRbase-arrayops.Rnw:497-501
###################################################
cell2entry(c(1,1,1,1), dim2222)
entry2cell(1, dim2222)
cell2entry(c(2,1,2,1), dim2222)
entry2cell(6, dim2222)


###################################################
### code chunk number 28: gRbase-arrayops.Rnw:526-528
###################################################
nextCell(c(1,1,2,1), dim2222)
nextCell(c(2,2,2,1), dim2222)


###################################################
### code chunk number 29: gRbase-arrayops.Rnw:539-541
###################################################
nextCellSlice(c(2,1,1,2),  sliceset=2, dim2323)
nextCellSlice(c(1,3,2,1),  sliceset=c(2,3), dim2323)


###################################################
### code chunk number 30: gRbase-arrayops.Rnw:559-560
###################################################
r1<-slice2entry(slicecell=c(1,2), sliceset=c(2,3), dim2222); r1


###################################################
### code chunk number 31: gRbase-arrayops.Rnw:566-567
###################################################
do.call(rbind, lapply(r1, entry2cell, dim2222))


###################################################
### code chunk number 32: gRbase-arrayops.Rnw:579-580
###################################################
p<-permuteCellEntries(perm=c(2,1), dim=c(2,3)); p


###################################################
### code chunk number 33: gRbase-arrayops.Rnw:586-590
###################################################
(A <- array(11:16, dim=c(2,3)))
Ap <- A[p]
dim(Ap) <- c(3,2)
Ap


###################################################
### code chunk number 34: gRbase-arrayops.Rnw:596-597
###################################################
aperm(A, c(2,1))


###################################################
### code chunk number 35: gRbase-arrayops.Rnw:614-617
###################################################
ff <- factGrid(dim2222)
head(ff, 4)
tail(ff, 4)


###################################################
### code chunk number 36: gRbase-arrayops.Rnw:623-625
###################################################
aa <- expand.grid(list(1:2,1:2,1:2,1:2))
head(aa, 4)


###################################################
### code chunk number 37: gRbase-arrayops.Rnw:630-631
###################################################
factGrid(dim2222, slicecell=c(1,2), sliceset=c(2,3))



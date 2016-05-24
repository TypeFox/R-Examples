### R code from vignette source 'planning.rnw'

###################################################
### code chunk number 1: Setup
###################################################
rm( list=ls() )

##source( "../../scripts/schedpack.r" )
##source( "../../scripts/schedpack.r" )

source( "../../scripts/from-weikko-jarros-schedpack.r" )

## this is the switch for the dan pickett problem, which can take considerable time
## this is the huge 'middle' matrix that contains the parameters for the 'model'
include.primary.dcm <- FALSE
include.primary.dcm <- TRUE

if(!require(xtable, quietly=TRUE)) install.packages("xtable")

##if(!require(glpk, quietly=TRUE)) 
##  install.packages("glpk", repos="http://cran.r-project.org")

##if(!require(e1071, quietly=TRUE)) 
##  install.packages("e1071", repos="http://cran.r-project.org", dependencies=TRUE)

##if(!require(maptools, quietly=TRUE)) install.packages("maptools")

if(!require(lattice, quietly=TRUE)) install.packages("lattice")

## you have to make sure this works as intended.
if(!require(glpkAPI, quietly=TRUE)) install.packages("glpkAPI", "--configure-args='--enable-gmp=no'")

options(width=65)


###################################################
### code chunk number 2: planning.rnw:415-420
###################################################

A <- 8     # A
T <- 6     # T
l <- 10    # \ell
P <- T * l # |T|


###################################################
### code chunk number 3: planning.rnw:423-424
###################################################
stand.acres <- c(5000,5000,5000,5000,30000,20000,12000,2000)


###################################################
### code chunk number 4: planning.rnw:434-436
###################################################
leusch.ylds <- read.table("../../data/leuschner.txt", 
                          header = TRUE)


###################################################
### code chunk number 5: planning.rnw:458-462
###################################################
library(glpkAPI)
leusch.lp <- initProbGLPK()
##library(glpk)
##leusch.lp <- lpx_create_prob()


###################################################
### code chunk number 6: planning.rnw:467-469
###################################################
##lpx_set_prob_name(leusch.lp, "leuschner")
setProbNameGLPK(leusch.lp, pname = "leuschner")


###################################################
### code chunk number 7: planning.rnw:477-479
###################################################
setObjDirGLPK(leusch.lp, GLP_MAX)
##lpx_set_obj_dir(leusch.lp, LPX_MAX)


###################################################
### code chunk number 8: planning.rnw:538-540
###################################################
##lpx_add_cols(leusch.lp, nrow( leusch.ylds ))
addColsGLPK(leusch.lp, nrow( leusch.ylds ))


###################################################
### code chunk number 9: planning.rnw:556-562
###################################################

leusch.ylds$dv <- paste("a", 
                        leusch.ylds$stand, 
                        leusch.ylds$period, 
                        sep="")
head(leusch.ylds)


###################################################
### code chunk number 10: planning.rnw:591-596
###################################################
for( t in 1:nrow(leusch.ylds)) {
  label <- leusch.ylds[t,]$dv
  ##lpx_set_col_name(leusch.lp, t, label)
  setColNameGLPK(leusch.lp, t, label)
}


###################################################
### code chunk number 11: planning.rnw:607-611
###################################################
for(t in 1:nrow(leusch.ylds)) {
##  lpx_set_col_bnds(leusch.lp, t, LPX_LO, 0.0, 0.0)
  setColBndGLPK(leusch.lp, t, GLP_LO, 0.0, 0.0) 
}


###################################################
### code chunk number 12: planning.rnw:629-634
###################################################
for(t in 1:nrow( leusch.ylds)) {
  ## set the objective coefficient in period t = 0.0
##  lpx_set_obj_coef(leusch.lp, t, leusch.ylds[t,]$vol)
  setObjCoefsGLPK(leusch.lp, t, leusch.ylds[t,]$vol)
}


###################################################
### code chunk number 13: planning.rnw:717-719
###################################################
##lpx_add_rows(leusch.lp, length(stand.acres))
addRowsGLPK(leusch.lp, length(stand.acres))


###################################################
### code chunk number 14: planning.rnw:727-730
###################################################
acre.consts <- t( kronecker(diag(length(stand.acres)), 
                            as.matrix( rep(1,T))))



###################################################
### code chunk number 15: planning.rnw:749-750
###################################################
val <- rep(1, T)  # this is the value of the constraint


###################################################
### code chunk number 16: planning.rnw:760-795
###################################################

## loop over the stands to generate 
## the idx is the index for the set of rows
## we are trying to set the coefficients for
for(i in 1:length(stand.acres)) {

  idx <- rep(0, T)  # this is the index of the col coeffs

  ## manually create the index using id
  ## and assign row i, column j stored as idx[id]
  id <- 1
  for(j in 1:ncol(acre.consts)) {
    if(acre.consts[i,j] != 0.0) {
      idx[id] <- j
      id <- id + 1
    }
  }  
 
  ## now, perform all three tasks within the same loop
  ## set the constraint row name
  ##lpx_set_row_name(leusch.lp, i, paste("s", i, sep=""))
  setRowNameGLPK(leusch.lp, i, paste("s", i, sep=""))

  ## set the upper bound on the acres for this stand
  ## to make sure no more than stand.acres[i] are cut
  ##lpx_set_row_bnds(leusch.lp, i, LPX_UP, 0.0, 
  ##                 stand.acres[i])
  setRowBndGLPK(leusch.lp, i, GLP_UP, 0.0, 
                     stand.acres[i])
  
  ## set matrix row coefficients (?glp_set_mat_row)
  ## idx = ???, row index = i, val = vector of 1's
  ##lpx_set_mat_row(leusch.lp, i, length(val), idx, val)  
  setMatRowGLPK(leusch.lp, i, length(val), idx, val)
}


###################################################
### code chunk number 17: planning.rnw:803-807
###################################################
## lpx_get_num_rows(leusch.lp)  
## lpx_get_num_cols(leusch.lp)
getNumRowsGLPK(leusch.lp)
getNumColsGLPK(leusch.lp)


###################################################
### code chunk number 18: planning.rnw:854-857
###################################################

target.acres <- rep(14000, T)



###################################################
### code chunk number 19: planning.rnw:866-868
###################################################
##lpx_add_rows(leusch.lp, length(target.acres))
addRowsGLPK(leusch.lp, length(target.acres))


###################################################
### code chunk number 20: planning.rnw:877-910
###################################################

val <- rep(1,8) ## a vector of 8 ones.

for(i in 1:length(target.acres)) {
  idx <- 
     as.numeric(rownames(subset(leusch.ylds, period == i)))
##  lpx_set_row_name(leusch.lp, 
##                   i + length(stand.acres), 
##                   paste("tac", i, sep="" ))
  setRowNameGLPK(leusch.lp,
                   i + length(stand.acres), 
                   paste("tac", i, sep="" ))
  ## lpx_set_row_bnds(leusch.lp, 
  ##                  i + length(stand.acres), 
  ##                  LPX_FX, 
  ##                  target.acres[i], 
  ##                  target.acres[i])
  setRowBndGLPK(leusch.lp,
                   i + length(stand.acres), 
                   GLP_FX,
                   target.acres[i], 
                   target.acres[i])
  ## lpx_set_mat_row(leusch.lp, 
  ##                 i + length(stand.acres), 
  ##                 length(val), 
  ##                 idx, 
  ##                 val)
  setMatRowGLPK(leusch.lp,
                  i + length(stand.acres), 
                  length(val), 
                  idx, 
                  val)
}


###################################################
### code chunk number 21: planning.rnw:926-930
###################################################
##lpx_get_num_rows(leusch.lp)  
##lpx_get_num_cols(leusch.lp)
getNumRowsGLPK(leusch.lp)
getNumColsGLPK(leusch.lp)


###################################################
### code chunk number 22: planning.rnw:953-955
###################################################
##lpx_simplex(leusch.lp)
solveSimplexGLPK(leusch.lp)


###################################################
### code chunk number 23: planning.rnw:1008-1011
###################################################
##leusch.obj <- lpx_get_obj_val(leusch.lp)
leusch.obj <- getObjValGLPK(leusch.lp)
leusch.obj 


###################################################
### code chunk number 24: planning.rnw:1038-1039
###################################################
leusch.col.rpt <- get.col.report(leusch.lp)


###################################################
### code chunk number 25: planning.rnw:1066-1070
###################################################
ac.per <- 
   matrix(as.numeric(as.matrix(I(leusch.col.rpt$activity))), 
          T, 8)
ac.per


###################################################
### code chunk number 26: planning.rnw:1082-1083
###################################################
rowSums(ac.per)


###################################################
### code chunk number 27: planning.rnw:1101-1103
###################################################
vol.cut.a <- c(199.4,197.4,203.4,211.4,217.6,195.4)



###################################################
### code chunk number 28: leuschner-area-control
###################################################

##vol.cut.a <- c(199.400,197.400,203.400,211.400,217.600,195.400)
vol.cut.b <- c(203.79, 203.79, 203.79, 203.79, 203.79, 203.79)

par(las = 1)
plot(vol.cut.a,  ylim=c(190,220), type="b", lwd=2, lty=1,
     main="Woodflow",
     xlab="Cutting period", ylab="Woodflow (cf x 10^3)")
par( new=TRUE )
plot( vol.cut.b, ylim=c(190,220), type="l", lwd=2, 
     lty=1, pch=3, main=NA, xlab=NA, ylab=NA )

##abline( h=mean( vol.cut.a ) - 0.05 * mean( vol.cut.a ), lty=2 )
##abline( h=mean( vol.cut.a ) + 0.05 * mean( vol.cut.a ), lty=2 )

leg <- cbind( problem=c("Strict Area Control","Strict Volume Control"), 
             color=c(1,1) )
leg.text <- leg[,1]
leg.col <- leg[,2]
legend(1, 220, 
       leg.text, 
       ##fill=leg.col,
       cex=1.0, 
       col=leg.col, 
       ##pt.bg=leg.col,
       lwd=c(2,2),
       lty=c(1,1),
       pch=c(1,NA),
       ncol=1,
       title="Plan Types" )



###################################################
### code chunk number 29: fig-leuschner-area-control
###################################################

##vol.cut.a <- c(199.400,197.400,203.400,211.400,217.600,195.400)
vol.cut.b <- c(203.79, 203.79, 203.79, 203.79, 203.79, 203.79)

par(las = 1)
plot(vol.cut.a,  ylim=c(190,220), type="b", lwd=2, lty=1,
     main="Woodflow",
     xlab="Cutting period", ylab="Woodflow (cf x 10^3)")
par( new=TRUE )
plot( vol.cut.b, ylim=c(190,220), type="l", lwd=2, 
     lty=1, pch=3, main=NA, xlab=NA, ylab=NA )

##abline( h=mean( vol.cut.a ) - 0.05 * mean( vol.cut.a ), lty=2 )
##abline( h=mean( vol.cut.a ) + 0.05 * mean( vol.cut.a ), lty=2 )

leg <- cbind( problem=c("Strict Area Control","Strict Volume Control"), 
             color=c(1,1) )
leg.text <- leg[,1]
leg.col <- leg[,2]
legend(1, 220, 
       leg.text, 
       ##fill=leg.col,
       cex=1.0, 
       col=leg.col, 
       ##pt.bg=leg.col,
       lwd=c(2,2),
       lty=c(1,1),
       pch=c(1,NA),
       ncol=1,
       title="Plan Types" )



###################################################
### code chunk number 30: planning.rnw:1184-1185
###################################################
leusch.col.rpt[1,] 


###################################################
### code chunk number 31: planning.rnw:1212-1215
###################################################
leusch.row.rpt <- get.row.report(leusch.lp)
leusch.row.rpt



###################################################
### code chunk number 32: planning.rnw:1248-1251
###################################################

leusch.row.rpt[1,]



###################################################
### code chunk number 33: planning.rnw:1303-1305
###################################################
##lpx_write_mps(leusch.lp, "leuschner.mps")
writeMPSGLPK(leusch.lp, GLP_MPS_FILE, "leuschner.mps")


###################################################
### code chunk number 34: write leush
###################################################
##lpx_write_cpxlp(leusch.lp, "leuschner.xlp")
writeLPGLPK(leusch.lp, "leuschner.xlp")


###################################################
### code chunk number 35: planning.rnw:1387-1390
###################################################
##lpx_delete_prob(leusch.lp)
delProbGLPK(leusch.lp)
leusch.lp


###################################################
### code chunk number 36: planning.rnw:1448-1450
###################################################
system("rm -fr package-Ch9")
package.skeleton(name = "package-Ch9")



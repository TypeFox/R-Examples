###################################################
### chunk number 1: Setup
###################################################
rm( list=ls() )

source( "../../scripts/schedpack.r" )

## this is the switch for the dan pickett problem, which can take considerable time
## this is the huge 'middle' matrix that contains the parameters for the 'model'
include.primary.dcm <- FALSE
include.primary.dcm <- TRUE

if(!require(xtable, quietly=TRUE)) install.packages("xtable")
if(!require(glpkAPI, quietly=TRUE)) 
  install.packages("glpkAPI", repos="http://cran.r-project.org")
if(!require(e1071, quietly=TRUE)) 
  install.packages("e1071", repos="http://cran.r-project.org", dependencies=TRUE)
if(!require(maptools, quietly=TRUE)) install.packages("maptools")
if(!require(lattice, quietly=TRUE)) install.packages("lattice")

options(width=65)


###################################################
### chunk number 2: 
###################################################

A <- 8     # A
T <- 6     # T
l <- 10    # \ell
P <- T * l # |T|



###################################################
### chunk number 3: 
###################################################

## \lambda(\mathcal{A}_{b})
stand.acres <- c(5000,5000,5000,5000,30000,20000,12000,2000)



###################################################
### chunk number 4: 
###################################################


leusch.ylds <- read.table("../../data/leuschner.txt", 
                          header = TRUE)



###################################################
### chunk number 5: 
###################################################

#library(glpk)
#leusch.lp <- lpx_create_prob()


#call glp_create_prob using glpkAPI initProbGLPK(ptrtype = "glpk_prob")

library(glpkAPI)
leusch.lp <- initProbGLPK()


###################################################
### chunk number 6: 
###################################################
#lpx_set_prob_name(leusch.lp, "leuschner")

#call glp_set_prob_name|assign (change) objective function name
#using GLPKAPI setProbNameGLPK(lp, pname = NULL)

setProbNameGLPK(leusch.lp, pname = "leushcner")


###################################################
### chunk number 7: 
###################################################

#lpx_set_obj_dir(leusch.lp, LPX_MAX)

#call glp_set_obj_dir
#Usage setObjDirGLPK(lp, lpdir)

setObjDirGLPK(leusch.lp, GLP_MAX)

###################################################
### chunk number 8: 
###################################################

#lpx_add_cols(leusch.lp, nrow( leusch.ylds ))

#Call GLPK function glp_add_cols. 
#Usage addColsGLPK(lp, ncols)

addColsGLPK(leusch.lp, nrow( leusch.ylds ))

###################################################
### chunk number 9: 
###################################################

leusch.ylds$dv <- paste("a", 
                        leusch.ylds$stand, 
                        leusch.ylds$period, 
                        sep="")
head(leusch.ylds)



###################################################
### chunk number 10: 
###################################################

for( t in 1:nrow(leusch.ylds)) {
  label <- leusch.ylds[t,]$dv
#  lpx_set_col_name(leusch.lp, t, label)
  setColNameGLPK(leusch.lp, t, label)
  }


#Call glp_set_col_name. 
#Usage setColNameGLPK(lp, j, cname = NULL)


###################################################
### chunk number 11: 
###################################################
for(t in 1:nrow(leusch.ylds)) {
#  lpx_set_col_bnds(leusch.lp, t, LPX_LO, 0.0, 0.0)
  setColBndGLPK(leusch.lp, t, GLP_LO, 0.0, 0.0) 
}

#Call GLPK function glp_set_col_bnds. 
#Usage setColBndGLPK(lp, j, type, lb, ub)


###################################################
### chunk number 12: 
###################################################

for(t in 1:nrow( leusch.ylds)) {
  ## set the objective coefficient in period t = 0.0
  #lpx_set_obj_coef(leusch.lp, t, leusch.ylds[t,]$vol)
  setObjCoefsGLPK(leusch.lp, t, leusch.ylds[t,]$vol)
} #$

#calls the GLPK function glp_set_obj_coef
#Usage setObjCoefsGLPK(lp, j, obj_coef)

###################################################
### chunk number 13: 
###################################################
#lpx_add_rows(leusch.lp, length(stand.acres))

addRowsGLPK(leusch.lp, length(stand.acres))

#Calls the GLPK function glp_add_rows. 
#Usage addRowsGLPK(lp, nrows)

###################################################
### chunk number 14: 
###################################################
acre.consts <- t( kronecker(diag(length(stand.acres)), 
                            as.matrix( rep(1,6))))



###################################################
### chunk number 15: 
###################################################

val <- rep(1, 6)  # this is the value of the constraint



###################################################
### chunk number 16: 
###################################################

## loop over the stands to generate 
## the idx is the index for the set of rows
## we are trying to set the coefficients for
for(i in 1:length(stand.acres)) {

  idx <- rep(0, 6)  # this is the index of the col coeffs

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
#  lpx_set_row_name(leusch.lp, i, paste("s", i, sep=""))
  
# Call GLPK function glp_set_row_name. 
#Usage setRowNameGLPK(lp, i, rname = NULL)
  setRowNameGLPK(leusch.lp, i, paste("s", i, sep=""))
  
  ## set the upper bound on the acres for this stand
  ## to make sure no more than stand.acres[i] are cut
  
#  lpx_set_row_bnds(leusch.lp, i, LPX_UP, 0.0, 
#                   stand.acres[i])
  # Call the GLPK function glp_set_row_bnds. 
  # Usage setRowBndGLPK(lp, i, type, lb, ub)
  
  setRowBndGLPK(leusch.lp, i, GLP_UP, 0.0, 
                     stand.acres[i])
  
  ## set matrix row coefficients (?glp_set_mat_row)
  ## idx = ???, row index = i, val = vector of 1's
#  lpx_set_mat_row(leusch.lp, i, length(val), idx, val) 
  
# Call GLPK function glp_set_mat_row. 
# Usage setMatRowGLPK(lp, i, len, ind, val)

  setMatRowGLPK(leusch.lp, i, length(val), idx, val)
}



###################################################
### chunk number 17: 
###################################################

# Call the GLPK function glp_get_num_rows. 
# Usage getNumRowsGLPK(lp)

#lpx_get_num_rows(leusch.lp)  
getNumRowsGLPK(leusch.lp)


# Call the GLPK function glp_get_num_rows. 
# Usage getNumColsGLPK(lp)

#lpx_get_num_cols(leusch.lp)
getNumColsGLPK(leusch.lp)


###################################################
### chunk number 18: 
###################################################

target.acres <- rep(14000, 6)



###################################################
### chunk number 19: 
###################################################


#Calls the GLPK function glp_add_rows. 
#Usage addRowsGLPK(lp, nrows)

#lpx_add_rows(leusch.lp, length(target.acres))
addRowsGLPK(leusch.lp, length(target.acres))


###################################################
### chunk number 20: 
###################################################

# Call GLPK function glp_set_row_name. 
#Usage setRowNameGLPK(lp, i, rname = NULL)

# Call the GLPK function glp_set_row_bnds. 
# Usage setRowBndGLPK(lp, i, type, lb, ub)

# Call GLPK function glp_set_mat_row. 
# Usage setMatRowGLPK(lp, i, len, ind, val)


val <- rep(1,8) ## a vector of 8 ones.

for(i in 1:length(target.acres)) {
  idx <- 
     as.numeric(rownames(subset(leusch.ylds, period == i)))
#  lpx_set_row_name(leusch.lp, 
  setRowNameGLPK(leusch.lp,
                   i + length(stand.acres), 
                   paste("tac", i, sep="" ))
#  lpx_set_row_bnds(leusch.lp, 
  setRowBndGLPK(leusch.lp,
                   i + length(stand.acres), 
#                   LPX_FX, 
                   GLP_FX,
                   target.acres[i], 
                   target.acres[i])
#  lpx_set_mat_row(leusch.lp, 
  setMatRowGLPK(leusch.lp,
                  i + length(stand.acres), 
                  length(val), 
                  idx, 
                  val)
}



###################################################
### chunk number 21: 
###################################################

# Call the GLPK function glp_get_num_rows. 
# Usage getNumRowsGLPK(lp)

#lpx_get_num_rows(leusch.lp)  
getNumRowsGLPK(leusch.lp)


# Call the GLPK function glp_get_num_rows. 
# Usage getNumColsGLPK(lp)

#lpx_get_num_cols(leusch.lp)
getNumColsGLPK(leusch.lp)




###################################################
### chunk number 22: 
###################################################

# Call the  GLPK function glp_simplex. 
# Usage solveSimplexGLPK(lp)

#lpx_simplex(leusch.lp)

solveSimplexGLPK(leusch.lp)


###################################################
### chunk number 23: 
###################################################

# Call the GLPK function glp_get_obj_val. 
# Usage getObjValGLPK(lp)

#leusch.obj <- lpx_get_obj_val(leusch.lp)
leusch.obj <- getObjValGLPK(leusch.lp)
leusch.obj 



###################################################
### chunk number 24: 
###################################################
leusch.col.rpt <- get.col.report(leusch.lp)


###################################################
### chunk number 25: 
###################################################
ac.per <- 
   matrix(as.numeric(as.matrix(I(leusch.col.rpt$activity))), 
          6, 8)
ac.per


###################################################
### chunk number 26: 
###################################################
rowSums(ac.per)


###################################################
### chunk number 27: 
###################################################
vol.cut.a <- c(199.4,197.4,203.4,211.4,217.6,195.4)



###################################################
### chunk number 28: leuschner-area-control
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
### chunk number 29: fig-leuschner-area-control
###################################################
ac.per <- NULL
coef.per <- NULL

ac.per <- 
  matrix(as.numeric(as.matrix(I(leusch.col.rpt$activity))), 
         6, 8)
ac.per

coef.per <- matrix(as.numeric(as.matrix(I(leusch.col.rpt$coef))), 6, 8)
coef.per

vol.cut.a <- c(rowSums(ac.per *  coef.per) / 1000)

#<- c(199.400,197.400,203.400,211.400,217.600,195.400)
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
### chunk number 30: 
###################################################
leusch.col.rpt[1,] 


###################################################
### chunk number 31: 
###################################################
leusch.row.rpt <- get.row.report(leusch.lp)
leusch.row.rpt



###################################################
### chunk number 32: 
###################################################

leusch.row.rpt[1,]



###################################################
### chunk number 33: 
###################################################

# Call the GLPK function glp_write_mps. 
# Usage writeMPSGLPK(lp, fmt, fname)

#lpx_write_mps(leusch.lp, "leuschner.mps")

writeMPSGLPK(leusch.lp, GLP_MPS_FILE, "leuschner.mps")

###################################################
### chunk number 34: write leush
###################################################

# Call the GLPK function glp_write_lp. 
# Usage writeLPGLPK(lp, fname)

#lpx_write_cpxlp(leusch.lp, "leuschner.xlp")

writeLPGLPK(leusch.lp, "leuschner.xlp")

###################################################
### chunk number 35: 
###################################################

# Call the GLPK function glp_delete_prob. Consult the GLPK documentation
# Usage delProbGLPK(lp)

#lpx_delete_prob(leusch.lp)

delProbGLPK(leusch.lp)
leusch.lp


###################################################
### chunk number 36: 
###################################################
system("rm -fr package-Ch9")
package.skeleton(name = "package-Ch9")



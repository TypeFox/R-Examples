### R code from vignette source 'cias.Rnw'

###################################################
### code chunk number 1: set_seed_chunk
###################################################
set.seed(0)


###################################################
### code chunk number 2: time_saver
###################################################
calc_from_scratch <- TRUE


###################################################
### code chunk number 3: cias.Rnw:82-85
###################################################
ignore <- require(multivator,quietly=TRUE)
ignore <- require(abind,quietly=TRUE)
ignore <- require(emulator,quietly=TRUE)


###################################################
### code chunk number 4: setupcias
###################################################
jj <- latin.hypercube(21,6)
colnames(jj) <- c("cias1","cias2","A_p1","B_p1","B_p2","C_p1")
rownames(jj) <- c(
                 paste("module_A_run",1:7,sep=""),
                 paste("module_B_run",1:7,sep=""),
                 paste("module_C_run",1:7,sep="")
                 )
jj[1:7  , 4:6  ] <- 0
jj[8:14 ,c(3,6)] <- 0
jj[15:21, 3:5  ] <- 0

real_cias_mdm <- mdm(jj, factor(rep(LETTERS[1:3],each=7)))


###################################################
### code chunk number 5: makedisplayableciasmdm
###################################################
                                        # have to create a version that looks good.  The
                                        # difference is that '0.000' displays as '0'.


jj <- xold(real_cias_mdm)
jj <- round(jj,3)
storage.mode(jj) <- 'character'
jj[nchar(jj)==3] <- paste(jj[nchar(jj)==3] ,'0',sep='')
jj[nchar(jj)==4] <- paste(jj[nchar(jj)==4] ,'0',sep='')
cias_mdm <- noquote(jj)


###################################################
### code chunk number 6: showcias
###################################################
cias_mdm


###################################################
### code chunk number 7: definemhp
###################################################
jjM <- matrix(1,3,3)
diag(jjM) <- 2

jjB <- matrix(0,6,6)
diag(jjB) <- 1
jjB <- abind(jjB,jjB,jjB,along=3)

cias_mhp <- mhp(M=jjM, B = jjB, levels=levels(real_cias_mdm),names=names(real_cias_mdm))


###################################################
### code chunk number 8: cheat
###################################################
cias_mdm <- real_cias_mdm


###################################################
### code chunk number 9: showsummary
###################################################
summary(cias_mhp)


###################################################
### code chunk number 10: cias.Rnw:187-194
###################################################
cias_LoF <- list(
                 A = function(x){ c(const=1,x[1:2],x[3  ]) },
                 B = function(x){ c(const=1,x[1:2],x[4:5]) },
                 C = function(x){ c(const=1,x[1:2],x[6  ]) }
                 )
cias_beta <- 1:13



###################################################
### code chunk number 11: showlof
###################################################
cias_LoF


###################################################
### code chunk number 12: dosomestuff
###################################################
cias_obs <- obs_maker(cias_mdm, cias_mhp, cias_LoF, cias_beta)
cias_expt <- experiment(cias_mdm , cias_obs)


###################################################
### code chunk number 13: definineunk
###################################################
jj <- cias_mdm[1:3,]
types(jj) <- levels(jj)
xold(jj)[,3] <- 0
xold(jj)[,1:2] <- 0.5
rownames(jj) <- paste("m",LETTERS[1:3],sep=".")
cias_unknown <- jj


###################################################
### code chunk number 14: showunk
###################################################
cias_unknown


###################################################
### code chunk number 15: usemultem
###################################################
multem(cias_unknown, cias_expt, cias_mhp, cias_LoF, give=TRUE)



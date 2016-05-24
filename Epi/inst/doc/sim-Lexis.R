### R code from vignette source 'sim-Lexis.rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: start
###################################################
options( width=90 )
library( Epi )
print( sessionInfo(), l=F )


###################################################
### code chunk number 2: Lexis
###################################################
data(DMlate)
dml <- Lexis( entry = list(Per=dodm, Age=dodm-dobth, DMdur=0 ),
               exit = list(Per=dox),
        exit.status = factor(!is.na(dodth),labels=c("DM","Dead")),
               data = DMlate )


###################################################
### code chunk number 3: cut
###################################################
dmi <- cutLexis( dml, cut = dml$doins,
                      pre = "DM",
                new.state = "Ins",
                new.scale = "t.Ins",
             split.states = TRUE )
summary( dmi )
str(dmi)


###################################################
### code chunk number 4: boxes
###################################################
boxes( dmi, boxpos = list(x=c(20,20,80,80),
                        y=c(80,20,80,20)),
            scale.R = 1000, show.BE = TRUE )


###################################################
### code chunk number 5: split
###################################################
Si <- splitLexis( dmi, 0:30/2, "DMdur" )
dim( Si )
print( subset( Si, lex.id==97 )[,1:10], digits=6 )


###################################################
### code chunk number 6: knots
###################################################
nk <- 5
( ai.kn <- with( subset(Si,lex.Xst=="Ins"),
                 quantile( Age+lex.dur  , probs=(1:nk-0.5)/nk ) ) )
( ad.kn <- with( subset(Si,lex.Xst=="Dead"),
                 quantile( Age+lex.dur  , probs=(1:nk-0.5)/nk ) ) )
( di.kn <- with( subset(Si,lex.Xst=="Ins"),
                 c(0,quantile( DMdur+lex.dur, probs=(1:(nk-1))/nk ) )) )
( dd.kn <- with( subset(Si,lex.Xst=="Dead"),
                 c(0,quantile( DMdur+lex.dur, probs=(1:(nk-1))/nk ) )) )
( ti.kn <- with( subset(Si,lex.Xst=="Dead(Ins)"),
                 c(0,quantile( t.Ins+lex.dur, probs=(1:(nk-1))/nk ) )) )


###################################################
### code chunk number 7: Poisson
###################################################
library( splines )
DM.Ins <- glm( (lex.Xst=="Ins") ~ Ns( Age  , knots=ai.kn ) +
                                  Ns( DMdur, knots=di.kn ) +
                                  I(Per-2000) + sex,
               family=poisson, offset=log(lex.dur),
               data = subset(Si,lex.Cst=="DM") )
DM.Dead <- glm( (lex.Xst=="Dead") ~ Ns( Age  , knots=ad.kn ) +
                                    Ns( DMdur, knots=dd.kn ) +
                                    I(Per-2000) + sex,
               family=poisson, offset=log(lex.dur),
               data = subset(Si,lex.Cst=="DM") )
Ins.Dead <- glm( (lex.Xst=="Dead(Ins)") ~ Ns( Age  , knots=ad.kn ) +
                                          Ns( DMdur, knots=dd.kn ) +
                                          Ns( t.Ins, knots=ti.kn ) +
                                          I(Per-2000) + sex,
               family=poisson, offset=log(lex.dur),
               data = subset(Si,lex.Cst=="Ins") )


###################################################
### code chunk number 8: prop-haz
###################################################
with( Si, table(lex.Cst) )
All.Dead <- glm( (lex.Xst %in% c("Dead(Ins)","Dead")) ~
                                          Ns( Age  , knots=ad.kn ) +
                                          Ns( DMdur, knots=dd.kn ) +
                                          lex.Cst +
                                          I(Per-2000) + sex,
               family=poisson, offset=log(lex.dur),
               data = Si )
round( ci.exp( All.Dead ), 3 )


###################################################
### code chunk number 9: get-dev
###################################################
what <- c("null.deviance","df.null","deviance","df.residual")
( rD <- unlist(  DM.Dead[what] ) )
( rI <- unlist( Ins.Dead[what] ) )
( rA <- unlist( All.Dead[what] ) )
round( c( dd <- rA-(rI+rD), "pVal"=1-pchisq(dd[3],dd[4]+1) ), 3 )


###################################################
### code chunk number 10: pr-array
###################################################
pr.rates <- NArray( list( DMdur = seq(0,12,0.1),
                          DMage = 4:7*10,
                          r.Ins = c(NA,0,2,5),
                          model = c("DM/Ins","All"),
                           what = c("rate","lo","hi") ) )
str( pr.rates )


###################################################
### code chunk number 11: sim-Lexis.rnw:481-482
###################################################
ci.pred


###################################################
### code chunk number 12: make-pred
###################################################
nd <- data.frame( DMdur = as.numeric( dimnames(pr.rates)[[1]] ),
                lex.Cst = factor( 1, levels=1:4,
                                  labels=levels(Si$lex.Cst) ),
                    sex = factor( 1, levels=1:2, labels=c("M","F")),
                lex.dur = 1000 )
for( ia in dimnames(pr.rates)[[2]] )
   {
dnew <- transform( nd, Age = as.numeric(ia)+DMdur,
                       Per = 1998+DMdur )
pr.rates[,ia,1,"DM/Ins",] <- ci.pred(  DM.Dead, newdata = dnew )
pr.rates[,ia,1,"All"   ,] <- ci.pred( All.Dead, newdata = dnew )
for( ii in dimnames(pr.rates)[[3]][-1] )
   {
dnew = transform( dnew, lex.Cst = factor( 2, levels=1:4,
                                          labels=levels(Si$lex.Cst) ),
                          t.Ins = ifelse( (DMdur-as.numeric(ii)) >= 0,
                                           DMdur-as.numeric(ii), NA ) )
pr.rates[,ia, ii ,"DM/Ins",] <- ci.pred( Ins.Dead, newdata = dnew )
pr.rates[,ia, ii ,"All"   ,] <- ci.pred( All.Dead, newdata = dnew )
    }
    }


###################################################
### code chunk number 13: mort-int
###################################################
par( mar=c(3,3,1,1), mgp=c(3,1,0)/1.6, las=1 )
plot( NA, xlim=c(40,82), ylim=c(5,300), bty="n",
      log="y", xlab="Age", ylab="Mortality rate per 1000 PY" )
abline( v=seq(40,80,5), h=outer(1:9,10^(0:2),"*"), col=gray(0.8) )
for( aa in 4:7*10 ) for( ii in 1:4 )
   matlines( aa+as.numeric(dimnames(pr.rates)[[1]]),
            cbind( pr.rates[,paste(aa),ii,"DM/Ins",],
                   pr.rates[,paste(aa),ii,"All"   ,] ),
            type="l", lty=1, lwd=c(3,1,1),
            col=rep(c("red","limegreen"),each=3) )


###################################################
### code chunk number 14: Tr
###################################################
Tr <- list( "DM" = list( "Ins"       = DM.Ins,
                         "Dead"      = DM.Dead  ),
           "Ins" = list( "Dead(Ins)" = Ins.Dead ) )


###################################################
### code chunk number 15: make-ini
###################################################
str( Si[NULL,1:9] )
ini <- subset(Si,FALSE,select=1:9)
str( ini )
ini <- subset(Si,select=1:9)[NULL,]
str( ini )


###################################################
### code chunk number 16: ini-fill
###################################################
ini[1:2,"lex.id"] <- 1:2
ini[1:2,"lex.Cst"] <- "DM"
ini[1:2,"Per"] <- 1995
ini[1:2,"Age"] <- 60
ini[1:2,"DMdur"] <- 5
ini[1:2,"sex"] <- c("M","F")
ini


###################################################
### code chunk number 17: simL
###################################################
Nsim <- 5000
system.time( simL <- simLexis( Tr,
                              ini,
                          t.range = 12,
                                N = Nsim ) )


###################################################
### code chunk number 18: sum-simL
###################################################
summary( simL, by="sex" )


###################################################
### code chunk number 19: Tr.p-simP
###################################################
Tr.p <- list( "DM" = list( "Ins"       = DM.Ins,
                           "Dead"      = All.Dead  ),
             "Ins" = list( "Dead(Ins)" = All.Dead ) )
system.time( simP <- simLexis( Tr.p,
                                ini,
                            t.range = 12,
                                  N = Nsim ) )
summary( simP, by="sex" )


###################################################
### code chunk number 20: Cox-dur
###################################################
library( survival )
Cox.Dead <- coxph( Surv( DMdur, DMdur+lex.dur,
                         lex.Xst %in% c("Dead(Ins)","Dead")) ~
                   Ns( Age-DMdur, knots=ad.kn ) +
                   I(lex.Cst=="Ins") +
                   I(Per-2000) + sex,
               data = Si )
round( ci.exp( Cox.Dead ), 3 )
round( ci.exp( All.Dead ), 3 )


###################################################
### code chunk number 21: TR.c
###################################################
Tr.c <- list( "DM" = list( "Ins"       = Tr$DM$Ins,
                           "Dead"      = Cox.Dead  ),
             "Ins" = list( "Dead(Ins)" = Cox.Dead ) )
system.time( simC <- simLexis( Tr.c,
                                ini,
                            t.range = 12,
                                  N = Nsim ) )
summary( simC, by="sex" )


###################################################
### code chunk number 22: nState
###################################################
system.time(
nSt <- nState( subset(simL,sex=="M"),
               at=seq(0,11,0.2), from=1995, time.scale="Per" ) )
nSt[1:10,]


###################################################
### code chunk number 23: pstate0
###################################################
pM <- pState( nSt, perm=c(1,2,4,3) )
head( pM )
par( mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(3,1,0)/1.6 )
plot( pM )
plot( pM, border="black", col="transparent", lwd=3 )
text( rep(as.numeric(rownames(pM)[nrow(pM)-1]),ncol(pM)),
      pM[nrow(pM),]-diff(c(0,pM[nrow(pM),]))/5,
      colnames( pM ), adj=1 )
box( col="white", lwd=3 )
box()


###################################################
### code chunk number 24: pstatex
###################################################
clr <- c("limegreen","orange")
# expand with a lighter version of the two chosen colors
clx <- c( clr, rgb( t( col2rgb( clr[2:1] )*2 + rep(255,3) ) / 3, max=255 ) )
par( mfrow=c(1,2), las=1, mar=c(3,3,4,2), mgp=c(3,1,0)/1.6 )
# Men
plot( pM, col=clx )
lines( as.numeric(rownames(pM)), pM[,2], lwd=3 )
mtext( "60 year old male, diagnosed 1990, aged 55", side=3, line=2.5, adj=0, col=gray(0.6) )
mtext( "Survival curve", side=3, line=1.5, adj=0 )
mtext( "DM, no insulin   DM, Insulin", side=3, line=0.5, adj=0, col=clr[1] )
mtext( "DM, no insulin", side=3, line=0.5, adj=0, col=clr[2] )
axis( side=4 )
axis( side=4, at=1:19/20, labels=FALSE )
axis( side=4, at=1:99/100, labels=FALSE, tcl=-0.3 )
# Women
pF <- pState( nState( subset(simL,sex=="F"),
                      at=seq(0,11,0.2),
                      from=1995,
                      time.scale="Per" ),
              perm=c(1,2,4,3) )
plot( pF, col=clx )
lines( as.numeric(rownames(pF)), pF[,2], lwd=3 )
mtext( "60 year old female, diagnosed 1990, aged 55", side=3, line=2.5, adj=0, col=gray(0.6) )
mtext( "Survival curve", side=3, line=1.5, adj=0 )
mtext( "DM, no insulin   DM, Insulin", side=3, line=0.5, adj=0, col=clr[1] )
mtext( "DM, no insulin", side=3, line=0.5, adj=0, col=clr[2] )
axis( side=4 )
axis( side=4, at=1:19/20, labels=FALSE )
axis( side=4, at=1:99/100, labels=FALSE, tcl=-0.3 )


###################################################
### code chunk number 25: pstatey
###################################################
par( mfrow=c(1,2), las=1, mar=c(3,3,4,2), mgp=c(3,1,0)/1.6 )
# Men
pM <- pState( nState( subset(simL,sex=="M"),
                      at=seq(0,11,0.2),
                      from=60,
                      time.scale="Age" ),
              perm=c(1,2,4,3) )
plot( pM, col=clx, xlab="Age" )
lines( as.numeric(rownames(pM)), pM[,2], lwd=3 )
mtext( "60 year old male, diagnosed 1990, aged 55", side=3, line=2.5, adj=0, col=gray(0.6) )
mtext( "Survival curve", side=3, line=1.5, adj=0 )
mtext( "DM, no insulin   DM, Insulin", side=3, line=0.5, adj=0, col=clr[1] )
mtext( "DM, no insulin", side=3, line=0.5, adj=0, col=clr[2] )
axis( side=4 )
axis( side=4, at=1:19/20, labels=FALSE )
axis( side=4, at=1:19/20, labels=FALSE, tcl=-0.4 )
axis( side=4, at=1:99/100, labels=FALSE, tcl=-0.3 )
# Women
pF <- pState( nState( subset(simL,sex=="F"),
                      at=seq(0,11,0.2),
                      from=60,
                      time.scale="Age" ),
              perm=c(1,2,4,3) )
plot( pF, col=clx, xlab="Age" )
lines( as.numeric(rownames(pF)), pF[,2], lwd=3 )
mtext( "60 year old female, diagnosed 1990, aged 55", side=3, line=2.5, adj=0, col=gray(0.6) )
mtext( "Survival curve", side=3, line=1.5, adj=0 )
mtext( "DM, no insulin   DM, Insulin", side=3, line=0.5, adj=0, col=clr[1] )
mtext( "DM, no insulin", side=3, line=0.5, adj=0, col=clr[2] )
axis( side=4 )
axis( side=4, at=1:9/10, labels=FALSE )
axis( side=4, at=1:19/20, labels=FALSE, tcl=-0.4 )
axis( side=4, at=1:99/100, labels=FALSE, tcl=-0.3 )


###################################################
### code chunk number 26: comp-0
###################################################
PrM  <- pState( nState( subset(simP,sex=="M"),
                        at=seq(0,11,0.2),
                        from=60,
                        time.scale="Age" ),
                perm=c(1,2,4,3) )
PrF  <- pState( nState( subset(simP,sex=="F"),
                        at=seq(0,11,0.2),
                        from=60,
                        time.scale="Age" ),
                perm=c(1,2,4,3) )
CoxM <- pState( nState( subset(simC,sex=="M"),
                        at=seq(0,11,0.2),
                        from=60,
                        time.scale="Age" ),
                perm=c(1,2,4,3) )
CoxF <- pState( nState( subset(simC,sex=="F"),
                        at=seq(0,11,0.2),
                        from=60,
                        time.scale="Age" ),
                perm=c(1,2,4,3) )

par( mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(3,1,0)/1.6 )
 plot(   pM, border="black", col="transparent", lwd=3 )
lines(  PrM, border="blue" , col="transparent", lwd=3 )
lines( CoxM, border="red"  , col="transparent", lwd=3 )
text( 60.5, 0.05, "M" )
box( lwd=3 )

 plot(   pF, border="black", col="transparent", lwd=3 )
lines(  PrF, border="blue" , col="transparent", lwd=3 )
lines( CoxF, border="red"  , col="transparent", lwd=3 )
text( 60.5, 0.05, "F" )
box( lwd=3 )


###################################################
### code chunk number 27: sim-Lexis.rnw:987-988
###################################################
options( keep.source=TRUE )


###################################################
### code chunk number 28: sim-Lexis.rnw:1004-1007
###################################################
cbind(
attr( ini, "time.scale" ),
attr( ini, "time.since" ) )


###################################################
### code chunk number 29: sim-Lexis.rnw:1032-1033
###################################################
simLexis


###################################################
### code chunk number 30: sim-Lexis.rnw:1050-1051
###################################################
Epi:::simX


###################################################
### code chunk number 31: sim-Lexis.rnw:1063-1064
###################################################
Epi:::sim1


###################################################
### code chunk number 32: sim-Lexis.rnw:1076-1077
###################################################
Epi:::lint


###################################################
### code chunk number 33: sim-Lexis.rnw:1087-1088
###################################################
Epi:::get.next


###################################################
### code chunk number 34: sim-Lexis.rnw:1097-1098
###################################################
Epi:::chop.lex


###################################################
### code chunk number 35: sim-Lexis.rnw:1115-1116
###################################################
nState


###################################################
### code chunk number 36: sim-Lexis.rnw:1125-1126
###################################################
pState


###################################################
### code chunk number 37: sim-Lexis.rnw:1130-1132
###################################################
plot.pState
lines.pState



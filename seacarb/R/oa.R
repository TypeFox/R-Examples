# Copyright (C) 2009 Jean-Pierre Gattuso and H?loise Lavigne
# with a most valuable contribution of Bernard Gentili <gentili@obs-vlfr.fr>
# and valuable suggestions from Jean-Marie Epitalon <epitalon@lsce.saclay.cea.fr>
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
#
oa <-
function(flag, var1, var2, pCO2f, pCO2s=1e6, S=35, T=25, P=0, Pt=0, Sit=0, k1k2='x', kf='x', ks="d", pHscale="T", plot=FALSE,  b="u74"){

n <- max(length(var1), length(var2), length(pCO2f), length(pCO2s), length(S), length(T), length(P), length(Pt), length(Sit), length(k1k2), length(kf), length(pHscale), length(ks))
if(length(flag)!=n){ flag <- rep(flag[1],n)}
if(length(var1)!=n){ var1 <- rep(var1[1],n)}
if(length(var2)!=n){ var2 <- rep(var2[1],n)}
if(length(pCO2f) !=n){ pCO2f <- rep(pCO2f[1],n)}
if(length(pCO2s) !=n){ pCO2s <- rep(pCO2s[1],n)}
if(length(S)!=n){ S <- rep(S[1],n)}
if(length(T)!=n){ T <- rep(T[1],n)}
if(length(P)!=n){ P <- rep(P[1],n)}
if(length(Pt)!=n){ Pt <- rep(Pt[1],n)}
if(length(Sit)!=n){ Sit <- rep(Sit[1],n)}
if(length(k1k2)!=n){ k1k2 <- rep(k1k2[1],n)}
if(length(kf)!=n){ kf <- rep(kf[1],n)}
if(length(ks)!=n){ ks <- rep(ks[1],n)}
if(length(pHscale)!=n){pHscale <- rep(pHscale[1],n)}
if(length(b)!=n){ b <- rep(b[1],n)}

# if the concentrations of total silicate and total phosphate are NA
# they are set at 0
Sit[is.na(Sit)] <- 0
Pt[is.na(Pt)] <- 0

# initial system
# --------------
SYSi <- carb(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Sit=Sit, Pt=Pt, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)

# CO2 bubbling
# -----------------------------
SYSfCO2bubbling <- carb(flag=24, var1=pCO2f, var2=SYSi$ALK, S=S, T=T, P=P, Sit=Sit, Pt=Pt, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)


# Seawater mixing
# ---------------------------------
SYSs <- carb(flag=24, var1=pCO2s, var2=SYSi$ALK, S=S, T=T, P=P, Sit=Sit, Pt=Pt, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)

SYSf <- carb(flag=24, var1=pCO2f, var2=SYSi$ALK, S=S, T=T, P=P, Sit=Sit, Pt=Pt, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)

DICi <- SYSi$DIC
DICs <- SYSs$DIC
DICf <- SYSf$DIC

ws <- (DICf-DICi)/(DICs-DICi)

ws <- round(ws, 8)
wi <- round(1-ws, 8)

SYSfSWmixing <- SYSf

# Addition of strong acid
# -----------------------------------

SYSfacid <- carb(flag=25, var1=pCO2f, var2=SYSi$DIC, S=S, T=T, P=P, Sit=Sit, Pt=Pt, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
ALKf <- SYSfacid$ALK
ALKi <- SYSi$ALK

hplus <- ALKi-ALKf

hplus <- round(hplus, 8)

## on va devoir ajouter hplus mole par kg de solution de HCl

# Addition of strong acid and carbonate
# ------------------------------------------------------
SYSfCarboAcid <- carb(flag=24, var1=pCO2f, var2=SYSi$ALK, S=S, T=T, P=P, Sit=Sit, Pt=Pt, k1k2=k1k2, kf=kf,  ks=ks, pHscale=pHscale, b=b)

deltaDIC <- SYSfCarboAcid$DIC - SYSi$DIC

## with CO3 addition
TAint <- SYSi$ALK + 2*deltaDIC
TAintCO3 <- TAint  ## we keep the value for the plot

DICint <- SYSfCarboAcid$DIC

SYSint <- carb(flag=15, var1=TAint, var2=DICint, S=S, T=T, P=P, Sit=Sit, Pt=Pt, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)

deltaTA <- SYSf$ALK-SYSint$ALK

CO3 <- round(deltaDIC, 8)
ACIDco3 <- round(-deltaTA, 8)

## with HCO3 addition
TAint <- SYSi$ALK + deltaDIC
TAintHCO3 <- TAint  ## we keep the value for the plot

DICint <- SYSfCarboAcid$DIC

SYSint <- carb(flag=15, var1=TAint, var2=DICint, S=S, T=T, P=P, Sit=Sit, Pt=Pt, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)

deltaTA <- SYSf$ALK-SYSint$ALK

HCO3 <- round(deltaDIC, 8)
ACIDhco3 <- round(-deltaTA, 8)


###############################################################################
## OUTPUTS
###############################################################################

OUT <- new.env()

## Summary table for sw state after each experiments
## -------------------------------------------------
co <- as.data.frame(c(rep("initial",n), rep("final-CO2 bubbling", n), rep("final-SW mixing", n), rep("final-acid", n), rep("final-HCO3-and-acid", n), rep("final-CO3-and-acid", n)))

State <- rbind(SYSi, SYSfCO2bubbling, SYSfSWmixing, SYSfacid ,SYSfCarboAcid, SYSfCarboAcid)

State <- cbind(co, State)

names(State)[1] <- "comment"

OUT$summary <- State


## Summary table for manipulations of each experiments
## ---------------------------------------------------

col1 <- c(rep("CO2 bubbling", n), rep("Mixing", n), rep("Addition of acid", n), rep("Addition of HCO3 and acid", n), rep("Addition of CO3 and acid", n))
col2 <- c(rep("Air pCO2 (uatm)", n), rep("Weight fraction high-CO2", n), rep("H+ (mol/kg)", n), rep("HCO3 (mol/kg)", n), rep("CO3 (mol/kg)", n))
col3 <- c(pCO2f,ws,hplus,HCO3,CO3)
col4 <- c(rep(NA,n*3),rep("H+ (mol/kg)", n),rep("H+ (mol/kg)", n)) 
col5 <- c(rep(NA,n*3), ACIDhco3, ACIDco3)

perturbation <- cbind(col1, col2, col3, col4, col5)
perturbation <- data.frame(perturbation)
names(perturbation) <- c("Method", "1st param.", "Value 1st param.", "2nd param.", "Value 2nd param.")
OUT$perturbation <- perturbation

## Description of manipulation
## ---------------------------------------------------

col1 <- c("CO2 bubbling", "Mixing", "Addition of acid", "Addition of HCO3 and acid", "Addition of CO3 and acid")
col2 <- c( paste("Bubble seawater with a gas with pCO2=", round(pCO2f, 3)[1], "uatm", sep=""),
paste("Mix", wi[1], "kg of the normal seawater with", ws[1] ,"kg of high-pCO2 seawater" ),
paste("Add", hplus[1], "mol of H+ to 1 kg of the initial seawater"),
paste("Add", HCO3[1], "mol of HCO3 and then", ACIDhco3[1], "mol of H+ to 1 kg of initial seawater" ),
paste("Add", CO3[1], "mol of CO3 and then", ACIDco3[1], "mol of H+ to 1 kg of initial seawater" ))
description <- cbind(col1, col2)
description <- data.frame(description)
names(description) <- c("Approach", "Description of manipulations")
OUT$description <- description

OUT <- as.list(OUT)


########### GRAPHIQUE

if(plot){
par(mai=c(1,1,1.6,1))
x1 <- SYSi$DIC[1] - 0.1*(SYSf$DIC[1]-SYSi$DIC[1])  # xlim inf
x2 <- SYSf$DIC[1] + 0.1*(SYSf$DIC[1]-SYSi$DIC[1])  # xlim sup
y1 <- SYSfacid$ALK[1] - 0.1*(TAintCO3[1]-SYSfacid$ALK[1])  # ylim inf
y2 <- TAintCO3[1] + 0.1*(TAintCO3[1]-SYSfacid$ALK[1])       # ylim sup
xa <- SYSi$DIC[1]   # starting point
ya <- SYSi$ALK[1]   # starting point
plot(xa ,ya , xlim=c(x1,x2), ylim=c(y1, y2),
	xlab=expression(paste("Dissolved inorganic carbon (mol ", kg^"-1", ")", sep="" )), 
	ylab=expression(paste("Total alkalinity (mol ", kg^"-1", ")", sep="" )),
	pch=19)

## Bubbling and sea watermixing
xb <- SYSf$DIC[1]
yb <- SYSf$ALK[1] 
arrows(xa, ya, xb, yb, col=4, lwd=2)

## strong acid
xb <- SYSfacid$DIC[1]
yb <- SYSfacid$ALK[1]
arrows(xa, ya, xb, yb, col=3, lwd=2)

## strong acid  + base
# CO3
xb <- SYSfCarboAcid$DIC[1]
yb <- TAintCO3[1]
arrows(xa, ya, xb, yb, col=2, lwd=2)
ybb <- SYSfCarboAcid$ALK[1]
arrows(xb, yb, xb, ybb, col=2, lwd=2)
# HCO3
xb <- SYSfCarboAcid$DIC[1]
yb <- TAintHCO3[1]
arrows(xa, ya, xb, yb, col="orange", lwd=2)
ybb <- SYSfCarboAcid$ALK[1]
arrows(xb, yb, xb, ybb, col="orange", lty=2, lwd=2)

alk <- seq(y1, y2,, 15)
dic <- seq(x1, x2,, 15)
z <- expand.grid(ALK=alk, DIC=dic)
cp <- carb(flag=15, var1=z$ALK, var2=z$DIC, T=T[1], P=P[1], S=S[1], Pt=Pt[1], Sit=Sit[1], pHscale=pHscale[1], k1k2=k1k2[1], kf=kf[1])
pCO2 <- cp$pCO2
M <- tapply(pCO2, z, unique)
contour(dic, alk, t(M),
	levels=c(seq(0, 10000, 100)),
	labcex=0.5,
	method="flattest",
	col=1,
	lwd=1,
	lty="dashed",
	add=TRUE,
	axes = FALSE, frame = FALSE
	)
	mtext(side=3, text=expression(paste(CO[2], " bubbling and seawater mixing")), col=4, line=5, adj=0)
	mtext(side=3, text="Addition of strong acid", col=3, line=3.5, adj=0)
	mtext(side=3, text=expression(paste("Addition of ", HCO[3]^"-", " and strong acid")), col="orange" , line=2, adj=0)
	mtext(side=3, text=expression(paste("Addition of ", CO[3]^"2-", " and strong acid")), col=2 , line=0.5, adj=0)
}
 return(OUT)
}

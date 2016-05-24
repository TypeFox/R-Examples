## ----include=FALSE, cache=FALSE------------------------------------------
library(knitr)
#opts_chunk$set(cache=FALSE,tidy=FALSE,highlight=FALSE)
opts_chunk$set(cache = FALSE, tidy = FALSE, out.height = "4.0in", out.width = "4.0in", fig.height = 5.2, fig.width = 5.2, fig.align = "center")
library(biogas)

## ----echo=FALSE----------------------------------------------------------
  options(width=75)

  #listing <- function(x, options) {
  #  paste("\\begin{lstlisting}[basicstyle=\\ttfamily,breaklines=true]\n",
  #    x, "\\end{lstlisting}\n", sep = "")
  #}
  #knit_hooks$set(source=listing, output=listing)

  # biogas functions
  files <- list.files('~/Dropbox/biogas_package/biogas/R', full.names = TRUE)
  for(i in files) source(i)

## ----install, eval = FALSE-----------------------------------------------
#    install.packages("biogas")

## ----load, eval = FALSE--------------------------------------------------
#    library(biogas)

## ----cellulose1----------------------------------------------------------
  predBg(form = "C6H10O5")

## ----cellulose1.5--------------------------------------------------------
  predBg("C6H10O5")

## ----args1---------------------------------------------------------------
  args(predBg)

## ----cellulose2----------------------------------------------------------
  predBg("C6H10O5", mass = 5.75)

## ----cellulose3----------------------------------------------------------
  predBg("C6H10O5", mass = 5.75, fs = 0)
  predBg("C6H10O5", mass = 5.75, fs = 0.05)
  predBg("C6H10O5", mass = 5.75, fs = 0.20)

## ----cellulose3.5--------------------------------------------------------
  predBg("C6H10O5", mass = 5.75, fs = c(0, 0.05, 0.2))

## ----afp1----------------------------------------------------------------
  predBg(c("CHOOH", "CH3COOH", "CH3CH2COOH"), mol = 1)
  predBg(c("CHOOH", "CH3COOH", "CH3CH2COOH"), mass = 1)

## ----fa1-----------------------------------------------------------------
  predBg("(CHOOH)0.5 (CH3COOH)0.5", mass = 1)
  predBg("(CHOOH)0.5 (CH3COOH)0.5", mol = 1)

## ----fa2-----------------------------------------------------------------
  50/predBg("(CHOOH)0.5 (CH3COOH)0.5", mol = 1)

## ----fa3-----------------------------------------------------------------
  predBg("(CHOOH)0.5 (CH3COOH)0.5", mol = 0.0036)

## ----afoutput1-----------------------------------------------------------
  predBg("(CHOOH)0.5 (CH3COOH)0.5", mol = 0.0036, value = "all")

## ----afoutput2-----------------------------------------------------------
  predBg("(CHOOH)0.5 (CH3COOH)0.5", mol = 0.0036, fs = c(0, 0.01, 0.05, 0.1), 
	 value = "all")

## ----afoutput3-----------------------------------------------------------
  predBg("(CHOOH)0.5 (CH3COOH)0.5", mol = 0.06, fs = c(0, 0.05, 0.08, 0.1, 0.15), 
	 value = "all")

## ----cod1----------------------------------------------------------------
  predBg(COD = 2.6)

## ----cod2----------------------------------------------------------------
  predBg(COD = 2.6, fd = 0.6, fs = 0.1)

## ----manure4-------------------------------------------------------------
  predBg(mcomp = c(carbohydrate = 0.682, protein = 0.158, lipid = 0.054, 
		   VFA = 0.031, lignin = 0.075), 
	 mass = 1, fd = 0.4, fs = 0.1)

## ----manure4.5-----------------------------------------------------------
  predBg(mcomp = c(carbohydrate = 0.682, protein = 0.158, lipid = 0.054, 
		   VFA = 0.031, lignin = 0.075), 
	 mass = 1, fd = 0.4, fs = 0.1, value = "all")

## ----manure4.6-----------------------------------------------------------
  predBg(mcomp = c(carbohydrate = 0.682, protein = 0.158, lipid = 0.054, 
		   VFA = 0.031, lignin = 0.075), 
	 mass = 1, fd = 0.4, fs = 0.1, shortform = FALSE, value = "all")

## ----manure4.7-----------------------------------------------------------
  predBg("C29.2H47O18.9N", mass = 1, fd = 0.4, fs = 0.1, shortform = FALSE, value = "all")

## ----manure5-------------------------------------------------------------
  predBg(mcomp = c(carbohydrate = 0.682, protein = 0.158, lipid = 0.054, 
		   VFA = 0.031, lignin = 0.075, 
		   C3H8O3 = 0.25),
	 mass = 1, fd = 0.4, fs = 0.1)

## ----manure5.5-----------------------------------------------------------
  predBg(mcomp = c(C29.2H47O18.9N = 0.8, C3H8O3 = 0.2), mass = 1, fd = 0.4, fs = 0.1)

## ----mcomp7--------------------------------------------------------------
  predBg(mcomp = c(C6H10O5 = 5, C54H100O7 = 1), mass = 1)

## ----xxx-----------------------------------------------------------------
  1/5*molMass("C6H10O5")/molMass("C54H100O7")

## ----mcomp8--------------------------------------------------------------
  predBg("(C6H10O5)1 (C54H100O7)0.037648")

## ----bgcomp1-------------------------------------------------------------
  predBg(mcomp = c(C6H10O5 = 5/6, C54H100O7 = 1/6), mass = 1, value = "all")

## ----subconc1------------------------------------------------------------
  predBg(mcomp = c(C6H10O5 = 5/6, C54H100O7 = 1/6), mass = 1, 
	 fd = 0.8, fs = 0.1, conc.sub = 50, pH = 7.5, temp = 35, 
	 value = "all")

## ----subconc2------------------------------------------------------------
  bg1 <- predBg(mcomp = c(C6H10O5 = 5/6, C54H100O7 = 1/6), mass = 1, 
		fd = 0.8, fs = 0.1, conc.sub = 50, pH = c(6.5 + 0:10*0.2), 
		temp = 35, value = "all")

## ----subconc2plot1-------------------------------------------------------
  plot(xCH4 ~ pH, data = bg1, type = 'o', col = "red")

## ----subconc2plot2-------------------------------------------------------
  plot(vBg ~ pH, data = bg1, type = 'o', ylim = c(0, max(bg1$vBg)), col = "blue")

## ----fsplot, echo = FALSE------------------------------------------------
  hrt <- 5:100
  fbd <- 0.8
  fsfa <- 0.06*((1 + (1 - fbd)*0.03*hrt)/(1 + 0.03*hrt))
  fscarb <- 0.28*((1 + (1 - fbd)*0.05*hrt)/(1 + 0.05*hrt))
  fssludge <- 0.11*((1 + (1 - fbd)*0.05*hrt)/(1 + 0.05*hrt))
  plot(hrt, fsfa, type = 'l', xlab = 'Solids retention time (d)', ylab = expression(f[s]~(fraction)), ylim = c(0, 0.25), col = 'blue', lwd = 2)
  lines(hrt, fscarb, col = 'red', lwd = 2)
  lines(hrt, fssludge, col = 'green', lwd = 2)
  grid(col = 'gray45')
  text(c(3, 3, 3), c(0.03, 0.1, 0.24), c('Fatty acids', 'Sludge', 'Carbohydrates'), pos = 4)


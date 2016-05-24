#####
# File:     Script - EPA09, Chapter 09 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 9 of the EPA Guidance Document
#
#           USEPA. (2009).  "Statistical Analysis of Ground Water 
#             Monitoring Data at RCRA Facilities, Unified Guidance."
#             EPA 530-R-09-007.
#             Office of Resource Conservation and Recovery, 
#             Program Information and Implementation Division.
#             March 2009.
#
# Author:   Steven P. Millard
#
# Last
# Updated:  August 20, 2012
#####

#NOTE:  Unless otherwise noted, differences between what is shown in the
#       EPA Guidance Document and results shown here are due to the 
#       EPA Guidance Document rounding intermediate results prior to 
#       the final result.

#######################################################################################

library(EnvStats)


# Example 9-1, pp. 9-3 to 9-5
#------------------------------

stats <- summaryStats(TCE.mg.per.L ~ Well, 
	data = EPA.09.Table.9.1.TCE.df,
	combine.groups = F, 
	digits = 3, stats.in.rows = TRUE)
stats
#        Well.1 Well.2
#N       14     13    
#Mean     0.063  0.118
#SD       0.079  0.02 
#Median   0.031  0.11 
#Min      0.004  0.099
#Max      0.25   0.17 
#NA's     1      2    
#N.Total 15     15   

Min.Max.mat <- stats[c("Min", "Max"), ]
Min.Max.mat
#    Well.1 Well.2
#Min  0.004  0.099
#Max  0.250  0.170

TCE            <- EPA.09.Table.9.1.TCE.df$TCE.mg.per.L
Date           <- EPA.09.Table.9.1.TCE.df$Date
Well           <- EPA.09.Table.9.1.TCE.df$Well
Data.Qualifier <- EPA.09.Table.9.1.TCE.df$Data.Qualifier

index.Well.1.Min <- Well == "Well.1" & !is.na(TCE) & TCE == Min.Max.mat["Min", "Well.1"]
index.Well.1.Max <- Well == "Well.1" & !is.na(TCE) & TCE == Min.Max.mat["Max", "Well.1"]
index.Well.2.Min <- Well == "Well.2" & !is.na(TCE) & TCE == Min.Max.mat["Min", "Well.2"]
index.Well.2.Max <- Well == "Well.2" & !is.na(TCE) & TCE == Min.Max.mat["Max", "Well.2"]

DQ.vec <- c(paste(Data.Qualifier[index.Well.1.Min], collapse = ","), 
	paste(Data.Qualifier[index.Well.1.Max], collapse = ","),
	paste(Data.Qualifier[index.Well.2.Min], collapse = ","), 
	paste(Data.Qualifier[index.Well.2.Max], collapse = ","))
DQ.mat <- matrix(DQ.vec, 2, 2)
DQ.mat
#     [,1]  [,2]
#[1,] "J,U" "U" 
#[2,] ""    "" 

df <- data.frame(Min.Max.mat[, 1], DQ.mat[, 1], Min.Max.mat[, 2], DQ.mat[, 2])
row.names(df) <- c("Min", "Max")
names(df) <- c("Well.1.TCE", "Well.1.DQ", "Well.2.TCE", "Well.2.DQ")
df
#    Well.1.TCE Well.1.DQ Well.2.TCE Well.2.DQ
#Min      0.004       J,U      0.099         U
#Max      0.250                0.170 

Well.index <- Well == "Well.1"
ND.index   <- Data.Qualifier == "U"
NA.index   <- is.na(TCE)
windows()
plot(Date[Well.index], TCE[Well.index], type = "n", 
	ylim = c(0, max(TCE, na.rm = T)),
	xlab = "Time", ylab = "Trichloroethene (mg/L)")
points(Date[Well.index & !ND.index], TCE[Well.index & !ND.index], 
	cex = 2, pch = 16, col = "blue")
points(Date[Well.index & ND.index], TCE[Well.index & ND.index], 
	cex = 2, pch = 1, col = "blue")
lines(Date[Well.index & !NA.index], TCE[Well.index & !NA.index], 
	lwd = 2, col = "blue")

Well.index <- Well == "Well.2"
points(Date[Well.index & !ND.index], TCE[Well.index & !ND.index], 
	cex = 2, pch = 15, col = "red")
points(Date[Well.index & ND.index], TCE[Well.index & ND.index], 
	cex = 2, pch = 14, col = "red")
lines(Date[Well.index & !NA.index], TCE[Well.index & !NA.index], 
	lwd = 2, lty = 2, col = "red")

legend(Date[1], 0.25, legend = c("Well 1", "Well 2"),
	col = c("blue", "red"), lty = c(1, 2), lwd = rep(2, 2),
	pch = c(16, 15), horiz = T)

title(main = "Figure 9-1. Time Series Plot of Trichloroethene Groudwater\nfor Wells 1 and 2 from 2005-2007")

rm(TCE, Date, Well, Data.Qualifier, stats, Min.Max.mat, 
	index.Well.1.Min, index.Well.1.Max, 
	index.Well.2.Min, index.Well.2.Max, 
	DQ.vec, DQ.mat, df, Well.index, ND.index, NA.index) 


###################################################################
	

# Example 9-2, pp. 9-7 to 9-8, Using Standard Boxplot
#----------------------------------------------------
# NOTE: Instructions for producing Figure 9-2 include
#       "Step 6. Draw the whiskers from each end of the 
#       box to the furthest data point to show the 
#       full range of the data."
#
#       This is not how boxplots are conventionally drawn,
#       so instead the standard boxplot is drawn here, where
#       the lower whisker extends to the smallest value within
#       1.5 x IQR of the 25th Percentile and the upper whisker
#       extends to to the largest value within 1.5 x IQR of 
#       the 75th Percentile.
#
#       To produce a boxplot with whiskers that exetend to
#       the extreme values, set the argument range=0 in the
#       call to the boxplot() function.

stats <- summaryStats(TCE.mg.per.L ~ Well, 
	data = EPA.09.Table.9.1.TCE.df,
	combine.groups = F, digits = 3, stats.in.rows = TRUE)
stats
#        Well.1 Well.2
#N       14     13    
#Mean     0.063  0.118
#SD       0.079  0.02 
#Median   0.031  0.11 
#Min      0.004  0.099
#Max      0.25   0.17 
#NA's     1      2    
#N.Total 15     15

stats.no.ND <- summaryStats(TCE.mg.per.L ~ Well, 
	data = EPA.09.Table.9.1.TCE.df, 
	subset = Data.Qualifier != "U",
	combine.groups = F, digits = 3, stats.in.rows = TRUE)

stats.no.ND
#        Well.1 Well.2
#N       11     10    
#Mean     0.079  0.124
#SD       0.082  0.019
#Median   0.05   0.118
#Min      0.004  0.107
#Max      0.25   0.17 
#NA's     1      2    
#N.Total 12     12

windows()
box.list <- boxplot(TCE.mg.per.L ~ Well, data = EPA.09.Table.9.1.TCE.df,
	boxwex = 0.3, varwidth = T,
	ylab = "Trichloroethene (mg/L)",
	main = "Figure 9-2. Boxplots of Trichloroethene Data for Wells 1 & 2")

points(1:2, stats["Mean",], pch = 8, cex = 0.7)

TCE            <- EPA.09.Table.9.1.TCE.df$TCE.mg.per.L
Well           <- EPA.09.Table.9.1.TCE.df$Well
limits <- sapply(split(TCE, Well), 
	function(x) enorm(x, ci=T)$interval$limits)
#Warning messages:
#1: In is.not.finite.warning(x) :
#  There were 1 nonfinite values in x : 1 NA's
#2: In enorm(x, ci = T) : 1 observations with NA/NaN/Inf in 'x' removed.
#3: In is.not.finite.warning(x) :
#  There were 2 nonfinite values in x : 2 NA's
#4: In enorm(x, ci = T) : 2 observations with NA/NaN/Inf in 'x' removed.

points(1:2, limits["UCL", ], pch = 94, cex = 1.2) 
points(1:2, limits["LCL", ], pch = 118) 

abline(h = 0.23, lty = 3)
text(0.5, 0.24, "PRG = 0.23 mg/L", adj = 0, cex = 0.8)

mtext(paste("FOD =", stats.no.ND["N", ], "/", stats["N", ]), 
	at = 1:2, side = 1, line = 3)

legend(1.25, 0.07, c("outlier > 1.5 x IQR", "95% UCL", "Mean", "95% LCL"), 
	pch = c(1, 94, 8, 118), pt.cex = c(1, 1.2, 0.7, 1),
	bty = "n")

rm(TCE, Well, stats, stats.no.ND, box.list, limits)


###################################################################


# Example 9-3, pp. 9-9 to 9-12
#-----------------------------

TCE  <- EPA.09.Table.9.1.TCE.df$TCE.mg.per.L
Well <- EPA.09.Table.9.1.TCE.df$Well

#####
# Frequency Histograms
#####

##
# Variation 1: Create histograms using the hist() function
windows()
par(mfrow = c(2, 1), mar = c(5, 4, 0, 1) + 0.1, oma = c(0, 0, 3, 0))
hist(TCE[Well == "Well.1"], col = "blue", ylim = c(0, 10),
	xlab = "Trichloroethene (mg/L) in Well 1", 
	ylab = "Frequency", main = "")
hist(TCE[Well == "Well.2"], col = "blue", ylim = c(0, 8),
	xlab = "Trichloroethene (mg/L) in Well 2", 
	ylab = "Frequency", main = "")
mtext("Figure 9-3. Frequency Histrograms of Trichloroethene by Well.", 
	outer=T, line = 1, cex = 1.25, font = 2)

# Variation 2: Force the x-axis to be the same for each histogram
windows()
xlim <- range(TCE, na.rm = T)
par(mfrow = c(2, 1), mar = c(5, 4, 0, 1) + 0.1, oma = c(0, 0, 3, 0))
hist(TCE[Well == "Well.1"], col = "blue", 
	xlim = xlim, ylim = c(0, 10),
	xlab = "Trichloroethene (mg/L) in Well 1", 
	ylab = "Frequency", main = "")
hist(TCE[Well == "Well.2"], col = "blue", 
	xlim = xlim, ylim = c(0, 8),
	xlab = "Trichloroethene (mg/L) in Well 2", 
	ylab = "Frequency", main = "")
mtext("Figure 9-3. Frequency Histrograms of Trichloroethene by Well.", 
	outer=T, line = 1, cex = 1.25, font = 2)


##
# Variation 3: Create histograms using the function 
# histogram in the lattice library.

library(lattice)
windows()
histogram(~ TCE.mg.per.L | Well, data = EPA.09.Table.9.1.TCE.df, 
	type = "count", layout = c(1, 2), 
	xlab = "Trichloroethene (mg/L)", ylab = "Frequency",
	main = "Figure 9-3. Frequency Histrograms of Trichloroethene by Well.")
	


#####
# Relative Frequency Histograms
#####

##
# Variation 1: Create relative frequency histograms using the hist() function
windows()
par(mfrow = c(2, 1), mar = c(5, 4, 0, 1) + 0.1, oma = c(0, 0, 3, 0))
hist(TCE[Well == "Well.1"], freq = F, 
	col = "blue", ylim = c(0, 15), 
	xlab = "Trichloroethene (mg/L) in Well 1", 
	ylab = "Relative Frequency (%)", main = "")
hist(TCE[Well == "Well.2"], freq = F, 
	col = "blue", ylim = c(0, 25),
	xlab = "Trichloroethene (mg/L) in Well 2", 
	ylab = "Relative Frequency (%)", main = "")
mtext("Figure 9-4. Relative Frequency Histrograms of Trichloroethene by Well.", 
	outer=T, line = 1, cex = 1.25, font = 2)

##
# Variation 2: Force the x-axis to be the same for each histogram
windows()
xlim <- range(TCE, na.rm = T)
par(mfrow = c(2, 1), mar = c(5, 4, 0, 1) + 0.1, oma = c(0, 0, 3, 0))
hist(TCE[Well == "Well.1"], freq = F, col = "blue", 
	xlim = xlim, ylim = c(0, 15),
	xlab = "Trichloroethene (mg/L) in Well 1", 
	ylab = "Relative Frequency (%)", main = "")
hist(TCE[Well == "Well.2"], freq = F, col = "blue", 
	xlim = xlim, ylim = c(0, 25),
	xlab = "Trichloroethene (mg/L) in Well 2", 
	ylab = "Relative Frequency (%)", main = "")
mtext("Figure 9-4. Relative Frequency Histrograms of Trichloroethene by Well.", 
	outer=T, line = 1, cex = 1.25, font = 2)

##
# Variation 3: Create relative frequency histograms using the function 
#              histogram in the lattice library.

library(lattice)
windows()
histogram(~ TCE.mg.per.L | Well, data = EPA.09.Table.9.1.TCE.df, 
	type = "percent", layout = c(1, 2), 
	xlab = "Trichloroethene (mg/L)", ylab = "Relative Frequency (%)",
	main = "Figure 9-4. Relative Frequency Histrograms of Trichloroethene by Well.")

rm(TCE, Well, xlim)


###################################################################


# Example 9-4, pp. 9-13 to 9-14
#------------------------------

As    <- EPA.09.Table.9.3.df$Arsenic.mg.per.L
As.DQ <- EPA.09.Table.9.3.df$Arsenic.Data.Qualifier

Hg    <- EPA.09.Table.9.3.df$Mercury.mg.per.L
Hg.DQ <- EPA.09.Table.9.3.df$Mercury.Data.Qualifier

ND.ind <- rep("", length(As))
ND.ind[As.DQ != "U" & Hg.DQ != "U"] <- "both detected"
ND.ind[As.DQ == "U" & Hg.DQ == "U"] <- "both non-detects"
ND.ind[As.DQ == "U" & Hg.DQ != "U"] <- "arsenic non-detect only"
ND.ind[As.DQ != "U" & Hg.DQ == "U"] <- "mercury non-detect only"
ND.ind <- factor(ND.ind, levels = c("both detected", 
	"both non-detects", "arsenic non-detect only", 
	"mercury non-detect only"))

# Summary statistics for Arsenic, setting non-detects to detection limit
summaryStats(As)
#    N Mean  SD Median Min Max
#As 15  0.1 0.1    0.1   0 0.3

# Summary statistics for Mercury, setting non-detects to detection limit
summaryStats(Hg)
#    N Mean  SD Median Min Max NA's n.Total
#Hg 14  0.1 0.1    0.1   0 0.3    1      15


# NOTE: there is one (As,Hg) pair (0.01, 0.02) where both are detects
#       and      one (As,Hg) pair (0.01, 0.02) where both are non-detects
#       so you can't use a solid symbol for one or both or they will
#       not both show up!

lim <- range(As, Hg, na.rm = T)
pch.vec <- c(1, 4, 6, 2)
col.vec <- c("black", "red", "blue", "blue")
levels.ND.ind <- levels(ND.ind)
n.levels <- length(levels(ND.ind))
windows()
plot(As, Hg, type = "n", xlim = lim, ylim = lim, 
	xlab = "Arsenic (mg/L)", ylab = "Mercury (mg/L)",
	main = "Figure 9-5. Sctterplot of Arsenic with Mercury from Well 3")
for(i in 1:n.levels) {
	index <- ND.ind == levels.ND.ind[i]
	points(As[index], Hg[index], lwd = 2,
		pch = pch.vec[i], col = col.vec[i])
}
abline(lm(Hg ~ As))
mtext(paste("Pearson Correlation =", 
	round(cor(As, Hg, use="complete.obs"), 2)), line = 0.5)

legend(0.01, 0.3, levels.ND.ind, pch = pch.vec, 
	col = col.vec, lty = 0, lwd = 2, bty = "n")

rm(As, As.DQ, Hg, Hg.DQ, ND.ind, lim, 
	pch.vec, col.vec, levels.ND.ind, n.levels,
	i, index)


###################################################################


# Example 9-5, pp. 9-15 to 9-16
#------------------------------
	

As    <- EPA.09.Table.9.3.df$Arsenic.mg.per.L
As.DQ <- EPA.09.Table.9.3.df$Arsenic.Data.Qualifier

Hg    <- EPA.09.Table.9.3.df$Mercury.mg.per.L
Hg.DQ <- EPA.09.Table.9.3.df$Mercury.Data.Qualifier

Sr    <- EPA.09.Table.9.3.df$Strontium.mg.per.L
Sr.DQ <- EPA.09.Table.9.3.df$Strontium.Data.Qualifier


lim <- range(As, Hg, Sr, na.rm = T)
pch.vec <- c(1, 4, 6, 2)

cor.vec <- rep(as.numeric(NA), 3)

# First plot Hg vs. As

ND.ind <- rep("", length(As))
ND.ind[As.DQ != "U" & Hg.DQ != "U"] <- "both detected"
ND.ind[As.DQ == "U" & Hg.DQ == "U"] <- "both non-detects"
ND.ind[As.DQ == "U" & Hg.DQ != "U"] <- "As non-detect only"
ND.ind[As.DQ != "U" & Hg.DQ == "U"] <- "Hg non-detect only"
ND.ind <- factor(ND.ind, levels = c("both detected", 
	"both non-detects", "As non-detect only", 
	"Hg non-detect only"))

levels.ND.ind <- levels(ND.ind)
n.levels <- length(levels(ND.ind))
windows()
plot(As, Hg, type = "n", xlim = lim, ylim = lim, 
	xlab = "mg/L", ylab = "mg/L",
	main = "Figure 9-6. Coded Sctterplot of Well 3 Arsenic, Mercury, and Strontium")
for(i in 1:n.levels) {
	index <- ND.ind == levels.ND.ind[i]
	points(As[index], Hg[index], lwd = 2,
		pch = pch.vec[i], col = "blue")
}
abline(lm(Hg ~ As), lwd = 2, col = "blue")
cor.vec[1] <- round(cor(As, Hg, use="complete.obs"), 2)


# Now plot Sr vs. As

ND.ind <- rep("", length(As))
ND.ind[As.DQ != "U" & Sr.DQ != "U"] <- "both detected"
ND.ind[As.DQ == "U" & Sr.DQ == "U"] <- "both non-detects"
ND.ind[As.DQ == "U" & Sr.DQ != "U"] <- "As non-detect only"
ND.ind[As.DQ != "U" & Sr.DQ == "U"] <- "Sr non-detect only"
ND.ind <- factor(ND.ind, levels = c("both detected", 
	"both non-detects", "As non-detect only", 
	"Sr non-detect only"))

levels.ND.ind <- levels(ND.ind)
n.levels <- length(levels(ND.ind))

for(i in 1:n.levels) {
	index <- ND.ind == levels.ND.ind[i]
	points(As[index], Sr[index], lwd = 2,
		pch = pch.vec[i], col = "red")
}
abline(lm(Sr ~ As), lwd = 2, col = "red")
cor.vec[2] <- round(cor(As, Sr, use="complete.obs"), 2)


# Now plot Sr vs. Hg

ND.ind <- rep("", length(Hg))
ND.ind[Hg.DQ != "U" & Sr.DQ != "U"] <- "both detected"
ND.ind[Hg.DQ == "U" & Sr.DQ == "U"] <- "both non-detects"
ND.ind[Hg.DQ == "U" & Sr.DQ != "U"] <- "Hg non-detect only"
ND.ind[Hg.DQ != "U" & Sr.DQ == "U"] <- "Sr non-detect only"
ND.ind <- factor(ND.ind, levels = c("both detected", 
	"both non-detects", "Hg non-detect only", 
	"Sr non-detect only"))

levels.ND.ind <- levels(ND.ind)
n.levels <- length(levels(ND.ind))

for(i in 1:n.levels) {
	index <- ND.ind == levels.ND.ind[i]
	points(Hg[index], Sr[index], lwd = 2,
		pch = pch.vec[i], col = "green")
}
abline(lm(Sr ~ Hg), lwd = 2, col = "green")
cor.vec[3] <- round(cor(As, Sr, use="complete.obs"), 2)


legend(0.01, 0.3, c(paste(
	"Mercury (Y) vs. Arsenic (X), Cor =", 
		cor.vec[1]),
	paste("Strontium (Y) vs. Arsenic (X), Cor =", 
		cor.vec[2]), 
	paste("Strontium (Y) vs. Mercury (X), Cor =",
		cor.vec[3])
	),  
	col = c("blue", "red", "green"), pch = 1, 
	lty = 1, cex = 0.8, lwd = 2, bty = "n")


rm(As, As.DQ, Hg, Hg.DQ, Sr, Sr.DQ, 
	ND.ind, lim, 
	pch.vec, levels.ND.ind, n.levels,
	i, index, cor.vec)



# NOTE:  It is easier to see what is going on by using 
#        matrix scatterplots, rather than
#        a coded scatterplot.

# Matrix scatterplots with:
# non-detects plotted at detection limit, 
# duplitcate points overlaid, 
# no regression lines, and
# no special symbols for non-detects.

dum.df <- EPA.09.Table.9.3.df[, 
	c("Arsenic.mg.per.L", "Mercury.mg.per.L", "Strontium.mg.per.L")]
names(dum.df) <- c("As (mg/L)", "Hg (mg/L)", "Sr (mg/L)")

lim <- range(dum.df, na.rm = T)

windows()
pairs(dum.df, xlim = lim, ylim = lim)


# Matrix scatterplots with:
# non-detects plotted at detection limit, 
# duplitcate points jittered, 
# regression lines included, but
# no special symbols for non-detects,
# points plotted in upper panels,
# Pearson correlations shown in lower panels.

panel <- function(x, y, factor.x = 10, factor.y = 10, ...) {
	points.w.dups(x, y, 
		factor.x = factor.x, 
		factor.y = factor.y, ...)
	abline(lm(y~x))
}

panel.cor <- function(x, y, digits=2, prefix="Cor = ", 
	cex = 2 * par("cex"), ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = "pairwise.complete.obs")
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    text(0.5, 0.5, txt, cex = cex)
}


windows()
set.seed(233)
pairs(dum.df, upper.panel = panel, 
	lower.panel = panel.cor, 
	xlim = lim, ylim = lim)


rm(dum.df, panel, panel.cor, lim)


#################################################################


# Example 9-6, pp. 9-17 to 9-20
#------------------------------

windows()
qqPlot(EPA.09.Table.9.4.nickel.vec, ylim = c(0, 1000),
	plot.pos.con = 0, add.line = T,
	ylab = "Nickel Concentration (ppb)",
	main = "Figure 9.7.  Normal Probability Plot")

windows()
qqPlot(EPA.09.Table.9.4.nickel.vec, dist = "lnorm", 
	plot.pos.con = 0, add.line = T,
	ylab = "log(nickel) Concentration log(ppb)",
	main = "Figure 9.8.  Probability Plot of Log Transformed Nickel Data")






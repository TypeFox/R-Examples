## ----setup, echo=FALSE, results="hide"--------------------------------------------------
library("knitr")
opts_chunk$set(fig.align="center", fig.width=7, fig.height=7)
options(width=90)

## ---------------------------------------------------------------------------------------
library("mountainplot")
data(singer, package = "lattice")
parts <- within(singer, {
section <- voice.part
section <- gsub(" 1", "", section)
section <- gsub(" 2", "", section)
section <- factor(section)
})
# Change levels to logical ordering
levels(parts$section) <- c("Bass","Tenor","Alto","Soprano")

## ----mtn--------------------------------------------------------------------------------
require(latticeExtra) # for ecdfplot
ecdfplot(~height|section, data = parts, groups=voice.part, type='l',
         layout=c(1,4),
         main="Empirical CDF",
         auto.key=list(columns=4), as.table=TRUE)

## ---------------------------------------------------------------------------------------
mountainplot(~height|section, data = parts,
             groups=voice.part, type='l',
             layout=c(1,4),
             main="Folded Empirical CDF",
             auto.key=list(columns=4), as.table=TRUE)

## ----finish, echo=FALSE, results="asis"-------------------------------------------------
toLatex(sessionInfo(), locale=FALSE)


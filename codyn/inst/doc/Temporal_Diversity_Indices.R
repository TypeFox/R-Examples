## ----echo=FALSE----------------------------------------------------------
library(knitr)
options(digits = 3)

## ----echo=FALSE----------------------------------------------------------
library(codyn)
library(knitr)
data(collins08)
kable(head(collins08))

## ----results='asis'------------------------------------------------------
KNZ_turnover <- turnover(df = collins08, 
                       time.var = "year", 
                       species.var = "species", 
                       abundance.var = "abundance", 
                       replicate.var = "replicate")
    kable(head(KNZ_turnover))

## ----results='asis'------------------------------------------------------
KNZ_turnover_agg <- turnover(df = collins08, 
                          species.var = "species",
                          time.var = "year",
                          abundance.var = "abundance",
                          replicate.var = NA)

## ----echo=FALSE----------------------------------------------------------
    kable(head(KNZ_turnover_agg))

## ----results='asis'------------------------------------------------------
KNZ_appearance <- turnover(df = collins08, 
                        time.var = "year",
                        species.var = "species",
                        abundance.var = "abundance",
                        replicate.var = "replicate",
                        metric = "appearance")

## ----results='asis'------------------------------------------------------
KNZ_disappearance <- turnover(df = collins08,
                        time.var = "year",
                        species.var = "species",
                        abundance.var = "abundance",
                        replicate.var = "replicate",
                        metric = "disappearance")

## ----fig.width = 7, fig.height = 4, echo=FALSE---------------------------
library(ggplot2)
ggplot(KNZ_turnover, aes(x=year, y=total, color=replicate)) + geom_line(size = 2) + theme_bw()

## ----fig.width = 7, fig.height = 4, echo=FALSE---------------------------
KNZ_appearance$metric<-"appearance"
names(KNZ_appearance)[1]="turnover"

KNZ_disappearance$metric<-"disappearance"
names(KNZ_disappearance)[1]="turnover"

KNZ_appdisapp<-rbind(KNZ_appearance, KNZ_disappearance)

ggplot(KNZ_appdisapp, aes(x=year, y=turnover, color=replicate)) + geom_line(size = 2) + theme_bw() + facet_wrap(~metric) + ggtitle("Species appearances and disappearances \n Annually burned vs unburned plots, Konza
          ")

## ------------------------------------------------------------------------

## Generate some sample data 
yr = 1977:2003
sp1 = .3*sin(yr) + rnorm(length(yr), 0, .1) + -.05*yr + 150
sp2 = .2*sin(yr) + rnorm(length(yr), 0, .1) + -.01*yr + 70
sp3 = .2*sin(yr) + rnorm(length(yr), 0, .1) + .01*yr + 30

dat = data.frame(year = rep(yr, 3), 
           species = gl(3, length(yr), labels = c("sp1","sp2","sp3")),
           abundance = c(sp1, sp2, sp3))

## ---- fig.width = 7, fig.height = 4--------------------------------------
ggplot(dat, aes(year, abundance, color = species)) + 
  geom_line(size = 2) + theme_bw()

## ---- fig.width = 5, fig.height = 5--------------------------------------
ggplot(dat, aes(year, abundance, color = species)) + 
  geom_line(size = 2) + coord_polar() + theme_bw() 

## ---- fig.width = 7, fig.height = 5--------------------------------------

aggdat <- aggregate(abundance ~ species * year * replicate, 
                    data = subset(collins08, 
                                    species == "andrgera" |
                                    species == "andrscop" | 
                                    species == "poaprat"| 
                                    species == "sorgnuta"), 
                    FUN = mean)

ggplot(aggdat, aes(year, abundance, color = species)) + 
  geom_line(size = 2) + coord_polar() + theme_bw() + facet_wrap(~replicate) +
  ggtitle("Dominant species abundances \n Annually burned vs unburned plots, Konza \n")

## ------------------------------------------------------------------------
KNZ_rankshift <- rank_shift(df=collins08,
                        time.var = "year",
                        species.var = "species",
                        abundance.var = "abundance", 
                        replicate.var = "replicate")
#Select the final time point from the returned time.var_pair
KNZ_rankshift$year <- as.numeric(substr(KNZ_rankshift$year_pair, 6,9))

## ----echo=FALSE----------------------------------------------------------
kable(head(KNZ_rankshift))

## ------------------------------------------------------------------------
KNZ_rankshift_agg <- rank_shift(df = collins08,
                        time.var = "year",
                        species.var = "species",
                        abundance.var = "abundance")

## ----fig.width = 7, fig.height = 4---------------------------------------
# Create a column with the final year from the returned time.var_pair
KNZ_rankshift$year <- as.numeric(substr(KNZ_rankshift$year_pair, 6, 9))

# Plot it
ggplot(KNZ_rankshift, aes(year, MRS, color=replicate)) + 
  geom_line(size= 2) + theme_bw() 

## ------------------------------------------------------------------------
    rate.res <- rate_change(collins08,  
                    time.var= "year", 
                    species.var= "species", 
                    abundance.var= "abundance", 
                    replicate.var = "replicate")

## ---- echo = F-----------------------------------------------------------
rownames(rate.res) = NULL

## ---- fig.width = 7, fig.height = 5--------------------------------------
kable(rate.res)


## ---- fig.width = 7, fig.height = 5--------------------------------------

#Use the rate_change_interval function to generate the full data frame of distances by time lag intervals
comm.res <- rate_change_interval(collins08, 
                              time.var = "year",
                              species.var = "species",
                              abundance.var = "abundance",
                              replicate.var = "replicate")

ggplot(comm.res, aes(interval, distance, color = replicate)) + facet_wrap(~replicate) + 
  geom_point() + theme_bw() + stat_smooth(method = "lm", se = F, size = 2)



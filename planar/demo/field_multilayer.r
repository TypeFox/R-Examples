
## ----load, echo=FALSE,results='hide'-----------------------------------------
library(knitr)
library(ggplot2)
opts_chunk$set(fig.path="fieldmultilayer/",
               warning=FALSE,error=FALSE,message=FALSE,tidy=FALSE)
library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))


## ----setup, results='hide'----------------------------------------------------
library(planar)
library(ggplot2)
require(reshape2)
require(plyr)


## ----simulation----------------------------------------------------------
m <- lfief(wavelength=633, angle=44*pi/180, polarisation='p',
           thickness=c(0, 50, 0), dmax=400, res=5000,
           epsilon=list(1.5^2, -12+1i, 1.0^2), 
           displacement=TRUE)
m$L1 <- factor(m$L1, labels=c("glass", "metal", "air"))
str(m)


## ----plot, fig.width=12--------------------------------------------------

limits <- ddply(m, .(L1), summarize, 
                xmin=min(x), xmax=max(x), ymin=-Inf, ymax=Inf)

  ggplot(m) +
  facet_grid(variable~., scales="free")+
  geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, 
                fill=factor(L1)),     data=limits, alpha=0.2) +
  geom_path(aes(x, value, colour=factor(L1))) +
  scale_x_continuous("x /nm",expand=c(0,0)) +
  scale_y_continuous("LFIEF",expand=c(0,0)) +
  geom_hline(yintercept=0) +
  scale_colour_brewer("", palette="Set1")+
  scale_fill_brewer("", palette="Pastel1") +
  theme_minimal()





## ----setup, include=FALSE, cache=FALSE-----------------------------------
# set global chunk options
library(knitr)
opts_chunk$set(fig.path='figure/minimal-', fig.align='center', fig.show='hold',dev='pdf',
warning=FALSE,dev.args=list(family="Palatino"), tidy.opts=list(keep.blank.line=FALSE, width.cutoff=60))
options(replace.assign=TRUE,width=90)
#render_listings()


## ----loadnCal, include=FALSE, cache=FALSE, tidy=TRUE, echo=TRUE----------
library(kyotil)
library(MASS)
set.seed(1)


## ----corplot, include=TRUE, cache=FALSE, tidy=TRUE, echo=TRUE, eval=TRUE, fig.width=8, fig.height=8.5, fig.pos="tbp", fig.cap="ncal graphical output, drm fit."----
dat=data.frame(mvrnorm(100,c(0,0),matrix(c(1,.5,.5,1),2,2)))
names(dat)=c("x","y")
corplot(y~x, dat, main="Some Title", add.diagonal.line=TRUE)



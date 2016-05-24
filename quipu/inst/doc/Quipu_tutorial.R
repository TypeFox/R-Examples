
## ----quipu,include=FALSE,echo=FALSE--------------------------------------
options(keep.source = TRUE, width = 60)
qp <- packageDescription("quipu")


## ------------------------------------------------------------------------

if(interactive()){
  quipu::runDemo()  
}


## ------------------------------------------------------------------------
library(quipu)
library(xtable)
data("potato.quipu", package="quipu")
dat = potato.quipu


## ----echo=FALSE, results='asis'------------------------------------------
xtable(head(dat))


## ------------------------------------------------------------------------
if(interactive()) rquipu(dat)


## ----figureexample1------------------------------------------------------
rquipu(dat, a.subset = "sample.1")


## ----figureexample2------------------------------------------------------
rquipu(dat, a.subset = "sample.1", 
       col.node = c("darkred", "orange","green","black"))


## ----figureexample3------------------------------------------------------
rquipu(dat, a.subset = "sample.2", 
       node.size = c(2, 1.5, 1, 0.5))


## ----figureexample4------------------------------------------------------
rquipu(dat, a.subset = "sample.1", 
       col.node = c("darkred", "orange","darkgreen","black"),
       node.size = c(0.9, 1.2, 1.5, 1.8))


## ----figureexample5------------------------------------------------------
rquipu(dat, a.subset = "sample.1", 
       grp.brks = c(0.25, 0.5, 0.75))


## ----figureexample6------------------------------------------------------
rquipu(dat, a.subset = "sample.1",
       col.node = c("darkred", "grey80","grey80","blue"),
       node.size = c(0.9, 1.2, 1.2, 1.8),
       grp.brks = c(0.25, 0.5, 0.75)
       )


## ------------------------------------------------------------------------
data("allele.freqs", package="quipu")


## ----echo=FALSE, results='asis'------------------------------------------
xtable(head(allele.freqs), digits = 4)


## ------------------------------------------------------------------------
allele.freqs[3,3] <- 0.0001


## ----echo=FALSE, results='asis'------------------------------------------
xtable(head(allele.freqs), digits = 4)


## ----figureexample7------------------------------------------------------
rquipu(dat, a.subset = "sample.1", obs.alls.frq = allele.freqs)


## ------------------------------------------------------------------------
xx = rbind(c("STG0001", 300, 0.01), allele.freqs)
xx$marker_size = as.integer(xx$marker_size)
xx$frequency <- as.numeric(xx$frequency)


## ----echo=FALSE, results='asis'------------------------------------------
xtable(head(xx), digits = 4)


## ----figureexample8------------------------------------------------------
rquipu(dat, a.subset = "sample.1", 
            obs.alls.frq = xx, 
            obs.alls.frq.ref = "PotDb v1.0")


## ----figureexample9, fig.width=6, fig.height=2---------------------------
par( mfrow = c( 1, 3 ) )
rquipu(dat, 
       layout = "no text",
       res = c(600, 400))



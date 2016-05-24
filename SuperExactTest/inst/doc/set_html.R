## ----setup, include=FALSE------------------------------------------------
library(knitr)

## ------------------------------------------------------------------------
library("SuperExactTest")
data("eqtls")

## ------------------------------------------------------------------------
str(cis.eqtls)

## ------------------------------------------------------------------------
(length.gene.sets=sapply(cis.eqtls,length))

## ------------------------------------------------------------------------
total=18196
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))

## ------------------------------------------------------------------------
(p=sapply(0:101,function(i) dpsets(i, length.gene.sets, n=total)))

## ------------------------------------------------------------------------
common.genes=intersect(cis.eqtls[[1]], cis.eqtls[[2]], cis.eqtls[[3]],
 cis.eqtls[[4]])
(num.observed.overlap=length(common.genes))
(FE=num.observed.overlap/num.expcted.overlap)

## ------------------------------------------------------------------------
dpsets(num.observed.overlap, length.gene.sets, n=total)

## ------------------------------------------------------------------------
cpsets(num.observed.overlap-1, length.gene.sets, n=total, lower.tail=FALSE)

## ------------------------------------------------------------------------
fit=MSET(cis.eqtls, n=total, lower.tail=FALSE)
fit$FE
fit$p.value

## ------------------------------------------------------------------------
res=supertest(cis.eqtls, n=total)

## ----fig1, fig.width = 5, fig.height = 5, fig.cap = "A circular plot visualizing all possible intersections and the corresponding statistics amongst *cis*-eQTL gene sets."----
plot(res, sort.by="size")

## ----fig2, fig.width = 9, fig.height = 4, fig.cap = "A bar chart illustrating all possible intersections among *cis*-eQTL gene sets in a matrix layout."----
plot(res, Layout="landscape", degree=2:4, sort.by="size")

## ------------------------------------------------------------------------
summary(res)

## ---- eval=FALSE---------------------------------------------------------
#  ?summary.msets

## ---- eval=FALSE---------------------------------------------------------
#  write.csv(summary(res)$Table, file="summary.table.csv", row.names=FALSE)


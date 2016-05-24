## load library
require("GMD")
require(cluster)

## load data 
data(ruspini)

## heatmap on a `dist' object
x <- gdist(ruspini)
main <- "Heatmap of Ruspini data"
dev.new(width=10,height=10)
heatmap.3(x, main=main, mapratio=1) # default with a title and a map in square!
heatmap.3(x, main=main, revC=TRUE)  # reverse column for a symmetric look
heatmap.3(x, main=main, kr=2, kc=2) # show partition by predefined number of clusters

## show partition by elbow
css.multi.obj <- css.hclust(x,hclust(x))
elbow.obj <- elbow.batch(css.multi.obj,ev.thres=0.90,inc.thres=0.05)
heatmap.3(x, main=main, revC=TRUE, kr=elbow.obj$k, kc=elbow.obj$k)
heatmap.3(x, main=main, sub=sub("\n"," ",attr(elbow.obj,"description")), cex.sub=1.25,
          revC=TRUE,kr=elbow.obj$k, kc=elbow.obj$k) # show elbow info as subtitle

## side plot for every row clusters
dev.new(width=10,height=10)
expr1 <- list(quote(plot(do.call(rbind,i.x),xlab="x",ylab="y",
                         xlim=range(ruspini$x),ylim=range(ruspini$y),)))
heatmap.3(x, main=main, revC=TRUE, kr=elbow.obj$k, kc=elbow.obj$k, trace="none",
          row.data=as.list(data.frame(t(ruspini))),
          plot.row.clusters=TRUE,plot.row.clusters.list=list(expr1))

## side plot for every col clusters
dev.new(width=10,height=10)
expr2 <- list(quote(plot(do.call(rbind,i.x),xlab="x",ylab="y",
                         xlim=range(ruspini$x),ylim=range(ruspini$y),)))
heatmap.3(x, main=main, revC=TRUE, kr=elbow.obj$k, kc=elbow.obj$k, trace="none",
          col.data=as.list(data.frame(t(ruspini))),
          plot.col.clusters=TRUE,plot.col.clusters.list=list(expr2))


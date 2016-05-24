# webpage-figs.R

#png("rpart.plot-example1.png", width=600, height=600)
tree <- rpart(survived ~ ., data=ptitanic, cp=.012)
prp(tree, type=4, extra=6, faclen=0, main="rpart.plot",
    # fam.main="serif", # not supported for postscript devices (ok with png files)
    cex.main=3, col.main="slategray4", col=1, max.auto.cex=1.6)
#dev.off()

#png("rpart.plot-example2.png", width=600, height=600)
tree <- rpart(survived ~ ., data=ptitanic, cp=.012)
prp(tree, branch.type=5, yesno=FALSE, faclen=0,
    main="rpart.plot:\nbranch width shows number of observations",
    # fam.main="serif", # not supported for postscript devices (ok with png files)
    cex.main=2, col.main="slategray4", col=1, max.auto.cex=2)
#dev.off()

library(earth)
data(ozone1)
# return the given node and all its ancestors (a vector of node numbers)
path.to.root <- function(node, ancestors=NULL)
{
    if(node == 1)   # root?
        c(1, ancestors)
    else            # recurse, %/% 2 gives the parent of node
        c(node, path.to.root(node %/% 2, ancestors))
}
fit.oz <- rpart(O3~., data=ozone1)
node <- 22 # 22 is our chosen node, arbitrary for this example
path <- path.to.root(node)
nodes <- as.numeric(row.names(fit.oz$frame))
cols <- ifelse(nodes %in% path, 1, "slategray4")
lwds <- ifelse(nodes %in% path, 2, 1)
lty  <- ifelse(nodes %in% path, 1, 2)
nn.box.cols <- ifelse(nodes %in% path, 1, "slategray4")
#png("rpart.plot-example3.png", width=600, height=600)
prp(fit.oz, type=2, clip.right.labs=F, nn=TRUE, tweak=1.4,
   nn.cex=1.2, nn.border=0, nn.col="white",
   main=paste("rpart.plot:\npath to node", node), cex.main=2, lwd=lwds, digits=4,
   col=cols, branch.col=cols, split.col=cols, nn.box=nn.box.cols, yesno=FALSE)
#dev.off()

#png("rpart.plot-example4.png", width=600, height=600)
old.bg <- par(bg="gray20") # this doesn't seem to work for postscript files, ok with png files
iris.tree <- rpart(Species~., data=iris)
prp(iris.tree, type=0, extra=104, main="Fisher's Iris Data",
    under=TRUE, yesno=FALSE, fallen.leaves=TRUE, branch=.2, lt=" < ",
    varlen=0, faclen=0, branch.col="wheat", split.col="wheat", under.col="wheat",
    col=c("orangered", "orange", "cyan1")[iris.tree$frame$yval], xcompact=F, ycompact=F,
    cex.main=3, col.main="white", max.auto.cex=2, split.font=1, font=2)
par(bg=old.bg)
#dev.off()

# user-manual-figs.R: figures and code in the rpart.plot vignette

library(rpart.plot)
library(earth)
data(ptitanic)
data(ozone1)
par(mfrow=c(3,3))
par(mgp = c(1.5, .5, 0))

#--- front.R ---

a <- rpart(O3~., data=ozone1, cp=.024)
y <- a$frame$yval
cols <- c("#EE8888", "#EEEE88", "#88CC88")
cols <- ifelse(y > 20, cols[1], ifelse(y < 15, cols[3], cols[2]))
prp(a, main="[front] An Example",
    box.col=cols, type=4, fallen=T, branch=.3, round=0, leaf.round=9,
    clip.right.labs=F, under.cex=1,
    prefix="ozone\n", branch.col="gray", branch.lwd=2,
    extra=101, under=T, lt=" < ", ge=" >= ", cex.main=1.5)

#--- example.R ---

tree <- rpart(survived~., data=ptitanic)

tree <- rpart(survived~., data=ptitanic, cp=.02)

prp(tree, main="[example] default prp\n(type = 0)")

prp(tree, main="type = 4, extra = 6\n", type=4, extra=6, cex.main=1.2)

cols <- ifelse(tree$frame$yval == 1, "darkred", "green4")
                         # green if survived

prp(tree, main="assorted arguments",
    branch.lwd=2,        # gyp for user manual
    extra=106,           # display prob of survival and percent of obs
    nn=TRUE,             # display the node numbers
    fallen.leaves=TRUE,  # put the leaves on the bottom of the page
    branch=.5,           # change angle of branch lines
    faclen=0,            # do not abbreviate factor levels
    trace=1,             # print the automatically calculated cex
    shadow.col="gray",   # shadows under the leaves
    branch.lty=3,        # draw branches using dotted lines
    split.cex=1.2,       # make the split text larger than the node text
    split.prefix="is ",  # put "is " before split text
    split.suffix="?",    # put "?" after split text
    col=cols, border.col=cols,   # show proportions greater than .5 in green
    split.box.col="lightgray",   # lightgray split boxes (default is white)
    split.border.col="darkgray", # darkgray border on split boxes
    split.round=.5)              # round the split box corners somewhat

# the old way for comparison
plot(tree, uniform=TRUE, compress=TRUE, branch=.2)
text(tree, use.n=TRUE, cex=.9, xpd=NA) # cex is a guess, depends on your window size
title("plot.rpart for comparison", cex=.9)

#--- type.R ---

a <- rpart(survived~., data=ptitanic, control=list(cp=.02))
prp(a, faclen=0, type=0, main="[type] type = 0\n(default)\n")
prp(a, faclen=0, type=1, main="type = 1\nlabel all nodes\n(like text.rpart all=TRUE)")
prp(a, faclen=0, type=2, main="type = 2\nsplit labels below node labels\n")
prp(a, faclen=0, type=3, main="type = 3\nleft and right split labels\n")
prp(a, faclen=0, type=4, main="type = 4\nlike type=3 but with interior labels\n(like text.rpart fancy=TRUE)")

#--- anova-extra.R ---

a <- rpart(Volume~., data=trees, cp=.5)
prp(a, clip.right.labs=F, split.cex=.8, tweak=1.3, type=4, extra=0,   yesno=F, faclen=0, under=T, cex.main=1.3, main="[anova-extra] extra = 0\ndefault\n")
prp(a, clip.right.labs=F, split.cex=.8, tweak=1.3, type=4, extra=1,   yesno=F, faclen=0, under=T, cex.main=1.3, main="extra = 1\nnbr of obs\n")
prp(a, clip.right.labs=F, split.cex=.8, tweak=1.3, type=4, extra=100, yesno=F, faclen=0, under=T, cex.main=1.3, main="extra = 100\npercentage of obs\n")
prp(a, clip.right.labs=F, split.cex=.8, tweak=1.3, type=4, extra=101, yesno=F, faclen=0, under=T, cex.main=1.3, main="extra = 101\nnbr and percentage\nof obs")

#--- class-extra.R ---

a <- rpart(survived~., data=ptitanic, cp=.03)
prp(a, under.cex=1, type=4, faclen=0, under=T, extra=0,   cex.main=1, main="[class-extra] extra = 0\ndefault\n")
prp(a, under.cex=1, type=4, faclen=0, under=T, extra=1,   cex.main=1, main="extra = 1\nnbr of obs per class\n")
prp(a, under.cex=1, type=4, faclen=0, under=T, extra=2,   cex.main=1, main="extra = 2\nclass rate\n")
prp(a, under.cex=1, type=4, faclen=0, under=T, extra=3,   cex.main=1, main="extra = 3\nmisclass rate\n")
prp(a, under.cex=1, type=4, faclen=0, under=T, extra=4,   cex.main=1, main="extra = 4\nprob per class\n(leaf probs sum to 1)")
prp(a, under.cex=1, type=4, faclen=0, under=T, extra=5,   cex.main=1, main="extra = 5\nprob per class,\nfitted class not displayed")
prp(a, under.cex=1, type=4, faclen=0, under=T, extra=6,   cex.main=1, main="extra = 6\nprob of 2nd class\n(useful for binary responses)")
prp(a, under.cex=1, type=4, faclen=0, under=T, extra=7,   cex.main=1, main="extra = 7\nprob of 2nd class,\nfitted class not displayed")
prp(a, under.cex=1, type=4, faclen=0, under=T, extra=8,   cex.main=1, main="extra = 8\nprob of fitted class\n")
prp(a, under.cex=1, type=4, faclen=0, under=T, extra=9,   cex.main=1, main="extra = 9\noverall prob\n(sum over all leaves is 1)")
prp(a, under.cex=1, type=4, faclen=0, under=T, extra=100, cex.main=1, main="extra = 100\npercent of obs\n")
prp(a, under.cex=1, type=4, faclen=0, under=T, extra=106, cex.main=1, main="extra = 106\nprob of 2nd class and\npercent of obs")

#--- fraction.R ---

tree <- rpart(survived~., data=ptitanic, cp=.02)
                         # cp=.02 because want small tree for demo
prp(tree, extra=7, prefix="fraction\n", tweak=1.2, main="[fraction]")


#--- dev.R ---

tree <- rpart(survived~., data=ptitanic, cp=.02)
                         # cp=.02 because want small tree for demo
my.node.fun <- function(x, labs, digits, varlen)
{
    paste("dev", x$frame$dev)
}
prp(tree, node.fun=my.node.fun, tweak=1.2, main="[dev]")


#--- standard-plus-dev.R ---

tree <- rpart(survived~., data=ptitanic, cp=.02)
my.node.fun2 <- function(x, labs, digits, varlen)
{
    paste(labs, "\ndev", x$frame$dev)
}
prp(tree, extra=6, node.fun=my.node.fun2, tweak=1.2, main="[standard-plus-dev]")


#--- double-newline-dev.R ---

tree <- rpart(survived~., data=ptitanic, cp=.02)
my.node.fun3 <- function(x, labs, digits, varlen)
{
    # identical to my.node.fun2 above but use \n\n instead of \n
    paste(labs, "\n\ndev", x$frame$dev)
}
prp(tree, extra=6, node.fun=my.node.fun3, tweak=1.2, main="[double-newline-dev]")

#--- prefix-col.R ---

tree <- rpart(survived~., data=ptitanic, cp=.02)
prp(tree, extra=6, tweak=1.2, main="prefix-col",
    box.col=ifelse(tree$frame$yval == 2, "palegreen3", "pink"))

#--- movie.R ---

# omitted, tested in test.prp.R

#--- depth-first-tree-construction.R ---

old.par <- par(no.readonly=TRUE)
par(mfrow=c(4,5))
data(ptitanic)
tree1 <- rpart(survived~., data=ptitanic, cp=.012)
par(mfrow=c(4,3), mar = c(0, 0, 2, 0))
for(iframe in 1:nrow(tree1$frame)) {
    cols <- ifelse(1:nrow(tree1$frame) <= iframe, "black", "gray")
    dev.hold()     # hold screen output to prevent flashing
    prp(tree1, compress=FALSE, faclen=3, varlen=-3,
       col=cols, branch.col=cols, split.col=cols,
       main=sprintf("\n%d", iframe), yesno=FALSE)
    dev.flush()
    # Sys.sleep(1) # wait one second
}
par(old.par)

#--- path-to-root.R ---

tree <- rpart(survived~., data=ptitanic, cp=0.02)

# return the given node and all its ancestors (a vector of node numbers)
path.to.root <- function(node)
{
    if(node == 1)   # root?
        node
    else            # recurse, %/% 2 gives the parent of node
        c(node, path.to.root(node %/% 2))
}

node <- 11          # 11 is our chosen node, arbitrary for this example
nodes <- as.numeric(row.names(tree$frame))
cols <- ifelse(nodes %in% path.to.root(node), "black", "darkgray")
prp(tree, nn=TRUE, col=cols, branch.col=cols, split.col=cols, nn.col=cols, tweak=1.2, main="[path-to-root]")

#--- gray-background.R ---

old.bg <- par(bg="gray50")
iris.tree <- rpart(Species~., data=iris)
prp(iris.tree, type=0, extra=8, main="[gray-background]",
    under=TRUE, yesno=FALSE,
    varlen=0, faclen=0, tweak=1.4,
    col=c("orangered", "orange", "wheat")[iris.tree$frame$yval])
par(bg=old.bg)

#--- snip.R ---

tree <- rpart(survived~., data=ptitanic, cp=.012)
new.tree <- prp(tree, main="[snip-part1]")$obj # interactively trim the tree
prp(new.tree, main="[snip-part2]")             # display the new tree

#--- heat-tree.R ---

par(mfrow=c(2,2))
heat.tree <- function (tree, low.is.green=FALSE, ...) {
    y <- tree$frame$yval
    if(low.is.green)
        y <- -y
    max <- max(y)
    min <- min(y)
    cols <- rainbow(99, end=.36)[
                ifelse(y >  y[1], (y-y[1]) * (99-50) / (max-y[1])+50,
                                  (y-min)  * (50-1)  / (y[1]-min)+1)]
    prp(tree, branch.col=cols, box.col=cols, ...)
}
data(ptitanic)
tree <- rpart(age ~ ., data=ptitanic)
heat.tree(tree, type=4, varlen=0, faclen=0, fallen.leaves=TRUE)
heat.tree(tree, low.is.green=TRUE)

#--- compress.R ---

par(mfrow=c(3,3))
a <- rpart(survived~., data=ptitanic, cp=.01)
prp(a, compress=FALSE, ycompress=FALSE, main="[compress]\ncompress=FALSE\nycompress=FALSE", cex.main=.9, trace=1)
text(.2, -.12, "calculated cex:\n0.44", cex=1, xpd=NA, font=3, col=rgb(.1,.2,.4))
prp(a, compress=TRUE, ycompress=FALSE, main="\ncompress=TRUE (default)\nycompress=FALSE", cex.main=.9, trace=1)
text(.2, .15, "calculated cex:\n0.69", cex=1, xpd=NA,  font=3, col=rgb(.1,.2,.4))
prp(a, compress=TRUE, ycompress=TRUE, main="\ncompress=TRUE (default)\nycompress=TRUE (default)", cex.main=.9, trace=1)
text(.2, .2, "calculated cex:\n0.82", cex=1, xpd=NA,   font=3,  col=rgb(.1,.2,.4))

#--- compact.R ---

plot(1,1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
plot(1,1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
plot(1,1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")

a <- rpart(survived~., data=ptitanic, cp=.1)
prp(a, ycompact=FALSE, xcompact=FALSE, faclen=0, trace=1, extra=6,
    main="[compact] xcompact=FALSE\nycompact=FALSE", cex.main=1, split.border.col=0, tweak=1.2)
rect(-.05, -.1, 1.05, 1.2, xpd=NA)
# par(mfg=c(1,1)) # for plotting on top of each other to see alignment
prp(a, faclen=0, trace=1, extra=6,
    main="\ndefault:\nxcompact=TRUE\nycompact=TRUE", cex.main=1, split.border.col=0, tweak=1.2)
rect(-.76, -1.96, 1.76, 3.35, xpd=NA) # , border=2)

#--- compatibility example ---

fit <- rpart(Kyphosis ~ Age + Number + Start, data=kyphosis) # from example(rpart)
oldpar <- par(mfrow=c(1,2), xpd=NA) # side by side comparison
plot(fit)
text(fit, use.n=TRUE)
prp(fit, extra=1, uniform=F, branch=1, yesno=F, border.col=0, xsep="/")
par(oldpar)

#--- split.fun example ---

a9 <- rpart(Price/1000 ~ Mileage + Type + Country, cu.summary)
split.fun <- function(x, labs, digits, varlen, faclen)
{
    gsub(" = ", ":\n", labs)
}
prp(a9, extra=100, under=T, yesno=F, split.fun=split.fun, main="a9 (split.fun example)")

#--- utilities.R ---

tree <- rpart(survived~., data=ptitanic, cp=.02)
print(as.numeric(row.names(tree$frame))) # node numbers in the order they appear in frame

node <- 2
print(node %/% 2)                               # parent of node

print(inode <- match(node, nodes))              # row index of node in frame

print(is.leaf <- tree$frame$var == "<leaf>")    # logical vec, indexed on row in frame

print(nodes[is.leaf])                           # the leaf node numbers

print(is.left <- nodes %% 2 == 0)               # logical vec, indexed on row in frame

print(ifelse(is.left, nodes+1, nodes-1))        # siblings of nodes

get.children <- function(node)          # node and all its children
    if(is.leaf[match(node, nodes)]) {
        node
    } else
        c(node,
          get.children(2 * node),       # left child
          get.children(2 * node + 1))   # right child

print(get.children(2))

#--- split-label12.R ---

tree.split.label12 <- rpart(Price/1000 ~ Mileage + Type + Country, cu.summary)
split.fun <- function(x, labs, digits, varlen, faclen)
{
    gsub(" = ", ":\n", labs)
}
prp(tree.split.label12, extra=1, branch=1, split.border.col=1, branch.col=0, split.yspace=1.5, main="split-label12")
prp(tree.split.label12, extra=1, branch=1, split.border.col=1, branch.col=0, , split.yspace=1.5, split.fun=split.fun)
prp(tree.split.label12, extra=100, under=T, yesno=F, split.fun=split.fun, tweak=.8)

#--- split-label13.R ---

tree.split.label13 <- rpart(Price/1000 ~ Mileage + Type + Country, cu.summary)
split.fun <- function(x, labs, digits, varlen, faclen)
{
    # replace commas with spaces (needed for strwrap)
    labs <- gsub(",", " ", labs)
    for(i in 1:length(labs)) {
        # split labs[i] into multiple lines
        labs[i] <- paste(strwrap(labs[i], width=25), collapse="\n")
    }
    labs
}
prp(tree, cex=1.05)
prp(tree.split.label13, split.fun=split.fun, main="split-label13")

# test.prp.R

library(rpart.plot)
data(ptitanic)
library(earth)
data(ozone1)
sessionInfo()
ititanic <- ptitanic
ititanic$survived <- as.integer(ititanic$survived == "survived")

if(!interactive())
    postscript(paper="letter", fonts=c("Helvetica", "NewCenturySchoolbook"))

example(rpart.plot)
example(prp)
print(citation("rpart.plot"))

# test format0 and formatf

x <- c(1.2345, 1.6, 1.23456, 12.345, 124.56,
       123, 123.456789012345, 1234, 9999, 12345, 123456, 1.234e6, 1.234e7,
       .123, .0123,
       1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6,
       1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6,
       .999, .999e1, .999e2, .999e3, .999e4, .999e5, .999e6,
       .999e-1, .999e-2, .999e-3, .999e-4, .999e-5, .999e-6)

x1 <- c(x, 1.2345e300, 1.2345e-300)

x  <- c(NA, Inf, -Inf, 0, x, -x)
x1 <- c(NA, Inf, -Inf, 0, x1, -x1)

options(digits=7)

cat("\n--- format0 ---\n")
s4  <- rpart.plot:::format0(x1, digits=4)
s2  <- rpart.plot:::format0(x1, digits=2)
s1  <- rpart.plot:::format0(x1, digits=1)
s0  <- rpart.plot:::format0(x1, digits=0)
sm3 <- rpart.plot:::format0(x1, digits=-3)
tab <- data.frame("digits=4"=s4, "digits=2"=s2, "digits=1"=s1, "digits=0"=s0, "digits=-3"=sm3, check.names=F)
row.names(tab) <- format(x1, digits=15)
print(tab)

cat("\n--- formatf ---\n")
s4  <- rpart.plot:::formatf(x, digits=4)
s2  <- rpart.plot:::formatf(x, digits=2)
s1  <- rpart.plot:::formatf(x, digits=1)
s0  <- rpart.plot:::formatf(x, digits=0)
tab <- data.frame("digits=4"=s4, "digits=2"=s2, "digits=1"=s1, "digits=0"=s0,
                  check.names=F)
row.names(tab) <- format(x, digits=15)
print(tab)

cat("\n--- formatf strip.leading.zeros ---\n")
s4  <- rpart.plot:::formatf(x, digits=4, strip.leading.zeros=TRUE)
s2  <- rpart.plot:::formatf(x, digits=2, strip.leading.zeros=TRUE)
s1  <- rpart.plot:::formatf(x, digits=1, strip.leading.zeros=TRUE)
s0  <- rpart.plot:::formatf(x, digits=0, strip.leading.zeros=TRUE)
tab <- data.frame("digits=4"=s4, "digits=2"=s2, "digits=1"=s1, "digits=0"=s0,
                  check.names=F)
row.names(tab) <- format(x, digits=15)
print(tab)

# examples from the vignette

fit <- rpart(survived~., data=ititanic)
cols <- ifelse(fit$frame$yval > .5, "palegreen", "pink")
par(mfrow=c(2,2))
prp(fit, box.col=cols, main="Page 3", prefix="probability\n", trace=1)

fit <- rpart(survived~., data=ititanic)
cp <- sort(unique(fit$frame$complexity))[4:5] # just do 2, for a quicker test
for(i in 1:length(cp)) {
    col <- ifelse(fit$frame$complexity >= cp[i], 1, "gray")
    lwd <- ifelse(fit$frame$complexity >= cp[i], 2, 1)
    prp(fit, type=1, col=col, branch.col=col, lwd=lwd,
           sub=sprintf("movie %g", i), col.s=2, trace=1)
}

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
prp(fit.oz, type=4, clip.right.labs=F, nn=TRUE, trace=3, # some niceties
   main=paste("Path to node", node), col.m=3, lwd=lwds, digits=4,
   col=cols, branch.col=cols, split.col=cols, nn.col=cols)

my.labs <- function(x, labs, digits, varlen)
{
    sprintf("ozone %.3g\ndev %.1f", x$frame$yval, x$frame$dev)
}
data(ozone1)
fit <- rpart(O3~., data=ozone1)
par(mfrow=c(2,2))
prp(fit, node.fun=my.labs, main="Page 4: my.labs", trace=1)

my.labs2 <- function(x, labs, digits, varlen)
{
    sprintf("%s\n%.2g",
            ifelse(x$frame$yval > .5, "survived", "died"),
            x$frame$yval)
}
fit2 <- rpart(survived~., data=ititanic)
prp(fit2, node.fun=my.labs2, main="my.labs2", trace=1)

my.labs3 <- function(x, labs, digits, varlen) # use passed in labs
{
    sprintf("%s\n\ndev %.1f", labs, x$frame$dev)
}
prp(fit2, node.fun=my.labs3, main="my.labs3\nextra=100", trace=1, extra=100, under=T)

# commented out for rpart.plot version 1.4-0 (user mode no longer supported)
# fit.user <- rpart(survived~., data=ptitanic, cp=.02)
# fit.user$method <- "user"
# fit.user$functions$text <- function (yval, dev, wt, ylevel, digits, n, use.n)
# {
#     nclass <- (ncol(yval) - 1L)/2
#     group <- yval[, 1L]
#     counts <- yval[, 1L + (1L:nclass)]
#     if (!is.null(ylevel))
#         group <- ylevel[group]
#     temp1 <- format(counts)
#     if (nclass > 1)
#         temp1 <- apply(matrix(temp1, ncol = nclass), 1, paste, collapse = "/")
#     cat("use.n=", use.n, "\n")
#     if (use.n)
#         out <- paste(group, "!\n", temp1, sep = "")
#     else
#         out <- format(group)
#     return(out)
# }
# prp(fit.user, node.fun=my.labs3, main="method=user\nmy.labs3 extra=100", trace=1, extra=100, under=T, prefix="result: ")

a20 <- rpart(survived~., data=ptitanic, control=list(cp=.02))
par(mfrow=c(3,3))

boxes.include.gap <- FALSE

prp(a20, type=4,
    main="Page 5: box positioning, type=4\n(1) extra=0\nprefix=0 suffix=0 split.suffix=0\nsplit.cex=1\n",
    cex.main=.9,
    under=F,
    extra=0,
    split.cex=1,
    faclen=0, trace=1,
    boxes.include.gap=boxes.include.gap)

prp(a20, type=4,
    main="(2) extra=1\nprefix=0 suffix=0 split.suffix=0\nsplit.cex=1\n",
    cex.main=.9,
    under=F,
    extra=1,
    split.cex=1,
    faclen=1, trace=1,
    branch=1,
    boxes.include.gap=boxes.include.gap)

prp(a20, type=4,
    main="(3) extra=1\nprefix=1 suffix=0 split.suffix=0\nsplit.cex=1\n",
    cex.main=.9,
    under=F,
    extra=1,
    prefix="prefix ",
    suffix="",
    split.suffix="",
    split.cex=1,
    faclen=0, trace=1,
    split.border.col=1,
    boxes.include.gap=boxes.include.gap)

prp(a20, type=4,
    main="(4) extra=4\nprefix=1 suffix=1 split.suffix=1\nsplit.cex=1\n",
    cex.main=.9,
    under=F,
    extra=4,
    prefix="prefix ",
    suffix="\nsuffix",
    split.suffix="\nsplit.suffix",
    split.cex=1,
    faclen=0, trace=1,
    boxes.include.gap=boxes.include.gap)

prp(a20, type=4,
    main="(5) extra=5\nprefix=0 suffix=0 split.suffix=0\nsplit.cex=1.4\n",
    cex.main=.9,
    under=F,
    extra=5,
    split.cex=1.4,
    faclen=0, trace=1,
    split.border.col=1,
    boxes.include.gap=boxes.include.gap)

prp(a20, type=4,
    main="(6) extra=101\nprefix=0 suffix=0 split.suffix=0\nsplit.cex=.7\n",
    cex.main=.9,
    under=F,
    extra=101,
    split.cex=.7,
    faclen=0, trace=1,
    branch=1,
    boxes.include.gap=boxes.include.gap)

prp(a20, type=4,
    main="(7) extra=1\nprefix=1 suffix=1 split.suffix=1\nsplit.cex=1.4\n",
    cex.main=.9,
    under=F,
    extra=1,
    prefix="prefix ",
    suffix="\nsuffix",
    split.suffix="\nsplit.suffix",
    split.cex=1.4,
    faclen=1, trace=1,
    split.border.col=1,
    boxes.include.gap=boxes.include.gap)

prp(a20, type=4,
    main="(8) extra=0\nprefix=1 suffix=1 split.suffix=1\nsplit.cex=.7\n",
    cex.main=.9,
    under=F,
    extra=0,
    prefix="prefix ",
    suffix="\nsuffix",
    split.suffix="\nsplit.suffix",
    split.cex=.7,
    faclen=0, trace=1,
    boxes.include.gap=boxes.include.gap)

# TODO split.yshift indexing is confusing
prp(a20, type=4, trace=2, nn=T, split.border.col=1,
    main="(9) manual yshift and split.yshift\n",
    prefix=c("up ", "", "up ", "", "", "", ""),
    yshift=c(2, 0, 2, 0, 0, 0, 0),
    branch=1,
    split.prefix=c("", "", "down ", "", "", "", ""),
    split.yshift=c(0, -3, 0, 0, 0, 0, 0))


a21 <- rpart(survived~., data=ptitanic, control=list(cp=.02))
par(mfrow=c(3,3))

prp(a21, type=1,
    main="Page 6: type=1\n(1) extra=0\nprefix=0 suffix=0 split.suffix=0\nsplit.cex=1\n",
    cex.main=.9,
    under=F,
    extra=0,
    split.cex=1,
    faclen=0, trace=1,
    boxes.include.gap=boxes.include.gap)

prp(a21, type=1,
    main="(2) extra=1\nprefix=0 suffix=0 split.suffix=0\nsplit.cex=1\n",
    cex.main=.9,
    under=F,
    extra=1,
    split.cex=1,
    faclen=1, trace=1,
    branch=1,
    boxes.include.gap=boxes.include.gap)

prp(a21, type=1,
    main="(3) extra=1\nprefix=1 suffix=0 split.suffix=0\nsplit.cex=1\n",
    cex.main=.9,
    under=F,
    extra=1,
    prefix="prefix ",
    suffix="",
    split.suffix="",
    split.cex=1,
    faclen=0, trace=1,
    boxes.include.gap=boxes.include.gap)

prp(a21, type=1,
    main="(4) extra=1\nprefix=1 suffix=1 split.suffix=1\nsplit.cex=1\n",
    cex.main=.9,
    under=F,
    extra=1,
    prefix="prefix ",
    suffix="\nsuffix",
    split.suffix="\nsplit.suffix",
    split.cex=1,
    faclen=1, trace=1,
    boxes.include.gap=boxes.include.gap)

prp(a21, type=1,
    main="(5) extra=0\nprefix=0 suffix=0 split.suffix=0\nsplit.cex=1.4\n",
    cex.main=.9,
    under=F,
    extra=0,
    split.cex=1.4,
    faclen=0, trace=1,
    branch=1,
    boxes.include.gap=boxes.include.gap)

prp(a21, type=1,
    main="(6) extra=1\nprefix=0 suffix=0 split.suffix=0\nsplit.cex=.7\n",
    cex.main=.9,
    under=F,
    extra=1,
    split.cex=.7,
    faclen=1, trace=1,
    split.border.col=1,
    boxes.include.gap=boxes.include.gap)

prp(a21, type=1,
    main="(7) extra=1\nprefix=1 suffix=1 split.suffix=1\nsplit.cex=1.4\n",
    cex.main=.9,
    under=F,
    extra=1,
    prefix="prefix ",
    suffix="\nsuffix",
    split.suffix="\nsplit.suffix",
    split.cex=1.4,
    faclen=0, trace=1,
    branch=1,
    split.border.col=1,
    boxes.include.gap=boxes.include.gap)

prp(a21, type=1,
    main="(8) extra=0\nprefix=1 suffix=1 split.suffix=1\nsplit.cex=.7\n",
    cex.main=.9,
    under=F,
    extra=0,
    prefix="prefix ",
    suffix="\nsuffix",
    split.suffix="\nsplit.suffix",
    split.cex=.7,
    faclen=1, trace=1,
    split.border.col=1,
    boxes.include.gap=boxes.include.gap)

prp(a21, type=1, trace=2, nn=0, ni=0,
    main="(9) manual yshift and split.yshift\n",
    split.border.col=1,
    cex.main=.9,
    prefix=c("up\n", "", "up\n", "up\n", "", "", ""),
    yshift=c(3, 0, 3, 3, 0, 0, 0),
    split.prefix=c("", "", "down\n", "", "", "", ""),
    split.yshift=c(0, -1, 0, 0, 0, 0, 0))


par(mfrow=c(3,3))
a8 <- rpart(survived~., data=ptitanic, control=list(cp=.02))
prp(a8, type=2, main="Page 7: type=2")
prp(a8, type=2, extra=4, main="extra=4")
prp(a8, type=2, extra=104, main="extra=104")

prp(a8, type=2, extra=0,   under=T, main="extra=0,   under=T")
prp(a8, type=2, extra=4,   under=T, main="extra=4,   under=T")
prp(a8, type=2, extra=104, under=T, main="extra=104, under=T")

prp(a8, type=2, extra=104, under=T, under.cex=.6, main="extra=104, under=T\nunder.cex=.6")
prp(a8, type=2, extra=104, under=T, under.cex=1, main="extra=104, under=T\nunder.cex=1.2")
prp(a8, type=2, extra=104, under=T, split.border.col=1, border.col=0, main="extra=104, under=T\nsplit.border=1, border=0")

par(mfrow=c(2,2))
fit3 <- rpart(survived~., data=ititanic, control=list(cp=.002))
prp(fit3, trace=3, nn=0, faclen=0, prefix="prob ", main="Page 8: ycompress")
prp(fit3, extra=100, trace=3, nn=TRUE, faclen=0, fallen.leaves=TRUE, main="fallen leaves")
prp(fit3, type=4, trace=3, nn=TRUE, clip.right.labs=0, split.border.col=1, main="type=4 ")
# use prefix below to force shifting of fallen leaves to test shifter
prp(fit3, type=4, branch=.5, extra=1, under=TRUE, trace=3, nn=FALSE, fallen.leaves=1, prefix="probability ", main="type=4, fallen leaves")

par(mfrow=c(2,2))
fit7 <- rpart(survived ~ ., data=ptitanic, cp=.01)
# this was wrong until I added check that a shift doesn't move nodes above the nodes for the level above
prp(fit7, extra=1, branch=1, trace=3, nn=1, main="Page 9: more ycompress")

par(mfrow=c(2,2))
prp(fit2, prefix=ifelse(fit2$frame$yval > .5, "survived\n", "died\n"), main="Page 10: miscellaneous 1",
    fam.main="NewCenturySchoolbook", cex.main=1.3, trace=1,
    border.col=0, split.border.col="steelblue3")
# test long names and big and small numbers
ptitanic1 <- ptitanic
ptitanic1$sibsp1234567890 <- 1e3 * ptitanic1$sibsp
ptitanic1$sibsp <- NULL
ptitanic1$age <- 1e-5 * ptitanic1$age
ptitanic1$parch <- 1e7 * ptitanic1$parch
fit2 <- rpart(survived~., data=ptitanic1)
prp(fit2, faclen=0, digits=4, trace=1,
    border.col=NA, split.border.col="steelblue3", split.round=1)

# test small tree, also tests xcompact and ycompact
fit.small <- rpart(survived~., data=ptitanic1, , control=list(cp=.1))
prp(fit.small, extra=100, faclen=0, main="small tree", trace=1)

par(mfrow=c(2,2))
fit4 <- rpart(survived~., data=ititanic, method="class", control=list(cp=.02))
prp(fit4, trace=2, cex=.8, tweak=1.1, main="Page 11: miscellaneous 2, tweak",
       xflip=TRUE, yflip=TRUE, type=1, extra=100,  yesno=FALSE)
# TODO wanna include family below, but postscript giving me grief
fit4.strange.method <- fit4
fit4.strange.method$method <- "unknown.method"
prp(fit4.strange.method, main="left=FALSE, fonts, user method", left=FALSE, font=c(1,2,3), split.cex=c(1, 1.2), branch=.5, trace=1, extra=1)
prp(fit4, main="uniform=FALSE", uniform=FALSE, trace=1)
data(ozone1)
fit.oz1 <- rpart(O3~., data=ozone1)
obj <- prp(fit.oz1, main="digits=7", digits=7, trace=1)
cat("obj returned by prp:\n")
print(obj)

# test extra and faclen etc. on anova model
a1 <- rpart(survived~., data=ititanic, control=list(cp=.03))
par(mfrow=c(3, 3))
plot(a1, unif=TRUE, branch=.3, main="Page 12: anova\nwith extra faclen etc.\n")
text(a1, fancy=T, fwidth=.35, fheight=0.3, use.n=TRUE, all=T, digits=3, xpd=NA, pretty=0)
prp(a1, extra=0, faclen=-3, varlen=2,             type=1, main="extra=0", trace=1)
prp(a1, extra=1, faclen=1,  varlen=-2,              type=4, main="extra=1", trace=1)
prp(a1, extra=100, faclen=3,             digits=3,  type=4, clip.right.labs=FALSE, facsep=" or ", main="extra=100", trace=1)
# plot(1, 1, type="n", axes=FALSE, xlab="", ylab="") # blank
# test xflip and left (note: left=FALSE cannot be used with type=4)
prp(a1, main="\n\nxflip",                xflip=TRUE,             extra=101, faclen=0, trace=1)
prp(a1, main="type=4",                                type=4, extra=101, faclen=0, trace=1)
prp(a1, main="xflip type=4",              xflip=TRUE, type=4, extra=101, faclen=0, trace=1)
prp(a1, main="xflip type=4 clip.right=FALSE",   xflip=TRUE, type=4, extra=101, faclen=0, clip.right=FALSE, trace=1)

# test extra and faclen etc. on class model
a2 <- rpart(survived~., data=ptitanic, control=list(cp=.02))
par(mfrow=c(3, 2))
plot(a2, unif=TRUE, branch=.3, main="Page 13: class")
text(a2, use.n=TRUE, all=T, digits=3, xpd=NA, pretty=0)
prp(a2, extra=0, eq=" eq ", lt=" lt ", ge=" ge ", facsep="|", xsep="/",
    type=4, main="extra=0", trace=3, split.border.col=1)
prp(a2, extra=1, type=1,  xsep=", ", main="extra=1", trace=1)
prp(a2, extra=100, type=3, clip.right.labs=FALSE, main="extra=100", trace=3, split.border.col=1, ycompress.cex=1)
prp(a2, extra=4, type=0, main="extra=4", faclen=0, trace=1, under=TRUE, col=2)
prp(a2, extra=104, type=0, main="extra=104", faclen=0, trace=1)

old.par <- par(mfrow=c(8,4), mar = c(4, 3, 2, 1), mgp = c(1.5, .5, 0))

a4 <- rpart(survived~., data=ptitanic, cp=.03)
plot(a4, unif=T, branch=.3); text(a4, use.n=1, cex=1, xpd=NA, pretty=0); title("Page 14: class extra\n", cex.main=.9)
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=0, cex.main=.9, main="extra=0\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=1, cex.main=.9, main="extra=1\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=4, cex.main=.9, main="extra=4\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=5, cex.main=.9, main="extra=5\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=6, cex.main=.9, main="extra=6\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, extra=7, under=F, cex.main=.9, main="extra=7\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=8, cex.main=.9, main="extra=8\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=9, cex.main=.9, main="extra=9\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=100, cex.main=.9, main="extra=100\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=101, cex.main=.9, main="extra=101\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=104, cex.main=.9, main="extra=104\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=105, cex.main=.9, main="extra=105\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=106, cex.main=.9, main="extra=106\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=107, cex.main=.9, main="extra=107\nunder=F")
prp(a4, type=1, yesno=T, faclen=-1, under=F, extra=109, cex.main=.9, main="extra=109\n under=F")

plot(0, 0, type="n", axes=FALSE, xlab="", ylab="")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=0, cex.main=.9, main="extra=0\nunder=T")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=1, cex.main=.9, main="extra=1\nunder=T")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=4, cex.main=.9, main="extra=4\nunder=T")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=5, cex.main=.9, main="extra=5\nunder=T")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=6, cex.main=.9, main="extra=6\nunder=T")
prp(a4, type=1, yesno=T, faclen=-1, extra=7, under=T, cex.main=.9, main="extra=7\nunder=T (ignored)")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=8, cex.main=.9, main="extra=8\nunder=T")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=9, cex.main=.9, main="extra=9\nunder=T")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=100, cex.main=.9, main="extra=100\nunder=T")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=101, cex.main=.9, main="extra=101\nunder=T")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=104, cex.main=.9, main="extra=104\nunder=T")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=105, cex.main=.9, main="extra=105\nunder=T")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=106, cex.main=.9, main="extra=106\nunder=T")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=107, cex.main=.9, main="extra=107\nunder=T")
prp(a4, type=1, yesno=T, faclen=-1, under=T, extra=109, cex.main=.9, main="extra=109\n under=T")

par(old.par)

par(mfrow=c(3,3))
prp(a4, type=1,          extra=2,   main="Page 15: class extra continued\nextra=2 (classification rate)")
prp(a4, type=1, under=T, extra=3,   main="extra=3 (misclassification rate)\nunder=T")
prp(a4, type=1,          extra=102, main="extra=102 (classification rate)\n")
prp(a4, type=1, under=T, extra=103, main="extra=103 (misclassification rate)\nunder=T")

# poisson
set.seed(8)
ozone2 <- ozone1
ozone2$O3a <- round(runif(330, 1, 10))
y <- cbind(ozone2$O3, ozone2$O3a)
a5 <- rpart(y~.-O3-O3a, data=ozone2, control=list(cp=.04))
par(mfrow=c(2, 3))
plot(a5, unif=TRUE, branch=.3, main="Page 16: poisson\n"); text(a5, use.n=TRUE, all=T, digits=3, xpd=NA, cex=1.1)
prp(a5, extra=0, digits=3, type=4, trace=1, main="extra=0\ntype=4")
prp(a5, extra=1, type=4, clip.right=FALSE, under=TRUE, main="extra=1: nbr of events, nbr of obs\ntype=4", trace=1, under.cex=1)
prp(a5, extra=2, trace=1, type=0, under=T, main="extra=2: nbr of events", under.cex=1)
prp(a5, extra=102, type=4, under=TRUE, xsep="/", main="extra=102\ntype=4", trace=1, under.cex=1)

# prefix, suffix, etc.
a7 <- rpart(survived~., data=ptitanic, control=list(cp=.02))
par(mfrow=c(2, 2))
# test many parameters, and their vectorization
prp(a7, main="Page 17: many params, prefix, suffix, etc.", Margin=.03,
    extra=4, under=T, prefix="res:", suffix=" (probs)", split.suffix="\n\nabc", faclen=0, trace=3,
    nn=1,
    under.col=c(2,3), under.font=c(3,2), under.ygap=c(.2,-.2), under.cex=c(1.1, .8),
    adj=c(0,.5), split.adj=c(.5,1), yshift=c(-.5,.5),
    shadow.col=c(1,2), split.shadow.col=c("pink","blue"),
    space=c(.8,.6), yspace=c(.5,.1), border.col=c("gray", "green3", "pink"),
    lty=c(1,2),
    shadow.offset=c(.4,1,2),
    split.shadow.offset=c(.4,.4,1),
    nn.font=c(1,3), nn.adj=c(1,0), nn.col=c(1,2), nn.border.col=c(0,1,2))

my.split.labs <- function(x, labs, digits, varlen, faclen)
{
    sprintf("my.split.lab\n%s", labs)
}
prp(a7, type=4, extra=4, under=T,
    faclen=0, trace=3,
    split.fun=my.split.labs,
    split.prefix="L[", split.suffix="]L",
    right.split.prefix="R[", right.split.suffix="]R",
    round=9, leaf.round=0,
    ycompress.cex=.8) # force ycompress for testing with type=4

prp(a7, type=4, extra=1, under=F, prefix="response:",
    suffix="\n\n (probs)", split.suffix="\n\nabc", faclen=0, trace=3)
par(mfrow=c(1, 1))

data(iris)
a.iris <- rpart(Species~., data=iris)
par(mfrow=c(2, 2))
old.bg <- par(bg="darkgray")
prp(a.iris, main="Page 18: gray background",
    type=4, extra=1, under=TRUE,
    col=c("orange", "green", "wheat")[a.iris$frame$yval], under.col="red",
    border.col=c(3,4), nn.col=c(2,3),
    split.border.col=5,
    shadow.col="black",
    split.shadow.col="lightgray",
    branch.col=c("orange4", "white"),
    branch.lwd=c(3,2), branch.lty=1:3)
par(bg=old.bg)
par(mfrow=c(1, 1))

par(mfrow=c(2, 3))
a <- rpart(survived~., data=ptitanic, control=list(cp=.01))
prp(a, uniform=T, branch=.4, compress=T, extra=104, trace=2, main="Page 19: test mar and xpd args")
prp(a, uniform=T, branch=.4, compress=T, extra=104, mar=c(1,2,3,4), trace=2, main="test mar=c(1,2,3,4)")
prp(a, uniform=T, branch=.4, compress=T, extra=104, mar=c(5,2,3,4), trace=2, main="test mar=c(5,2,3,4)")
prp(a, uniform=T, branch=.4, compress=T, extra=104, xpd=T, trace=2, prefix="123456789", cex=1, main="test xpd=T, par=1")
prp(a, uniform=T, branch=.4, compress=T, extra=104, xpd=F, trace=2, prefix="123456789", cex=1, main="test xpd=F, par=1")
par(mfrow=c(1, 1))

# shadows

a <- rpart(pclass ~ ., data=ptitanic, control=rpart.control(cp=.01))
par(mfrow=c(2,3))
prp(a, type=0, faclen=0, extra=1, under=F, shadow.col="darkgray", nn=T, split.shadow.col="darkgray", main="Page 20: shadows")
prp(a, type=1, faclen=0, extra=1, under=F, shadow.col="darkgray", nn=T, main="type=1")
prp(a, type=1, faclen=0, extra=2, under=T, shadow.col="darkgray", nn=T, main="type=1")
prp(a, type=2, faclen=0, extra=3, under=F, shadow.col="darkgray", nn=T, split.shadow.col="darkgray", main="type=2")
prp(a, type=3, faclen=0, extra=4, under=T, shadow.col="darkgray", nn=T, split.shadow.col="darkgray", main="type=3")
prp(a, type=4, faclen=0, extra=101, under=T, shadow.col="darkgray", nn=T, split.shadow.col="darkgray", main="type=4")
par(mfrow=c(1,1))

# misc.

# test that do.par correctly restores eveything, also test do.par=FALSE
a <- rpart(pclass ~ ., data=ptitanic, cp=.005)
par(mfrow=c(3,3))
old.par <- par(no.readonly=TRUE)
prp(a, trace=2, main="Page 21: do.par=default") # trace=2 so can see the grid
    # set par settings that can legally change to NULL for comparison
    old.par$usr <- old.par$fig <- old.par$mfg <- old.par$xaxp <- old.par$yaxp <- NULL
    par <- par(no.readonly=TRUE)
    par$usr <- par$fig <- par$mfg <- par$xaxp <- par$yaxp <- NULL
    stopifnot(isTRUE(all.equal(old.par, par)))
prp(a, trace=2, main="do.par=FALSE", do.par=FALSE)
    par <- par(no.readonly=TRUE)
    par$usr <- par$fig <- par$mfg <- par$xaxp <- par$yaxp <- NULL
    stopifnot(isTRUE(all.equal(old.par, par)))
par(mfrow=c(1,1))

# different branch types
a <- rpart(pclass ~ ., data=ptitanic, cp=.02)
par(mfrow=c(2,3))
prp(a, branch.type=5, main="Page 22: branch.type=5\nwt")
prp(a, branch.type=1, main="branch.type=1\ndev")
prp(a, branch.type=2, main="branch.type=2\nsqrt(dev)\nuniform=FALSE", uniform=FALSE)
prp(a, branch.type=6, fallen.leaves=T, main="branch.type=6\ncomplexity\nfallen.leaves")
prp(a, branch.type=7, fallen.leaves=T, main="branch.type=7\nabs(yval)\nfallen.leaves")
prp(a, branch.type=8, main="branch.type=8\nyval - min(yval)")

par(mfrow=c(2,3))
# continuous response
a.age <- rpart(age~., data=ptitanic, cp=.04)
prp(a.age, branch.type=7, branch.col="pink", main="Page 23: branch.type=7\ncontinuous response")

# test different types with branch.type
# prp(a, type=1, branch.type=5, branch.col="slategray3", main="type=1\nbranch.type=5") # already tested
prp(a, type=2, branch.type=5, branch.col="slategray3", main="type=2\nbranch.type=5\n")
# prp(a, type=3, branch.type=5, branch.col="slategray3", main="type=1\nbranch.type=5") # not yet supported
# prp(a, type=4, branch.type=5, branch.col="slategray3", main="type=1\nbranch.type=5") # not yet supported
prp(a, type=2, branch.type=3, branch=0, branch.col="slategray3", main="type=2\nbranch.type=3\nbranch=0")
prp(a, type=2, branch.type=4, branch=1, main="type=2\nbranch.type=4\nbranch=1",
    branch.col=c("slategray","slateblue2","slateblue")[a$frame$yval])

branch.fun1 <- function(x)
{
    width <- x$frame$wt
}
root <- rpart(survived ~ ., data=ptitanic, cp=.5)
prp(a, branch.type=branch.fun1, branch.col="slategray3", main="branch.fun1")

par(mfrow=c(2,3))
prp(root, branch.type=5, main="Page 24: branch.type=5\nsingle node tree")
prp(a, branch=0, branch.type=5, branch.tweak=1.5, branch.col="slategray3",
    branch.fill=2, main="branch.type=5\nbranch args")

par(mfrow=c(4,4))
set.seed(1924)
root <- rpart(survived ~ ., data=ptitanic, cp=.5)
temp <- prp(root, main="Page 25: single node tree")
print(temp)
prp(root, type=1, main="type=1")
prp(root, type=2, extra=1, main="type=2, extra=1")
prp(root, type=3, extra=2, under=T, main="type=3, extra=4, under=T")
prp(root, type=4, extra=3, main="type=4, extra=4",
    prefix="l[", suffix="]r",
    split.prefix="L[", split.suffix="]L",
    right.split.prefix="R[", right.split.suffix="]R",
    round=9, leaf.round=0)
prp(root, branch.type=5, main="branch.type=5")
par(mfrow=c(1,1))

source("user-manual-figs.R")
par(mfrow=c(1,1))

set.seed(1924)
# use.prp <- FALSE
# source("code.in.rpart.report.with.prp.R")
use.prp <- TRUE
source("code.in.rpart.report.with.prp.R")

# clip.left.labs and clip.right.labs
par(mfrow=c(3,3))
prp(tree, type=4, clip.left.labs=F, clip.right.labs=F, main="clip.left.labs=F, clip.right.labs=F")
# prp(tree, type=4, clip.left.labs=F, clip.right.labs=T, main="clip.left.labs=F, clip.right.labs=T") #default
prp(tree, type=3, clip.left.labs=T, clip.right.labs=F, main="clip.left.labs=T, clip.right.labs=F")
prp(tree, type=3, clip.left.labs=T, clip.right.labs=T, main="clip.left.labs=T, clip.right.labs=T")

prp(tree, type=3, xflip=T, clip.left.labs=F, clip.right.labs=F, main="clip.left.labs=F, clip.right.labs=F\n                 xflip=T")
prp(tree, type=3, xflip=T, clip.left.labs=F, clip.right.labs=T, main="clip.left.labs=F, clip.right.labs=T\n                 xflip=T")
prp(tree, type=4, xflip=T, clip.left.labs=T, clip.right.labs=F, main="clip.left.labs=T, clip.right.labs=F\n                 xflip=T")
prp(tree, type=4, xflip=T, clip.left.labs=T, clip.right.labs=T, main="clip.left.labs=T, clip.right.labs=T\n                 xflip=T")

prp(tree, type=4,          clip.left.labs=c(T, F, T), clip.right.labs=c(T, F, F), main="clip.labs vectorization")
prp(tree, type=4, xflip=T, clip.left.labs=c(T, F, T), clip.right.labs=F,          main="xflip=T\nclip.labs vectorization")
par(mfrow=c(1,1))

# TODO mvpart is no longer on CRAN
#
# # mvpart, must be last because it changes plot.rpart, text.rpart, etc.
# library(mvpart)
# data(spider)
# par(mfrow=c(3,3))
# a <- mvpart(data.matrix(spider[,1:12])~twigs+water,spider, legend=FALSE, all=TRUE)
# prp(a, fallen=T, branch=1, under=T, type=0, extra=0, main="mvpart page 1\nnresp=12, extra=0")
# prp(a, fallen=T, under=T, type=1, extra=2, main="nresp=12, extra=2, under=T", under.cex=1)
# a <- mvpart(data.matrix(spider[,1:3])~twigs+water,spider, legend=FALSE, all=TRUE)
# prp(a, under=T, type=1, extra=101, main="extra=101")
# prp(a, under=T, type=2, extra=102, main="extra=102")
# prp(a, under=T, type=4, extra=3,   main="extra=3, under=F")
# prp(a, under=T, type=1, extra=4,   main="extra=4")
# prp(a, under=T, type=1, extra=105, main="extra=105")
#
# prp(a, under=F, type=4, extra=106, main="mvpart page 2\nextra=106, under=F")
# prp(a, under=T, type=4, extra=107, main="extra=107")
# prp(a, under=T, type=1, extra=8,   main="extra=8")
# prp(a, under=F, type=2, extra=109, main="extra=109, under=F")
# prp(a, under=T, type=3, extra=110, main="extra=110")
# prp(a, under=T, type=4, extra=111, main="extra=111")
# par(mfrow=c(1,1))

# # TODO this seems to not work with the new version of rpart (4.0.2)
# library(rpart.plot)
# library(rpartOrdinal)
# library(rpartScore)
# data(lowbwt)
# lowbwt <- lowbwt[1:80,]
# lowbwt$Category.s <-
#     ifelse(lowbwt$bwt <= 2500, 3,
#     ifelse(lowbwt$bwt <= 3000, 2,
#     ifelse(lowbwt$bwt <= 3500, 1,
#                                0)))
# # Gives error
# a <- rpartScore(Category.s ~ age + lwt + race + smoke +
#                 ptl + ht + ui + ftv, data = lowbwt)
# prp(a, extra=100, main="rpartScore\nextra=100", under=TRUE)

# Cannot install rpartOrdinal: package 'rpartOrdinal' is not available (for R version 3.2.0)
# library(rpartOrdinal)
# data(lowbwt)
# lowbwt$Category <- factor(
#     ifelse(lowbwt$bwt<=2500,3,
#     ifelse(lowbwt$bwt<=3000,2,
#     ifelse(lowbwt$bwt<=3500,1,
#                             0))),ordered=TRUE)
# a <- rpart(Category~age+lwt+race+smoke+ptl+ht+ui+ftv,data=lowbwt,method=ordinal)
# prp(a, main="rpartOrdinal\ntype=1, extra=0", type=1, extra=0, faclen=0)

# TODO mvpart is no longer on CRAN
#
# #--- appendix mvpart.R  ---
#
# library(mvpart)
# library(rpart.plot)
# data(spider)
# set.seed(1)
# response <- data.matrix(spider[,1:3, drop=F])
# tree1 <- mvpart(response~herbs+reft+moss+sand+twigs+water, data=spider,
#             legend=F, method="mrt", plot.add=F, xv="min")
#
# old.par <- par(par(mfrow=c(4,4)), mar = c(3, 3, 3, 1), mgp = c(1.5, .5, 0))
# prp1 <- function(tree1, extra, main, type=1, under=T, col=1, yesno=F, tweak=1,
#                  col.main="skyblue4", cex.main=1, ...)
# {
#     prp(tree1, type=type, extra=extra, main=main,
#         under=under, col=col, yesno=yesno, tweak=tweak,
#         col.main=col.main, cex.main=cex.main, ...)
# }
# prp1(tree1, extra=0, main="extra = 0\ndev", tweak=.8)
# prp1(tree1, extra=1, type=3, main="extra = 1 (type=3)\ndev,  n")
# prp1(tree1, extra=2, main="extra = 2\ndev,  frac", tweak=1.2)
# prp1(tree1, extra=3, main="extra = 3\ndev,  frac / sum(frac)")
# prp1(tree1, extra=4, main="extra = 4\nsqrt(dev)")
# prp1(tree1, extra=5, main="extra = 5\nsqrt(dev),  n")
# prp1(tree1, extra=6, main="extra = 6\nsqrt(dev),  frac", tweak=1.2)
# prp1(tree1, extra=7, main="extra = 7\nsqrt(dev),  frac / sum(frac)", tweak=1.1)
# prp1(tree1, extra=8, main="extra = 8\npredom species",   tweak=.8)
# prp1(tree1, extra=9, main="extra = 9\npredom species,  n", tweak=1)
# prp1(tree1, extra=10, main="extra = 10\npredom species,  frac", tweak=1.2)
# prp1(tree1, extra=11, main="extra = 11\npredom spec,  frac / sum(frac)", tweak=1.15)
# par(old.par)

par(mfrow=c(2,2))
source("webpage-figs.R")

# test rpart.plot version 1.5.3 (deal with situation where user has
# a variable named text in the current environment).
# Also test use of FUN argument.

# test that we got an error as expected from a try() call
expect.err <- function(object, expected.msg="")
{
    if(class(object)[1] == "try-error") {
        msg <- attr(object, "condition")$message[1]
        if(length(grep(expected.msg, msg)))
            cat("Got error as expected from ",
                deparse(substitute(object)), "\n", sep="")
        else
            stop(sprintf("Expected \"%s\"\n  but got \"%s...\"",
                         expected.msg, substr(msg, 1, 120)))
    } else
        stop("did not get expected try error")
}
cat("\ntest rpart.plot version 1.5.3\n")
par(mfrow=c(3,3))
a100 <- rpart(survived ~ ., data=ptitanic, cp=.02)
title("a100a", cex=.6)
prp(a100)
title("a100b", cex=.6)
text <- "this is not the text function"
prp(a100) # graph should be identical to the one on its left
title("a100c", cex=.6)
expect.err(try(prp(a100, FUN=function(xbad, y1, labels, ...) text(xbad, y1, labels, ...))),
           "the FUN argument to the prp function needs the following arguments")
title("a100d", cex=.6)
# user specified FUN only has to match up to the dots
prp(a100, FUN=function(x, y1, labels, ...) text(x, y1, labels, ...))

my.text <- function(x, y, labels, ...) text(x, y, labels, ...)
prp(a100, FUN=my.text)
title("a100e", cex=.6)
my.bad.text <- function(xbad, y, labels, ...) text(xbad, y, labels, ...)
expect.err(try(prp(a100, FUN=my.bad.text)),
           "the FUN argument to the prp function needs the following arguments")
title("a100f", cex=.6)
# define the function text in the global environment and use that
text <- function(x, y, labels, ...) graphics::text(x, y, paste0("my-", labels), ...)
# TODO boxes below aren't sized correctly for the user generated text
prp(a100, FUN=text)
title("a100g", cex=.6)
remove(text)

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}

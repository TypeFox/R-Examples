library(phylobase)
library(ape)

set.seed(1)

data(geospiza)

## make sure geospiza is properly formatted
if(is.character(checkval <- checkPhylo4(geospiza)))
  stop(checkval)


geospiza0 <-
  list(geospiza.tree=as(geospiza,"phylo"),geospiza.data=tipData(geospiza))
## push data back into list form as in geiger

t1 <-  try(p1 <- phylo4d(geospiza0$geospiza.tree,geospiza0$geospiza.data))
## Error in checkData(res, ...) :
##   Tip data names are a subset of tree tip labels.

p2 <- as(geospiza0$geospiza.tree,"phylo4")
plot(p2)

lab1 <- tipLabels(p2)
lab2 <- rownames(geospiza0$geospiza.data)

lab1[!lab1 %in% lab2]  ## missing data
lab2[!lab2 %in% lab1]  ## extra data (none)
p1 <- phylo4d(p2,geospiza0$geospiza.data, missing.data="warn")
p1 <- phylo4d(p2,geospiza0$geospiza.data, missing.data="OK")

plot(p1)
plot(p1,show.node.label=TRUE)
## one way to deal with it:

p1B <- prune(p1,tip="olivacea")

## or ...
p1C <- stats::na.omit(p1)

labels(p1C, "all") <- tolower(labels(p1C, "all"))

## trace("prune",browser,signature="phylo4d")
r1 <- read.tree(text="((t4:0.3210275554,(t2:0.2724586465,t3:0.2724586465):0.0485689089):0.1397952619,(t5:0.07551818331,t1:0.07551818331):0.385304634);")

## trace("phylo4d", browser, signature = "phylo")
## untrace("phylo4d", signature = "phylo")
tipdat <- data.frame(a=1:5, row.names=r1$tip.label)
q1 <- phylo4d(r1,tip.data=tipdat, node.data=data.frame(a=6:9), match.data=FALSE)
q2 <- prune(q1,1)
summary(q2)

tipdat2 <- tipdat
row.names(tipdat2)[1] <- "s1"
t1 <- try(q1 <- phylo4d(r1,tip.data=tipdat2))

plot(q2)
plot(q2,type="cladogram")
## plot(p2,type="dotchart",labels.nodes=nodeLabels(p2))
## trace("plot", browser, signature = c("phylo4d","missing"))
tipLabels(q1) <- paste("q",1:5,sep="")
nodeLabels(q1) <- paste("n",1:4,sep="")
p3 <- phylo4d(r1,tip.data=tipdat,node.data=data.frame(b=6:9), match.data=FALSE)
summary(p3)

plot(p1)

plot(subset(p1,tips.include=c("fuliginosa","fortis","magnirostris",
            "conirostris","scandens")))
## better error?
## Error in phy$edge[, 2] : incorrect number of dimensions

if(dev.cur() == 1) get(getOption("device"))()
plot(subset(p2,tips.include=c("fuliginosa","fortis","magnirostris",
            "conirostris","scandens")))

plot(p2,show.node.label=TRUE)

tree.owls <- read.tree(text="(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);")

z <- as(tree.owls,"phylo4")

example("phylo4d")
obj1 <- obj2 <- obj3 <- phylo4d(z, data.frame(wing=1:4,color=factor(c("b","w","b","b")), tail=runif(4)*10), match.data=FALSE)

obj2@data <- as.data.frame(obj2@data[,1])
obj3@data <- cbind(obj1@data,obj2@data)
obj4 <- obj1
obj4@data[2,3] <- NA
obj4@data[1,1] <- NA

nodeLabels(obj4) <- character(0)

obj5 <- obj1
tipData(obj4) <- subset(tipData(obj4),select=sapply(tipData(obj4),class)=="numeric")

treePlot(obj4)

E <- matrix(c(
    8,  9,
    9, 10,
   10,  1,
   10,  2,
    9,  3,
    9,  4,
    8, 11,
   11,  5,
   11,  6,
   11,  7,
    0,  8), ncol=2,byrow=TRUE)

P2 <- phylo4(E)

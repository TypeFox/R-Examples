## ----include=FALSE-------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(fig.width=7, fig.height=7, comment="")

## ------------------------------------------------------------------------
library(repo)

## ------------------------------------------------------------------------
rp <- repo_open(tempdir(), force=T)

## ------------------------------------------------------------------------
rp$attach("index.Rmd", "Source code for Repo vignette", c("source","Rmd"))

## ------------------------------------------------------------------------
myiris <- scale(as.matrix(iris[,1:4]))

## ------------------------------------------------------------------------
rp$put(
    obj = myiris,
    name = "myiris", 
    description = paste(
        "A normalized version of the iris dataset coming with R.",
        "Normalization is made with the scale function",
        "with default parameters."
    ),
    tags = c("dataset", "iris", "repodemo"), 
    src = "index.Rmd"
)

## ------------------------------------------------------------------------
rp$put(iris$Species, "irisLabels", "The Iris class lables.",
         c("labels", "iris", "repodemo"), "index.Rmd")

## ------------------------------------------------------------------------
irispca <- princomp(myiris)
iris2d <- irispca$scores[,c(1,2)]
plot(iris2d, main="2D visualization of the Iris dataset",
     col=rp$get("irisLabels"))

## ------------------------------------------------------------------------
fpath <- file.path(rp$root(), "iris2D.pdf")
pdf(fpath)
plot(iris2d, main="2D visualization of the Iris dataset",
     col=rp$get("irisLabels"))
invisible(dev.off())
rp$attach(fpath, "Iris 2D visualization obtained with PCA.",
            c("visualization", "iris", "repodemo"),
              "index.Rmd", to="myiris")

## ---- eval=FALSE---------------------------------------------------------
#  rp$sys("iris2D.pdf", "evince")

## ---- eval=FALSE---------------------------------------------------------
#  rp$sys("index.Rmd", "evince")

## ------------------------------------------------------------------------
plot(irispca)

## ------------------------------------------------------------------------
fpath <- file.path(rp$root(), "irisPCA.pdf")
pdf(fpath)
plot(irispca)
invisible(dev.off())
rp$attach(fpath, "Variance explained by the PCs of the Iris dataset",
            c("visualization", "iris", "repodemo"),
              "index.Rmd", to="iris2D.pdf")

## ------------------------------------------------------------------------
kiris <- kmeans(myiris, 5)$cluster
rp$put(kiris, "iris_5clu", "Kmeans clustering of the Iris data, k=5.",
         c("metadata", "iris", "kmeans", "clustering", "repodemo"),
           "index.Rmd", depends="myiris", T)

## ------------------------------------------------------------------------
plot(iris2d, main="Iris dataset kmeans clustering", col=kiris)

## ------------------------------------------------------------------------
fpath <- file.path(rp$root(), "iris2Dclu.pdf")
pdf(fpath)
plot(iris2d, main="Iris dataset kmeans clustering", col=kiris)
invisible(dev.off())
rp$attach(fpath, "Iris K-means clustering.",
	c("visualization", "iris", "clustering", "kmeans", "repodemo"),
	"index.Rmd", to="iris_5clu")

## ------------------------------------------------------------------------
res <- table(rp$get("irisLabels"), kiris)
rp$put(res, "iris_cluVsSpecies",
         paste("Contingency table of the kmeans clustering versus the",
               "original labels of the Iris dataset."),
         c("result", "iris","validation", "clustering", "repodemo", "hide"),
         "index.Rmd", c("myiris", "irisLabels", "iris_5clu"), T)

## ------------------------------------------------------------------------
rp$info()

## ------------------------------------------------------------------------
rp ## resolves to print(rp)

## ------------------------------------------------------------------------
print(rp, all=T)

## ------------------------------------------------------------------------
print(rp, tags="clustering", all=T)

## ------------------------------------------------------------------------
print(rp, tags="attachment", all=T)

## ------------------------------------------------------------------------
print(rp, tags="hide", all=T)

## ------------------------------------------------------------------------
rp$print(show="st")

## ------------------------------------------------------------------------
rp$find("clu", all=T)

## ------------------------------------------------------------------------
rp$pies()

## ---- eval=F-------------------------------------------------------------
#  rp$cpanel()

## ------------------------------------------------------------------------
rp$check()

## ------------------------------------------------------------------------
depgraph <- rp$dependencies(plot=F)
rownames(depgraph) <- colnames(depgraph) <- basename(rownames(depgraph))
library(knitr)
kable(depgraph)

## ------------------------------------------------------------------------
rp$dependencies()

## ------------------------------------------------------------------------
rp$dependencies(generated=F)

## ------------------------------------------------------------------------
x <- rp$get("myiris")

## ------------------------------------------------------------------------
rp$info("myiris")

## ------------------------------------------------------------------------
kiris2 <- kmeans(myiris, 5)$cluster
rp$put(kiris, "iris_5clu",
         "Kmeans clustering of the Iris data, k=5. Today's version!",
         c("metadata", "iris", "kmeans", "clustering", "repodemo"),
           "index.Rmd", depends="myiris", replace="addversion")

## ------------------------------------------------------------------------
rp

## ------------------------------------------------------------------------
rp$info("iris_5clu")

## ------------------------------------------------------------------------
rp$print(all=T)

## ------------------------------------------------------------------------
rp$info("iris_5clu#1")

## ----include=FALSE-------------------------------------------------------
dorun <- FALSE
result <- "This took 10 seconds to compute"
rp$stash(result)

## ------------------------------------------------------------------------
if(dorun) {
    Sys.sleep(10)
    result <- "This took 10 seconds to compute"
    rp$stash(result)
} else result <- rp$get("result")

## ------------------------------------------------------------------------
rp$info("result")

## ------------------------------------------------------------------------
expr <- expression({
	Sys.sleep(3)
	result <- "This took 3 seconds to compute"
})
	
system.time(rp$lazydo(expr)) # first run
system.time(rp$lazydo(expr)) # second run

## ------------------------------------------------------------------------
rp$put("Local content", "item1",
	"This points to big data you may want to download",
	"tag", URL="http://www.francesconapolitano.it/repo/remote")
print(rp$get("item1"))
rp$pull("item1", replace=T)
print(rp$get("item1"))

## ------------------------------------------------------------------------
h <- rp$handlers()
names(h)

## ------------------------------------------------------------------------
print(h$iris_cluVsSpecies())

## ------------------------------------------------------------------------
h$iris_cluVsSpecies("tag", "onenewtag")
h$iris_cluVsSpecies("info")

## ------------------------------------------------------------------------
h <- repo_open(rp$root())$handlers()

## ------------------------------------------------------------------------
h$repo

## ------------------------------------------------------------------------
h <- h$repo$handlers()

## ---- eval=FALSE---------------------------------------------------------
#  help(repo)

## ---- eval=FALSE---------------------------------------------------------
#  help(repo_func)

## ---- include=F----------------------------------------------------------
unlink(rp$root(), recursive=T)


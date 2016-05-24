## ----eval=TRUE,echo=FALSE,results='hide'--------------------------------------
library(knitr)
opts_chunk$set(fig.path='figs/OutbreakTools-', fig.keep='last', dev='pdf', fig.width=7, fig.height=7,
               tidy=FALSE, warning=FALSE, fig.show='asis', fig.align='center', cache=FALSE,
               out.width=".6\\textwidth")
options(width=80)

## ----results='hide'-----------------------------------------------------------
library(OutbreakTools)

## -----------------------------------------------------------------------------
getClassDef("obkData")

## -----------------------------------------------------------------------------
new("obkData")

## -----------------------------------------------------------------------------
data(ToyOutbreak)
class(ToyOutbreak)
slotNames(ToyOutbreak)
head(ToyOutbreak)
summary(ToyOutbreak)

## -----------------------------------------------------------------------------
head(ToyOutbreak@individuals)
head(ToyOutbreak@records$Fever)
ToyOutbreak@trees

## -----------------------------------------------------------------------------
class(ToyOutbreak@dna)
ToyOutbreak@dna
slotNames(ToyOutbreak@dna)
is.list(ToyOutbreak@dna@dna)
names(ToyOutbreak@dna@dna)
ToyOutbreak@dna@dna$gene1
class(ToyOutbreak@dna@dna$gene1)
class(ToyOutbreak@dna@meta)
head(ToyOutbreak@dna@meta)

## -----------------------------------------------------------------------------
cf <- c("a", "b", "a", "c", "d")
ct <- c("b", "c", "c", "d", "b")
oc.static <- new("obkContacts", cf, ct, directed=FALSE)
slotNames(oc.static)
oc.static

## ----graphstat,out.width=".7\\textwidth"--------------------------------------
plot(oc.static, main="Static contact network")

## -----------------------------------------------------------------------------
onset <- c(1, 2, 3, 4, 5)
terminus <- c(1.2, 4, 3.5, 4.1, 6)
oc.dynamic <- new("obkContacts",cf,ct, directed=FALSE,
                  start=onset, end=terminus)
slotNames(oc.dynamic)
oc.dynamic

## -----------------------------------------------------------------------------
as.data.frame(oc.dynamic)

## ----dynNet,out.width=".9\\textwidth"-----------------------------------------
par(mfrow=c(2,2))
plot(oc.dynamic@contacts,main="oc.dynamic - collapsed graph",
     displaylabels=TRUE)
plot(get.contacts(oc.dynamic, from=0, to=2),
     main="oc.dynamic - time 0--2", displaylabels=TRUE)
plot(get.contacts(oc.dynamic, from=2, to=4),
     main="oc.dynamic - time 2--4", displaylabels=TRUE)
plot(get.contacts(oc.dynamic, from=4, to=6),
     main="oc.dynamic - time 4--6", displaylabels=TRUE)

## -----------------------------------------------------------------------------
new("obkData")

## -----------------------------------------------------------------------------
data(ToyOutbreakRaw)
class(ToyOutbreakRaw)
names(ToyOutbreakRaw)

## -----------------------------------------------------------------------------
head(ToyOutbreakRaw$individuals)

## -----------------------------------------------------------------------------
lapply(ToyOutbreakRaw$records, head)

## -----------------------------------------------------------------------------
head(ToyOutbreakRaw$contacts)
head(ToyOutbreakRaw$contacts.start)
head(ToyOutbreakRaw$contacts.end)

## -----------------------------------------------------------------------------
ToyOutbreakRaw$dna

## -----------------------------------------------------------------------------
ToyOutbreakRaw$trees

## -----------------------------------------------------------------------------
attach(ToyOutbreakRaw)

x <- new ("obkData", individuals=individuals, records=records,
          contacts=contacts, contacts.start=contacts.start,
          contacts.end=contacts.end, dna=dna,
          dna.individualID=dna.info$individualID,
          dna.date=dna.info$date, sample=dna.info$sample, trees=trees)

detach(ToyOutbreakRaw)

head(x)
summary(x)

## ----eval=TRUE,echo=FALSE,results='hide'--------------------------------------
opts_chunk$set(eval=FALSE, echo=FALSE)

## ----eval=TRUE,echo=FALSE,results='hide'--------------------------------------
opts_chunk$set(eval=TRUE, echo=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  myFunction(x, y="foo")

## ----eval=FALSE---------------------------------------------------------------
#  myFunction(x, y="bar")

## -----------------------------------------------------------------------------
data(ToyOutbreak)
set.seed(1)
toKeep <- sample(get.nindividuals(ToyOutbreak),5)
toKeep
x <- subset(ToyOutbreak, individuals=toKeep)
summary(x)

## -----------------------------------------------------------------------------
get.nindividuals(x)
get.nindividuals(x, "records")
get.nindividuals(x, "dna")
get.nindividuals(x, "contacts")

## -----------------------------------------------------------------------------
get.individuals(ToyOutbreak, "contacts")

## -----------------------------------------------------------------------------
get.nlocus(x)
get.locus(x)

## -----------------------------------------------------------------------------
get.nsequences(x)
get.nsequences(x, "bylocus")
get.sequences(x)

## -----------------------------------------------------------------------------
get.trees(x)

## -----------------------------------------------------------------------------
get.dna(x)

## -----------------------------------------------------------------------------
get.dna(x, locus=2)

## -----------------------------------------------------------------------------
get.dna(x, id=c("311","222"))

## -----------------------------------------------------------------------------
get.sequences(x)
identical(get.dna(x, id=c("311","222")), get.dna(x, id=c(2,1)))

## -----------------------------------------------------------------------------
get.ncontacts(ToyOutbreak)
get.individuals(ToyOutbreak@contacts)
get.individuals(x)
get.ncontacts(x)

## ----getData1-----------------------------------------------------------------
get.data(x,"temperature")
get.data(x,"temperature", showSource=TRUE)

## -----------------------------------------------------------------------------
get.data(x, "Sex")

## -----------------------------------------------------------------------------
get.data(x, c("Sex","Age","infector"))

## -----------------------------------------------------------------------------
get.data(x, c("Sex","Age","infector"), showSource=TRUE)

## -----------------------------------------------------------------------------
get.data(x, "date")

## -----------------------------------------------------------------------------
get.data(x, "date", showSource=TRUE)

## -----------------------------------------------------------------------------
get.data(x, "date", where="records", showSource=TRUE)

## ----sugarman, warning=TRUE---------------------------------------------------
get.data(x, "sugarman")

## ----out.width==".8\\textwidth"-----------------------------------------------
x <- subset(ToyOutbreak, individuals=1:10)
get.ncontacts(x)
plot(x@contacts, main="Contacts in x", label.cex=1.25, vertex.cex=2)

## -----------------------------------------------------------------------------
as.matrix(x@contacts)
as.matrix(x@contacts, "edgelist")

## -----------------------------------------------------------------------------
as.data.frame(x@contacts)

## ----eval=FALSE, tidy=FALSE---------------------------------------------------
#  subset(x, individuals=NULL, locus=NULL, sequences=NULL,
#         date.from=NULL, date.to=NULL, date.format=NULL, ...)

## ----subset1------------------------------------------------------------------
data(ToyOutbreak)
x1 <- subset(ToyOutbreak, individuals=1:10)
x2 <- subset(ToyOutbreak, get.individuals(ToyOutbreak)[1:10])
identical(x1,x2)

## -----------------------------------------------------------------------------
data(FluH1N1pdm2009)
attach(FluH1N1pdm2009)

x <- new("obkData", individuals = individuals, dna = FluH1N1pdm2009$dna,
      dna.individualID = samples$individualID, dna.date = samples$date,
      trees = FluH1N1pdm2009$trees)

detach(FluH1N1pdm2009)

range(get.data(x, "date"))

## ----subsetdate---------------------------------------------------------------
min.date <- min(get.dates(x))
min.date
min.date+31
x1 <- subset(x, date.to=min.date+31)
summary(x)
summary(x1)

## ----lastsubset---------------------------------------------------------------
temp <- get.data(x, "location", showSource=TRUE)
head(temp)
toKeep <- temp$individualID[temp$location=="Europe"]
x.summerEur <- subset(x, date.from="01/06/2009", date.to="31/08/2009",
                      indiv=toKeep)
summary(x.summerEur)
head(x.summerEur)

## -----------------------------------------------------------------------------
x.summerEur@trees <- NULL
get.nsequences(x.summerEur)

## ----makephylo, fig.keep="all"------------------------------------------------
x2 <- make.phylo(x.summerEur)
summary(x2)

## ----out.width=".75\\textwidth"-----------------------------------------------
library(ape)
plot(get.trees(x2)[[1]])
axisPhylo()

## -----------------------------------------------------------------------------
plot(x2, "phylo")

## ----tree1,fig.keep="last",out.width=".75\\textwidth"-------------------------
x3 <- make.phylo(x.summerEur, locus=1, ask=FALSE, model="K80")
plot(get.trees(x3)[[1]])
axisPhylo()

## ----fig.width=10, out.width="\\textwidth"------------------------------------
set.seed(1)
x <- simuEpi(N=50, D=20, beta=0.01,plot=TRUE,makePhylo=TRUE)
summary(x)
x$dynamics
summary(x$x)

## -----------------------------------------------------------------------------
plot(x$x, "contacts", main="Transmission tree")

## -----------------------------------------------------------------------------
plot(x$x, "phylo")

## -----------------------------------------------------------------------------
data(HorseFlu)
summary(HorseFlu)

## ----plottime1----------------------------------------------------------------
plot(HorseFlu,'timeline')

## -----------------------------------------------------------------------------
args(plotIndividualTimeline)

## -----------------------------------------------------------------------------
plot(HorseFlu,'timeline', what="Vac")

## -----------------------------------------------------------------------------
plotIndividualTimeline(HorseFlu, what="dna", colorBy="yardID", orderBy="yardID",plotNames=TRUE)

## ----plotfirst20--------------------------------------------------------------
plot(HorseFlu,selection=1:20, colorBy="yardID", orderBy="yardID", size=5)

## -----------------------------------------------------------------------------
data(ToyOutbreak)
head(ToyOutbreak@individuals)

## ----plotgeo1,dev='png',results='hide',message=FALSE--------------------------
plot(ToyOutbreak,'geo', location=c('lon','lat'), zoom=14)

## ----dev='png',results='hide',message=FALSE-----------------------------------
plot(ToyOutbreak,'geo', location=c('lon','lat'), zoom=15,
     colorBy='Sex', center='11')

## ----plotmst1, cache=TRUE, fig.keep="last"------------------------------------
data(HorseFlu)
plot(HorseFlu,'mst')

## -----------------------------------------------------------------------------
plot(HorseFlu,'mst',individualID=42)

## -----------------------------------------------------------------------------
data(FluH1N1pdm2009)
attach(FluH1N1pdm2009)

x <- new("obkData", individuals = individuals, dna = FluH1N1pdm2009$dna,
      dna.individualID = samples$individualID, dna.date = samples$date,
      trees = FluH1N1pdm2009$trees)

detach(FluH1N1pdm2009)

summary(x)

## -----------------------------------------------------------------------------
get.trees(x)
tre <- get.trees(x)[[1]]
tre

## ----pdh1n1tree1, out.width="0.8\\textwidth"----------------------------------
plot(get.trees(x)[[1]], show.tip=FALSE)

## ----pdh1n1,dev='png'---------------------------------------------------------
plot(x, colorBy="location", orderBy="location")

## ----fig.keep="last", out.width="0.8\\textwidth"------------------------------
plotggphy(x)

## ----pdh1n1tree2, out.width="0.8\\textwidth"----------------------------------
p <- plotggphy(x, ladderize = TRUE,  branch.unit = "year")

## -----------------------------------------------------------------------------
head(x@individuals)

## ----pdh1n1tree3, out.width="\\textwidth"-------------------------------------
p <- plotggphy(x, ladderize = TRUE, branch.unit = "year",
               tip.color = "location", tip.size = 3, tip.alpha = 0.75)


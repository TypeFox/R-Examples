## ----setup, include=FALSE------------------------------------------------
library(knitr)
library(RevEcoR)
opts_chunk$set(fig.width=8, fig.height=5)
set.seed(60823316) 

## ----eval=FALSE----------------------------------------------------------
#  install.packages("RevEcoR")

## ----eval=FALSE----------------------------------------------------------
#  if (!require(devtools)
#    install.packages("devtools")
#  devtools::install_github("yiluheihei/RevEcoR")

## ----eval=TRUE-----------------------------------------------------------
library(RevEcoR) 

## ----eval=FALSE----------------------------------------------------------
#  ## download sample metabolic data from remote KEGG database
#  buc <- getOrgMetabolicData("buc")
#  data(kegg_buc)
#  head(buc)

## ----eval=TRUE, htmlcap="Figure 1 Reconstruction metabolic network of *Buchnera aphidicola APS*", fig.lp="Figure 1", fig.width=8, fig.height=8----
## species in KEGG 
buc.net <- reconstructGsMN(kegg_buc, RefData = NULL) 
igraph::print.igraph(buc.net) 
igraph::plot.igraph(buc.net, vertex.label=NA, vertex.size=5, edge.arrow.size=0.1)
  
## ko annotation profile species detected in a human microbiome in IMG (not in KEGG) 
annodir <- system.file("extdata/koanno.tab",package = "RevEcoR") 
metabolic.data <- read.delim(annodir,stringsAsFactors=FALSE) 
##load the reference metabolic data 
data(RefDbcache) 
g2 <- reconstructGsMN(metabolic.data, RefData = RefDbcache) 

## ----eval=TRUE, htmlcap="Figure 2The node colored with red represents the species' seed set",fig.lp="Figure 2", fig.width=8, fig.height=8----
## seed set prediction
seed.set <- getSeedSets(buc.net, 0.2) 
show(seed.set) 
head(seed.set@seeds)
## The node colored with red represents the species' seed set
nodes  <- igraph::V(buc.net)$name
seeds  <- unlist(seed.set@seeds)
seed.index  <- match(seeds,nodes)
node.color <- rep("SkyBlue2",length(nodes))
node.color[seed.index]  <- "red"
igraph::plot.igraph(buc.net, 
          vertex.label=NA, vertex.size=5, edge.arrow.size=0.1,
          vertex.color = node.color)

## ------------------------------------------------------------------------
# ptr metabolic network 
data(kegg_ptr) 
##ptr.net <- reconstructGsMN(getOrgMetabolicData("ptr")) 
ptr.net <- reconstructGsMN(kegg_ptr) 
# cooperation analysis between buc and ptr 
cooperation.index <- calculateCooperationIndex(buc.net,ptr.net) 
cooperation.index 

## ---- eval = FALSE, echo=TRUE--------------------------------------------
#  ##metabolic network reconstruction of these seven species
#  net  <- lapply(anno.species, reconstructGsMN)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  ## caculate interactions among vious species
#  interactions  <- calculateCooperationIndex(net, p = TRUE)
#  ## competition index
#  $competition.index
#            Aa        Ao        Fn        Pg        Sg        So        Va
#  Aa 1.0000000 0.4736842 0.3157895 0.2280702 0.4210526 0.4385965 0.2456140
#  Ao 0.4736842 1.0000000 0.3684211 0.3333333 0.4736842 0.4736842 0.2456140
#  Fn 0.5000000 0.5833333 1.0000000 0.4166667 0.5833333 0.5555556 0.4166667
#  Pg 0.4193548 0.6129032 0.4838710 1.0000000 0.6129032 0.5161290 0.3870968
#  Sg 0.5454545 0.6136364 0.4772727 0.4318182 1.0000000 0.9090909 0.3863636
#  So 0.5813953 0.6046512 0.4651163 0.3720930 0.9302326 1.0000000 0.3953488
#  Va 0.4827586 0.4827586 0.5172414 0.4137931 0.5862069 0.5862069 1.0000000
#  ## p value of competition index
#  $competition.index.p
#  Aa    Ao    Fn    Pg    Sg    So    Va
#  Aa 0.000 0.001 0.001 0.001 0.001 0.001 0.001
#  Ao 0.001 0.000 0.001 0.001 0.001 0.001 0.001
#  Fn 0.001 0.001 0.000 0.001 0.001 0.001 0.001
#  Pg 0.001 0.001 0.001 0.000 0.001 0.001 0.001
#  Sg 0.001 0.001 0.001 0.001 0.000 0.001 0.001
#  So 0.001 0.001 0.001 0.001 0.001 0.000 0.001
#  Va 0.001 0.001 0.001 0.001 0.001 0.001 0.000
#  ## complementarity index
#  $complementarity.index
#   Aa        Ao        Fn         Pg        Sg         So        Va
#  Aa 0.0000000 0.1052632 0.1228070 0.07017544 0.0877193 0.08771930 0.1228070
#  Ao 0.1403509 0.0000000 0.1403509 0.07017544 0.1228070 0.12280702 0.1403509
#  Fn 0.1944444 0.1666667 0.0000000 0.16666667 0.1111111 0.11111111 0.1388889
#  Pg 0.2258065 0.2258065 0.1612903 0.00000000 0.1612903 0.19354839 0.2258065
#  Sg 0.2272727 0.1818182 0.1590909 0.09090909 0.0000000 0.04545455 0.1590909
#  So 0.1860465 0.1395349 0.1860465 0.09302326 0.0000000 0.00000000 0.1395349
#  Va 0.2068966 0.1724138 0.1379310 0.17241379 0.1379310 0.13793103 0.0000000
#  ## p value of complementarity index
#  $complementarity.index.p
#   Aa    Ao    Fn    Pg    Sg    So    Va
#  Aa 0.000 0.001 0.001 0.001 0.001 0.001 0.001
#  Ao 0.001 0.000 0.001 0.001 0.001 0.001 0.001
#  Fn 0.001 0.001 0.000 0.001 0.001 0.001 0.001
#  Pg 0.001 0.001 0.001 0.000 0.001 0.001 0.001
#  Sg 0.001 0.001 0.001 0.001 0.000 0.001 0.001
#  So 0.001 0.001 0.001 0.001 0.001 0.000 0.001
#  Va 0.001 0.001 0.001 0.001 0.001 0.001 0.000

## ------------------------------------------------------------------------
data(gut_microbiome) 
## summary(gut_microbiome) 

## ---- eval = FALSE, echo = TRUE------------------------------------------
#  gut.nets <- lapply(gut_microbiome,reconstructGsMN)
#  seed.sets <- lapply(gut.nets,getSeedSets)
#  ## Since calculation is on large scale, species interactions prediction may take several hours
#  gut.interactions <- calculateCooperationIndex(gut.nets)
#  competition.index <- gut.interactions$competition.index
#  complementarity.index <- gut.interactions$complementarity.index

## ---- eval = TRUE, echo = TRUE-------------------------------------------
occurrence.score <- read.delim(system.file("extdata/occurrence.tab",
  package = "RevEcoR"),stringsAsFactors = FALSE, quote = "")

## ---- eval=FALSE,echo=TRUE-----------------------------------------------
#  competition.index <- (competition.index + t(competition.index))/2
#  complementarity.index <- (complementarity.index + t(complementarity.index))/2

## ---- eval=FALSE,echo=TRUE-----------------------------------------------
#  ## upper triangles, which is used to calculate the correlation
#  competition.upper <- competition.index[upper.tri(competition.index)]
#  occurrence.upper <- occurrence.score[upper.tri(occurrence.score)]
#  complementarity.upper <- complementarity.index[upper.tri(complementarity.index)]
#  
#  ## calculate the spearman correlation betwwen co-occurrence scores and two
#  ## interactions indices
#  competition.cor <- cor(competition.upper,occurrence.upper,method="spearman")
#  complementarity.cor <- cor(complementarity.upper,occurrence.upper,method="spearman")
#  
#  ## permutation-based mantel test. Random permutation the co-occurance score
#  ## 10000 times, P value is the fraction of correlations as high as or higher
#  ## than the original
#  if (require(magrittr)){
#    null.stat <- replicate(10000,
#      sample(1:116) %>% occurrence.score[.,.] %>%
#        .[upper.tri(.)]
#    )
#    competition.null <- cor(competition.upper,null.stat)
#    complementarity.null <- cor(complementarity.upper,null.stat)
#    length(which(competition.null >= competition.cor)) ## 0 p.competition < 0.00001
#    length(which(complementarity.null <= complementarity.cor)) ## 0 p.complementarity< 0.00001
#  }

## ---- eval=TRUE----------------------------------------------------------
sessionInfo() 


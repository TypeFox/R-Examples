## ----egg_data, cache=TRUE------------------------------------------------
library(rotl)
library(fulltext)

## if (suppressMessages(require(gdata))) {
##     doi <- "10.1111/jeb.12282"
##     xl_file <- ft_get_si(doi, 1, save.name="egg.xls")
##     egg_data <- read.xls(xls = xl_file)
## } else {
    egg_data <- read.csv(system.file("extdata", "egg.csv", package = "rotl"))
#}
head(egg_data)

## ----eggs_in_a_funnel, fig.width=6, fig.height=3-------------------------
plot(1/sqrt(egg_data$VZr), egg_data$Zr, pch=16,
     ylab="Effect size (Zr)",
     xlab="Precision (1/SE)",
     main="Effect sizes for sex bias in egg size among 51 brid species" )

## ---- clean_eggs---------------------------------------------------------
egg_data <- as.data.frame(egg_data)
egg_data$animal <- tolower(egg_data$animal)

## ---- birds, cache=TRUE--------------------------------------------------
taxa <- tnrs_match_names(unique(egg_data$animal), context="Animals")
head(taxa)

## ----bird_map------------------------------------------------------------
taxon_map <- structure(taxa$search_string, names=taxa$unique_name)

## ----odd_duck------------------------------------------------------------
taxon_map["Anser caerulescens"]

## ----birds_in_a_tree, fig.width=5, fig.height=5, fig.align='center'------
tr <- tol_induced_subtree(taxa$ott_id)
plot(tr, show.tip.label=FALSE)

## ----tip_lab-------------------------------------------------------------
tr$tip.label[1:4]

## ----clean_tips----------------------------------------------------------
otl_tips <- strip_ott_ids(tr$tip.label, remove_underscores=TRUE)
tr$tip.label <- taxon_map[ otl_tips ]

## ----remove_nodes--------------------------------------------------------
tr$node.label <- NULL

## ----model---------------------------------------------------------------
library(MCMCglmm, quiet=TRUE)
set.seed(123)

pr<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=0.002))
)

model <- MCMCglmm(Zr~1,random=~animal,
                       pedigree=tr,
                       mev=egg_data$VZr,
                       prior=pr,
                       data=egg_data,
                       verbose=FALSE)

## ----PhyH----------------------------------------------------------------
var_comps <- colMeans(model$VCV )
var_comps["animal"] / sum(var_comps)


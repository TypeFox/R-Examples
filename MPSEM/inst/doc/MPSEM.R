## ----load_package--------------------------------------------------------
library(MPSEM)

## ----load_data-----------------------------------------------------------
data(perissodactyla,package="caper")

## ----plot_phylogeny,echo=FALSE,fig.height=4.5----------------------------
plot(perissodactyla.tree)

## ----data_table,results="asis",echo=FALSE--------------------------------
library(xtable)
xtable(perissodactyla.data[,c(1L,2L,4L)])

## ----droping_species-----------------------------------------------------
spmatch <- match(perissodactyla.tree$tip.label,
                 perissodactyla.data[,1L])
perissodactyla.tree <- drop.tip(perissodactyla.tree,
                  perissodactyla.tree$tip.label[is.na(spmatch)])

## ----check_order---------------------------------------------------------
cbind(perissodactyla.tree$tip.label,perissodactyla.data[,1L])

## ----re-order_species----------------------------------------------------
spmatch <- match(perissodactyla.tree$tip.label,
                 perissodactyla.data[,1L])
perissodactyla.data <- perissodactyla.data[spmatch,]
all(perissodactyla.tree$tip.label==perissodactyla.data[,1L])

## ----change_rownames-----------------------------------------------------
rownames(perissodactyla.data) <- perissodactyla.data[,1L]
perissodactyla.data <- perissodactyla.data[,-1L]

## ----re-arranged_data,results="asis",echo=FALSE--------------------------
xtable(perissodactyla.data[,c(1L,3L)])

## ----training_testing_datasets-------------------------------------------
perissodactyla.train <- perissodactyla.data[-1L,,drop=FALSE]
perissodactyla.test <- perissodactyla.data[1L,,drop=FALSE]
perissodactyla.tree.train <- drop.tip(perissodactyla.tree,
                             tip="Dicerorhinus sumatrensis")

## ----display_weighting,echo=FALSE,fig.height=5---------------------------
par(mar=c(4.5,4.5,1,7)+0.1)
d <- seq(0,2,length.out=1000)
a <- c(0,0.33,0.67,1,0.25,0.75,0)
psi <- c(1,1,1,1,0.65,0.65,0.4)
cc <- c(1,1,1,1,1,1,1)
ll <- c(1,2,2,2,3,3,3)
trial <- cbind(a,psi)
colnames(trial) <- c("a","psi")
ntrials <- nrow(trial)
nd <- length(d)
w <- matrix(NA,ntrials,nd,dimnames=list(paste("a=",trial[,"a"],", psi=",trial[,"psi"],sep=""),
                                        paste("d=",round(d,4),sep="")))
for(i in 1:ntrials)
  w[i,] <- MPSEM::PEMweights(d,trial[i,"a"],trial[i,"psi"])
plot(NA,xlim=c(0,2),ylim=c(0,1.6),ylab=expression(paste(italic(w[list(italic(a),psi)]),~(phi))),
  xlab=expression(paste("Distance (",italic(phi),")",sep="")),axes=FALSE)
axis(1,at=seq(0,2,0.5),label=seq(0,2,0.5))
axis(2,las=1)
text(expression(paste(~~~a~~~~~~~psi)),x=2.2,y=1.57,xpd=TRUE,adj=0)
for(i in 1:ntrials) {
  lines(x=d,y=w[i,],col=cc[i],lty=ll[i])
  text(paste(sprintf("%.2f",trial[i,1]),sprintf("%.2f",trial[i,2]),sep="  "),
       x=rep(2.2,1),y=w[i,1000],xpd=TRUE,adj=0)
}
rm(d,a,psi,cc,ll,trial,ntrials,nd,w,i)

## ----convert_to_graph----------------------------------------------------
perissodactyla.pgraph <- 
               Phylo2DirectedGraph(perissodactyla.tree.train)

## ----graph_storage,echo=FALSE,size="tiny"--------------------------------
str(perissodactyla.pgraph)

## ----tree_labelled,fig.height=5------------------------------------------
tree <- perissodactyla.tree.train
tree$node.label <- paste("N",1L:tree$Nnode)
plot(tree,show.node.label=TRUE)
edgelabels(1L:nrow(tree$edge),
           edge=1L:nrow(tree$edge),bg="white",cex=0.75)
rm(tree)

## ----set_param-----------------------------------------------------------
steepness <- rep(0,attr(perissodactyla.pgraph,"ev")[1L])
evol_rate <- rep(1,attr(perissodactyla.pgraph,"ev")[1L])
steepness[15L:21] <- 0.25
evol_rate[15L:21] <- 2
steepness[9L:13] <- 0.8
evol_rate[9L:13] <- 0.5

## ----calculate_PEM-------------------------------------------------------
perissodactyla.PEM <- PEM.build(perissodactyla.pgraph,
                                d="distance",sp="species",
                                a=steepness,psi=evol_rate)

## ----Eigenvector_example,fig.height=3.5,fig.width=4.5--------------------
layout(matrix(c(1,1,1,2,2,3,3),1L,7L))
par(mar=c(5.1,2.1,4.1,2.1))
plot(perissodactyla.tree.train,x.lim=60,cex=0.75)
plot(y = 1L:nrow(perissodactyla.train), ylab="", xlab = "Loading",
     x = as.data.frame(perissodactyla.PEM)[,1L], xlim=0.5*c(-1,1),
     axes=FALSE, main = expression(bold(v)[1]))
axis(1) ; abline(v=0)
plot(y = 1L:nrow(perissodactyla.train), ylab="", xlab = "Loading",
     x = as.data.frame(perissodactyla.PEM)[,5L], xlim=0.5*c(-1,1),
     axes=FALSE, main = expression(bold(v)[5]))
axis(1) ; abline(v=0)

## ----PEM_opt1------------------------------------------------------------
perissodactyla.PEM_opt1 <- PEM.fitSimple(
                     y = perissodactyla.train[,"log.neonatal.wt"],
                     x = NULL,
                     w = perissodactyla.pgraph,
                     d = "distance", sp="species",
                     lower = 0, upper = 1)

## ----PEM_opt2------------------------------------------------------------
perissodactyla.PEM_opt2 <- PEM.fitSimple(
                     y = perissodactyla.train[,"log.neonatal.wt"],
                     x = perissodactyla.train[,"log.female.wt"],
                     w = perissodactyla.pgraph,
                     d = "distance", sp="species",
                     lower = 0, upper = 1)

## ----build_PEM_models----------------------------------------------------
lm1 <- lmforwardsequentialAICc(
                     y = perissodactyla.train[,"log.neonatal.wt"],
                     object = perissodactyla.PEM_opt1)
summary(lm1)
lm2 <- lmforwardsequentialAICc(
            y = perissodactyla.train[,"log.neonatal.wt"],
            x = perissodactyla.train[,"log.female.wt",drop=FALSE],
            object = perissodactyla.PEM_opt2)
summary(lm2) 

## ----make_prediction-----------------------------------------------------
perissodactyla.loc <- getGraphLocations(perissodactyla.tree,
                              targets="Dicerorhinus sumatrensis")
pred <- predict(object=perissodactyla.PEM_opt2,
                targets=perissodactyla.loc,
                lmobject=lm2,
                newdata=perissodactyla.test,
                "prediction",0.95)

## ----cross-validation----------------------------------------------------
perissodactyla.data <- data.frame(perissodactyla.data,
                         predictions = NA, lower = NA, upper = NA)
jackinfo <- list()
for(i in 1L:nrow(perissodactyla.data)) {
  jackinfo[[i]] <- list()
  jackinfo[[i]][["loc"]] <- getGraphLocations(perissodactyla.tree,
                       targets = rownames(perissodactyla.data)[i])

  jackinfo[[i]][["PEM"]] <- PEM.fitSimple(
                    y = perissodactyla.data[-i,"log.neonatal.wt"],
                    x = perissodactyla.data[-i,"log.female.wt"],
                    w = jackinfo[[i]][["loc"]]$x)
  jackinfo[[i]][["lm"]] <- lmforwardsequentialAICc(
           y = perissodactyla.data[-i,"log.neonatal.wt"],
           x = perissodactyla.data[-i,"log.female.wt",drop=FALSE],
           object = jackinfo[[i]][["PEM"]])
  predictions <- predict(object = jackinfo[[i]][["PEM"]],
      targets = jackinfo[[i]][["loc"]],
      lmobject = jackinfo[[i]][["lm"]],
      newdata = perissodactyla.data[i,"log.female.wt",drop=FALSE],
      "prediction",0.95)
  perissodactyla.data[i, c("predictions", "lower", "upper")] <-
                      unlist(predictions)
} ; rm(i, predictions)

## ----plot_pred_obs,echo=FALSE,fig.height=4,fig.width=4-------------------
par(mar=c(5,5,2,2)+0.1)
rng <- range(perissodactyla.data[,"log.neonatal.wt"], perissodactyla.data[,c("predictions","lower","upper")])
plot(NA, xlim = rng, ylim = rng, xlab = "Predicted", ylab = "observed", asp = 1, las = 1)
points(x = perissodactyla.data[,"predictions"], y = perissodactyla.data[,"log.neonatal.wt"])
abline(0,1)
arrows(x0 = perissodactyla.data[,"lower"],x1 = perissodactyla.data[,"upper"],
       y0 = perissodactyla.data[,"log.neonatal.wt"],
       y1 = perissodactyla.data[,"log.neonatal.wt"],
       length = 0.05,angle = 90,code = 3)

## ----influence_matrix,echo=TRUE,eval=FALSE-------------------------------
#  res <- PEMInfluence(perissodactyla.pgraph)

## ----PEM_updater,echo=TRUE,eval=FALSE------------------------------------
#  res <- PEM.updater(object = perissodactyla.PEM, a = 0, psi = 1)

## ----forcedSimple,echo=TRUE,eval=FALSE-----------------------------------
#  res <- PEM.forcedSimple(
#                      y = perissodactyla.train[,"log.neonatal.wt"],
#                      x = perissodactyla.train[,"log.female.wt"],
#                      w = perissodactyla.pgraph,
#                      a = steepness, psi = evol_rate)

## ----get_scores,echo=TRUE,eval=FALSE-------------------------------------
#  scores <- Locations2PEMscores(object = perissodactyla.PEM_opt2,
#                                gsc = perissodactyla.loc)


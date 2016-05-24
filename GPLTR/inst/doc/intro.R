### R code from vignette source 'intro.Rnw'

###################################################
### code chunk number 1: intro.Rnw:41-46
###################################################
options(continue = "  ", width = 60)
options(SweaveHooks=list(fig=function() par(mar = c(4.1, 4.1, 0.8, 1.1))))
pdf.options(pointsize = 10)
par(xpd = NA)  #stop clipping
library(GPLTR)


###################################################
### code chunk number 2: cart
###################################################
getOption("SweaveHooks")[["fig"]]()
data(burn)
head(burn, n = 10)

rpart.burn <- rpart(D2 ~ Z1  + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8 + Z9 
                  + Z10 + Z11, data = burn, method = "class")

print(rpart.burn)

#par(mar = rep(0.1, 4))

plot(rpart.burn)

text(rpart.burn, xpd = TRUE, cex = .6, use.n = TRUE)


###################################################
### code chunk number 3: gpltr
###################################################
getOption("SweaveHooks")[["fig"]]()
## use example(GPLTR) to have a brief overview of the contain of the package.

## fit the PLTR model after adjusting on gender (Z2) using the proposed method
## ?GPLTR or ?pltr.glm to access the help section 

## setting the parameters

args.rpart <- list(minbucket = 10, maxdepth = 4, cp = 0, 
                       maxcompete = 0, maxsurrogate = 0)
family <- "binomial"
X.names = "Z2"
Y.name = "D2"
G.names = c('Z1','Z3','Z4','Z5','Z6','Z7','Z8','Z9','Z10','Z11')

## Build the maximal tree with an adjustment  on gender (Z2)

pltr.burn <- pltr.glm(burn, Y.name, X.names, G.names, args.rpart = 
                      args.rpart, family = family,iterMax =8, iterMin = 6,
                        verbose = TRUE)

## Prunned back the maximal tree using either the BIC or the AIC criterion

pltr.burn_prun <- best.tree.BIC.AIC(xtree = pltr.burn$tree, burn, Y.name, 
                        X.names, family = family, verbose = FALSE)

## Summary of the selected tree by a BIC criterion

summary(pltr.burn_prun$tree$BIC)

## Summary of the final selected pltr model

summary(pltr.burn_prun$fit_glm$BIC)

## Plot the maximal tree and the BIC prunned tree

par(mfrow = c(2,1))

plot(pltr.burn$tree, main = '', margin = 0.05)

text(pltr.burn$tree, xpd = TRUE, cex = .4, col = 'blue')

plot(pltr.burn_prun$tree$BIC, branch = .5, main = '', margin = 0.05)

text(pltr.burn_prun$tree$BIC, xpd = TRUE, cex = .4, col = 'blue' )


###################################################
### code chunk number 4: intro.Rnw:264-279
###################################################
set.seed(150)
pltr.burn_CV <- best.tree.CV(pltr.burn$tree, burn, Y.name, X.names, 
                  G.names, family = family, args.rpart = args.rpart,
                  epsi = 0.001, iterMax = 15, iterMin = 8, ncv = 10,
                  verbose = FALSE) 

pltr.burn_CV$CV_ERROR

Bic_size <- sum(pltr.burn_prun$tree$BIC$frame$var == '<leaf>')

## Bic_size <- tree_select$best_index[[1]]

CV_ERROR_BIC <- pltr.burn_CV$CV_ERROR[[2]][Bic_size]

CV_ERROR_BIC


###################################################
### code chunk number 5: intro.Rnw:287-300 (eval = FALSE)
###################################################
## ## Use only one worker on a window plateform.
## 
## args.parallel = list(numWorkers = 10, type = "PSOCK")
## 
## index = Bic_size
## 
## p_value <- p.val.tree(xtree = fit_pltr$tree, data_pltr, Y.name, X.names,
##             G.names, B = 1000, args.rpart = args.rpart, epsi = 1e-3, 
##             iterMax = 15, iterMin = 8, family = family, LB = FALSE, 
##             args.parallel = args.parallel, index = index, verbose =
##             FALSE)
## 
## p_value$P.value


###################################################
### code chunk number 6: treesbag
###################################################
getOption("SweaveHooks")[["fig"]]()
 ##  ?bagging.pltr
set.seed(250)

Bag.burn <-  bagging.pltr(burn, Y.name, X.names, G.names, family, 
              args.rpart,epsi = 0.01, iterMax = 4, iterMin = 3, 
              Bag = 20, verbose = FALSE, doprune = FALSE)

## The thresshold values used

Bag.burn$CUT

##The set of PLTR models in the bagging procedure

PLTR_BAG.burn <- Bag.burn$Glm_BAG

##The set of trees in the bagging procedure

TREE_BAG.burn <- Bag.burn$Tree_BAG

## Look for the variability of trees in the bagging sequence

par(mfrow = c(3,2))

plot(TREE_BAG.burn[[1]], branch = .5, main = '', margin = 0.05)
text(TREE_BAG.burn[[1]], xpd = TRUE, cex = .6 )

plot(TREE_BAG.burn[[2]], branch = .5, main = '', margin = 0.05)
text(TREE_BAG.burn[[2]], xpd = TRUE, cex = .6 )

plot(TREE_BAG.burn[[3]], branch = .5, main = '', margin = 0.05)
text(TREE_BAG.burn[[3]], xpd = TRUE, cex = .6 )

plot(TREE_BAG.burn[[4]], branch = .5, main = '', margin = 0.05)
text(TREE_BAG.burn[[4]], xpd = TRUE, cex = .6 )

plot(TREE_BAG.burn[[5]], branch = .5, main = '', margin = 0.05)
text(TREE_BAG.burn[[5]], xpd = TRUE, cex = .6 )

plot(TREE_BAG.burn[[6]], branch = .5, main = '', margin = 0.05)
text(TREE_BAG.burn[[6]], xpd = TRUE, cex = .6 )


###################################################
### code chunk number 7: intro.Rnw:383-398
###################################################
## Use the bagging procedure to predict new features
# ?predict_bagg.pltr

Pred_Bag.burn <- predict_bagg.pltr(Bag.burn, Y.name, newdata = burn, 
                type = "response", thresshold = seq(0, 1, by = 0.1))

## The confusion matrix for each thresshold value using the majority vote

Pred_Bag.burn$CONF1

## The prediction error for each thresshold value

Pred_err.burn <- Pred_Bag.burn$PRED_ERROR1

Pred_err.burn


###################################################
### code chunk number 8: varimp
###################################################
getOption("SweaveHooks")[["fig"]]()
Var_Imp_BAG.burn <- VIMPBAG(Bag.burn, burn, Y.name)

## Importance score using the permutaion method for each thresshold value

Var_Imp_BAG.burn$PIS

par(mfrow=c(1,3))

barplot(Var_Imp_BAG.burn$PIS$CUT5, main = 'PIS', horiz = TRUE, las = 1,
        cex.names = .8, col = 'lightblue')

barplot(Var_Imp_BAG.burn$DIS, main = 'DIS', horiz = TRUE, las = 1,
        cex.names = .8, col = 'grey') 

barplot(Var_Imp_BAG.burn$DDIS, main = 'DDIS', horiz = TRUE, las = 1,
        cex.names = .8, col = 'purple')


###################################################
### code chunk number 9: rocoob
###################################################
getOption("SweaveHooks")[["fig"]]()
auc_BAG_oob <- bag.aucoob(Bag.burn, burn, Y.name)

## AUC of the predictor on OOB samples

auc_BAG_oob$AUCOOB

## Plot the ROC curve of the predictor based on OOB samples

par(mfrow=c(1, 1))

plot(auc_BAG_oob$FPR, auc_BAG_oob$TPR, type = 'b', lty = 3, col = 'blue', 
     xlab = 'false positive rate', ylab = 'true positive rate')

legend(0.7, 0.3, sprintf('%3.3f', auc_BAG_oob$AUCOOB), lty = c(1, 1), 
       lwd = c(2.5, 2.5), col = 'blue', title = 'AUC')



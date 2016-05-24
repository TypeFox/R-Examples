### R code from vignette source 'isocir.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: isocir.Rnw:80-81 (eval = FALSE)
###################################################
## library("isocir")


###################################################
### code chunk number 2: isocir.Rnw:444-446 (eval = FALSE)
###################################################
## data(cirdata)
## cirdata


###################################################
### code chunk number 3: isocir.Rnw:451-452 (eval = FALSE)
###################################################
## orderGroups <- c(1, 1, 1, 2, 2, 3, 4, 4)


###################################################
### code chunk number 4: isocir.Rnw:456-457 (eval = FALSE)
###################################################
## example1CIRE <- CIRE(cirdata, groups = orderGroups, circular = TRUE)


###################################################
### code chunk number 5: isocir.Rnw:461-462 (eval = FALSE)
###################################################
## example1CIRE


###################################################
### code chunk number 6: isocir.Rnw:499-501 (eval = FALSE)
###################################################
## data(datareplic) 
## orderGroups2 <- c(1:8)


###################################################
### code chunk number 7: isocir.Rnw:506-508 (eval = FALSE)
###################################################
## example2test <- cond.test(datareplic,groups=orderGroups2,biasCorrect=TRUE)
## example2test


###################################################
### code chunk number 8: isocir.Rnw:517-518 (eval = FALSE)
###################################################
## round(unlist(example2test$cirmeans), digits = 3)


###################################################
### code chunk number 9: isocir.Rnw:583-601 (eval = FALSE)
###################################################
##  data("cirgenes")
##  kappas <- c(2.64773, 3.24742, 2.15936, 4.15314, 4.54357,
##                29.07610, 6.51408, 14.19445, 5.66920, 11.12889)
## 
##  allresults <- list()
##  resultIsoCIRE <- matrix(ncol = ncol(cirgenes), nrow = nrow(cirgenes))
## 
##  SCEs <- vector(mode = "numeric", length = nrow(cirgenes))
##  pvalues <- vector(mode = "numeric", length = nrow(cirgenes))
## 
##  for (i in 1 : nrow(cirgenes)) {
##     k <- kappas[i]
##     genes <- as.numeric(cirgenes[i, !is.na(cirgenes[i, ]) ])
##     allresults[[i]] <- cond.test(genes, kappa = k)
##     resultIsoCIRE[ i, !is.na(cirgenes[i, ]) ] <- unlist(allresults[[i]]$CIRE)
##     SCEs[i] <- allresults[[i]]$SCE
##     pvalues[i] <- allresults[[i]]$pvalue
##  }



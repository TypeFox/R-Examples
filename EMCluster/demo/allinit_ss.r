library(EMCluster, quiet = TRUE)
set.seed(1234)

x <- da1$da
TC <- da1$class
lab <- da1$class

n <- nrow(x)
p <- ncol(x)
k <- 10

nk <- floor(k * 2 / 3)
tmp <- names(sort(table(lab), decreasing = TRUE))[1:nk]
lab[!(lab %in% tmp)] <- 0

for(i in 1:nk){
  tmp.id <- lab == tmp[i]
  id <- sample(which(tmp.id), sum(tmp.id) * 0.5, replace = FALSE)
  tmp.id[id] <- FALSE
  lab[tmp.id] <- 0
}
index.lab <- sort(unique(lab))
lab <- sapply(lab, function(i) which(index.lab == i) - 1)

ret.em <- init.EM(x, nclass = k, lab = lab, method = "em.EM")
ret.Rnd <- init.EM(x, nclass = k, lab = lab, method = "Rnd.EM",
                   EMC = .EMC.Rnd)
ret.Rndp <- init.EM(x, nclass = k, lab = lab, method = "Rnd.EM",
                    EMC = .EMC.Rndp)

par(mfrow = c(2, 2))
plotem(ret.em, x, main = "em")
plotem(ret.Rnd, x, main = "Rnd")
plotem(ret.Rndp, x, main = "Rnd+")


ret.all <-
cbind(
  c(ret.em$llhdval, ret.Rnd$llhdval, ret.Rndp$llhdval),
  c(RRand(ret.em$class, TC, lab = lab)$adjRand,
    RRand(ret.Rnd$class, TC, lab = lab)$adjRand,
    RRand(ret.Rndp$class, TC, lab = lab)$adjRand)
)
rownames(ret.all) <- c("em", "Rnd", "Rnd+")
colnames(ret.all) <- c("logL", "adjR")
ret.all

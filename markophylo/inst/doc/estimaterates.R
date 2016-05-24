## ----eval=FALSE----------------------------------------------------------
#  install.packages("markophylo", dependencies = TRUE, repos = "http://cran.r-project.org")

## ----eval=FALSE----------------------------------------------------------
#  install.packages(c("Rcpp","RcppArmadillo","ape","phangorn",
#  "numDeriv","knitr"), repos = "http://cran.r-project.org")
#  
#  install.packages("markophylo_1.0.2.tar.gz", repos = NULL, type = "source")

## ----eval=FALSE----------------------------------------------------------
#  install.packages(c("Rcpp","RcppArmadillo","ape","phangorn",
#  "numDeriv","knitr"), repos = "http://cran.r-project.org")
#  
#  install.packages("markophylo_1.0.2.tar.gz", repos = NULL, type = "source")

## ---- fig.show='asis',fig.align = 'center'-------------------------------
library(markophylo)
data(simdata1)

## ---- fig.show='asis',fig.align = 'center'-------------------------------
ape::plot.phylo(simdata1$tree, edge.width = 2, show.tip.label = FALSE, no.margin = TRUE)
ape::nodelabels(frame = "circle", cex = 0.7)
ape::tiplabels(frame = "circle", cex = 0.7)
print(simdata1$Q)

## ------------------------------------------------------------------------
print(table(simdata1$data))

## ------------------------------------------------------------------------
model1 <- estimaterates(usertree = simdata1$tree, userphyl = simdata1$data, 
                        alphabet = c(1, 2), rootprob = "equal", 
                        modelmat = matrix(c(NA, 1, 2, NA), 2, 2))
print(model1)

## ----eval=FALSE----------------------------------------------------------
#  model1 <- estimaterates(usertree = simdata1$tree, userphyl = simdata1$data,
#                          alphabet = c(1, 2), rootprob = "equal",
#                          modelmat = "ARD")
#  print(model1)

## ------------------------------------------------------------------------
print(class(model1))

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(model1$w, row.names = NA, col.names = c("Number of Times Pattern Observed") )

## ------------------------------------------------------------------------
filterall1 <- which(apply(simdata1$data, MARGIN = 1, FUN = 
                            function(x) isTRUE(all.equal(as.vector(x), c(1, 1, 1, 1)))))
filterall2 <- which(apply(simdata1$data, MARGIN = 1, FUN = 
                            function(x) isTRUE(all.equal(as.vector(x), c(2, 2, 2, 2)))))
filteredsimdata1 <- simdata1$data[-c(filterall1, filterall2), ]
model1_f <- estimaterates(usertree = simdata1$tree, userphyl = filteredsimdata1,
                          alphabet = c(1, 2), rootprob = "equal", 
                          modelmat = "ARD")
print(model1_f)

## ------------------------------------------------------------------------
model1_f_corrected <- estimaterates(usertree = simdata1$tree, userphyl = filteredsimdata1, 
                                    unobserved = matrix(c(1, 1, 1, 1, 2, 2, 2, 2), nrow = 2, 
                                                        byrow = TRUE), alphabet = c(1, 2), 
                        rootprob = "equal", modelmat = "ARD")
print(model1_f_corrected)

## ------------------------------------------------------------------------
data(simdata2)
print(simdata2$Q)
print(table(simdata2$data))

## ------------------------------------------------------------------------
model2 <- estimaterates(usertree = simdata2$tree, userphyl = simdata2$data, 
                        alphabet = c(1, 2), bgtype = "ancestornodes", bg = c(7),
                        rootprob = "equal", modelmat = "ARD")
print(model2)

## ---- fig.show='asis',fig.align = 'center'-------------------------------
ape::plot.phylo(simdata2$tree, edge.width = 2, show.tip.label = FALSE, no.margin = TRUE)
ape::nodelabels(frame = "circle", cex = 0.7)
ape::tiplabels(frame = "circle", cex = 0.7)

## ------------------------------------------------------------------------
model2_2 <- estimaterates(usertree = simdata2$tree, userphyl = simdata2$data, 
                        alphabet = c(1, 2), bgtype = "listofnodes", bg = list(c(3,4,7),c(1,2,5,6,7)),
                        rootprob = "equal", modelmat = "ARD")
print(model2_2)

## ---- fig.show='asis',fig.align = 'center'-------------------------------
plottree(model2, colors=c("blue", "darkgreen"), edge.width = 2, show.tip.label = FALSE, 
         no.margin = TRUE)
ape::nodelabels(frame = "circle", cex = 0.7)
ape::tiplabels(frame = "circle", cex = 0.7)

## ------------------------------------------------------------------------
data(simdata3)
print(dim(simdata3$data))
print(table(simdata3$data))
model3 <- estimaterates(usertree = simdata3$tree, userphyl = simdata3$data, 
                        alphabet = c("a", "c", "g", "t"), rootprob = "equal", 
                        partition = list(c(1:2500), c(2501:5000)), 
                        modelmat = "ER")
print(model3)

## ------------------------------------------------------------------------
data(simdata4)
print(table(simdata4$data))
model4 <- estimaterates(usertree = simdata4$tree, userphyl = simdata4$data, 
                        alphabet = c("a", "c", "g", "t"), rootprob = "maxlik",
                        ratevar = "discgamma", nocat = 4, 
                        modelmat = "ER")
print(model4)

## ------------------------------------------------------------------------
filteralla <- which(apply(simdata4$data, MARGIN = 1, FUN = 
                            function(x) isTRUE(all.equal(as.vector(x), 
                                                         c("a", "a", "a", "a")))))
filterg3 <- which(apply(simdata4$data, MARGIN = 1, FUN = 
                          function(x) table(x)["g"] >= 3) )
filteredsimdata4 <- simdata4$data[-c(filteralla, filterg3), ]
dim(simdata4$data)
dim(filteredsimdata4)

## ---- results = "hide"---------------------------------------------------
alphabet <- c("a", "c", "g", "t")
allpatt <- expand.grid(alphabet, alphabet, alphabet, alphabet) #all possible combinations
unob_patt_index_1 <- which(apply(allpatt, MARGIN = 1, FUN = function(x) table(x)["g"] >= 3) )
unob_patt_index_2 <- which(apply(allpatt, MARGIN = 1, FUN = function(x) 
  isTRUE(all.equal(as.vector(x), c("a", "a", "a", "a"))) ) )
unob_patt_indices <- sort(union(unob_patt_index_1, unob_patt_index_2)) #Ordered indices.
unob_patt <- allpatt[unob_patt_indices, ] #matrix of unique patterns

## ------------------------------------------------------------------------
model4_f_corrected <- estimaterates(usertree = simdata4$tree, userphyl = filteredsimdata4, 
                                    unobserved = unob_patt,
                        alphabet = c("a", "c", "g", "t"), rootprob = "maxlik", 
                        ratevar = "discgamma", nocat = 4, 
                        modelmat = "ER")
print(model4_f_corrected)

## ------------------------------------------------------------------------
list((1:15)[-seq(3, 15, by = 3)], seq(3, 15, by = 3) )

## ------------------------------------------------------------------------
data(simdata5)
print(dim(simdata5$data))
print(table(simdata5$data))
model5 <- estimaterates(usertree = simdata5$tree, userphyl = simdata5$data, 
                        alphabet = c("a", "c", "g", "t"), rootprob = "maxlik", 
                        partition = list((1:6000)[-seq(3, 6000, by = 3)], 
                                         seq(3, 6000, by = 3) ),
                        ratevar = "partitionspecificgamma", nocat = 4, 
                        modelmat = "ER")
print(model5)

## ------------------------------------------------------------------------
mp_data <- read.table("http://info.mcmaster.ca/~udang/mp_example_data", header = TRUE)
mp_tree <- ape::read.tree("http://info.mcmaster.ca/~udang/mp_example_tree")
mp_model <- estimaterates(usertree = mp_tree, userphyl = mp_data,
                        alphabet = c(1, 2, 3, 4), rootprob = "equal",
                        modelmat = "BD")
print(mp_model)

## ------------------------------------------------------------------------
mp_model_2 <- estimaterates(usertree = mp_tree, userphyl = mp_data,
                        alphabet = c(1, 2, 3, 4), rootprob = "maxlik",
                        modelmat = "BD")
print(mp_model_2)

## ------------------------------------------------------------------------
matrix(c(NA, 1, 1, 1, 1, NA, 1, 1, 1, 1, NA, 1, 1, 1, 1, NA
), nrow = 4, ncol = 4)

## ------------------------------------------------------------------------
matrix(c(NA, 1, 2, 3, 1, NA, 4, 5, 2, 4, NA, 6, 3, 5, 6, NA
), 4,4)

## ------------------------------------------------------------------------
matrix(c(NA, 1, 2, 3, 4, NA, 5, 6, 7, 8, NA, 9, 10, 11, 12, 
NA), 4, 4)

## ------------------------------------------------------------------------
matrix(c(NA, 1, 0, 0, 2, NA, 1, 0, 0, 2, NA, 1, 0, 0, 2, NA
), 4, 4)

## ------------------------------------------------------------------------
matrix(c(NA, 1, 0, 0, 1, NA, 1, 0, 0, 1, NA, 1, 0, 0, 1, NA
), 4, 4)

## ------------------------------------------------------------------------
matrix(c(NA, 1, 0, 0, 1, NA, 2, 0, 0, 2, NA, 3, 0, 0, 3, NA
), 4, 4)

## ------------------------------------------------------------------------
matrix(c(NA, 1, 0, 0, 4, NA, 2, 0, 0, 5, NA, 3, 0, 0, 6, NA
), 4, 4)

## ------------------------------------------------------------------------
m <- matrix(c(NA, 1, 2, 0, NA, 1, 0, 0, NA), 3, 3)
rownames(m) <- colnames(m) <- c("Absent", "1 copy", "2 copies")
print(m)

## ------------------------------------------------------------------------
m <- matrix(0, 5, 5)
diag(m) <- NA
diag(m[-1, ]) <- 1
diag(m[, -1]) <- 2
print(m)


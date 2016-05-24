library(gnm)
set.seed(1)

library(nnet)
.incidence <- class.ind(backPain$pain)
.counts <- as.vector(t(.incidence))
.rowID <- factor(t(row(.incidence)))
backPain <- backPain[.rowID, ]
backPain$pain <- C(factor(rep(levels(backPain$pain), nrow(.incidence)),
                          levels = levels(backPain$pain), ordered = TRUE),
                   treatment)

noRelationship <- gnm(.counts ~ pain, eliminate = .rowID,
                      family = "poisson", data = backPain)

oneDimensional <- update(noRelationship,
                         ~ . + Mult(pain, x1 + x2 + x3))

print(oneDimensional$deviance, digits=10)
print(oneDimensional$df)

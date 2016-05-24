library(Rarity)
data(spid.occ)

# Vector
o <- spid.occ[1:50, 1]
names(o) <- rownames(spid.occ)[1:50]
rWeights(o)

# with NAs
o[3:5] <- NA
rWeights(o)

# matrix
o <- spid.occ[1:50, ]
rWeights(o)

# with NAs
o[3, 1] <- NA
o[4, 2] <- NA
o[5, ] <- NA
rWeights(o)
rWeights(o, extended = TRUE)
rWeights(o, wMethods = c("W", "invQ", "oldW"), extended = TRUE)
assemblages.matrix <- cbind(assemblage.1 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.2 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.3 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.4 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.5 = sample(c(0, 1), 50, replace = TRUE))
rownames(assemblages.matrix) <- rownames(o)
rWeights(o, wMethods = c("W", "invQ", "oldW"),  rCutoff = "Leroy", assemblages = assemblages.matrix, extended = TRUE)
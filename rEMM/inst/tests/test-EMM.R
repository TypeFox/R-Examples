library("rEMM")
library("testthat")
data("16S")

data <- Mollicutes16S+1
test <- Mollicutes16S[2:10,]+1

## create two EMMs for different data
emm <- EMM("Kullback", threshold=0.1, data=data)
emm



## TRACDS
expect_identical(nstates(emm), length(states(emm)))

expect_is(current_state(emm), "character")
transitions(emm)
rare_transitions(emm, 1)

## TRACDS
cluster_counts(emm)
expect_equivalent(nrow(cluster_centers(emm)), nclusters(emm))
expect_equivalent(nclusters(emm), nstates(emm))
expect_equivalent(clusters(emm), rownames(cluster_centers(emm)))
expect_is(last_clustering(emm), "character")
expect_equivalent(nrow(data), length(last_clustering(emm)))

rare_clusters(emm, 1)
find_clusters(emm, test)

## score, predict et al
expect_equivalent(nrow(transition_table(emm, test)), nrow(test)-1L)
transition_table(emm, test, prior=FALSE)

score(emm, test)
score(emm, test, prior=FALSE)
score(emm, test, method="sum")
score(emm, test, method="sum", prior=FALSE)

predict(emm, "1")
p <- predict(emm, "1", probabilities=TRUE)
p[p>0]
p <- predict(emm, "1", probabilities=TRUE, prior=FALSE)
table(p)

## reset
reset(emm)
expect_identical(current_state(emm), NA_character_)

## copy: these need to be false!
expect_true(!identical(emm@tnn_d, copy(emm)@tnn_d))
expect_true(!identical(emm@tracds_d, copy(emm)@tracds_d))


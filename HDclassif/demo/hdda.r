#Supervised classification of the "wine" dataset.
#The data is scaled using the command scaling=TRUE.
#The learning is done using a random sample of 40 individuals
#whereas the testing is done with the 138 remaining individuals.
#The graph of the choice of the intrinsic dimensions is shown.

data(wine)
w <- wine[, -1]
cls <- wine[, 1]
set.seed(1)
ind <- sample(178, 40)
prms <- hdda(w[ind, ], cls[ind], scaling = TRUE, model = "all", graph = TRUE)
res <- predict(prms, w[-ind, ], cls[-ind])
plot(prms)
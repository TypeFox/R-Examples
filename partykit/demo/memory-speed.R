### packages and data
library("rpart")
library("RWeka")
library("partykit")
data("Shuttle", package = "mlbench")

### fit rpart and J48 trees
rp <- rpart(Class ~ ., data = Shuttle)
j48 <- J48(Class ~ ., data = Shuttle)

### convert to party
system.time(party_rp <- as.party(rp))
system.time(party_j48 <- as.party(j48))

### check depth/width
depth(party_rp)
width(party_rp)
depth(party_j48)
width(party_j48)

### compare object sizes
osize <- function(x) print(object.size(x), units = "Kb")
osize(rp)                    ## rpart representation
osize(party_rp)              ## full party (with terms, fitted values)
osize(node_party(party_rp))  ## only the raw partynode
osize(j48)                   ## J48 tree in external Java pointer
osize(party_j48)             ## full party (with terms, fitted values)
osize(node_party(party_j48)) ## only the raw partynode
osize(Shuttle)               ## learning data (not stored in any tree)

### set-up large prediction sample
set.seed(1)
nd <- Shuttle[sample(1:nrow(Shuttle), 1e6, replace = TRUE), ]

### compare predictions (speed and accuracy)
system.time(p_rp <- predict(rp, newdata = nd, type = "prob"))
system.time(p_party_rp <- predict(party_rp, newdata = nd, type = "prob"))
all.equal(p_rp, p_party_rp)

system.time(p_j48 <- predict(j48, newdata = nd))
system.time(p_party_j48 <- predict(party_j48, newdata = nd))
all.equal(p_j48, p_party_j48, check.attributes = FALSE)

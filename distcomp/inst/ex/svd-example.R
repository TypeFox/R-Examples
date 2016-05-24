library(distcomp)

print(availableComputations())

svdDef <- data.frame(compType = names(availableComputations())[2],
                     rank = 2L,
                     ncol = 5L,
                     id = "SVD",
                     stringsAsFactors = FALSE)

library(opencpu)

set.seed(12345)

## Three sites
nSites <- 3
siteData <- lapply(seq.int(nSites), function(i) matrix(rnorm(100), nrow=20))
sites <- lapply(seq.int(nSites),
                function(x) list(name = paste0("site", x),
                                 url = opencpu$url()))

## Upload definition
ok <- Map(uploadNewComputation, sites,
          lapply(seq.int(nSites), function(i) svdDef),
          siteData)

stopifnot(all(as.logical(ok)))

##
master <- SVDMaster$new(defnId = svdDef$id, k = svdDef$rank)

for (site in sites) {
  master$addSite(name = site$name, url = site$url)
}

result <- master$run()

print(result)

## Compare with:

x <- do.call(rbind, siteData)

print(result$d)
print(svd(x)$d)

print(result$v)
print(svd(x)$v[, 1:2])

## All  singular values (takes much longer!)
result <- master$run(k=5)
print(result$d)
print(svd(x)$d)

print(result$v)
print(svd(x)$v)


sessionInfo()

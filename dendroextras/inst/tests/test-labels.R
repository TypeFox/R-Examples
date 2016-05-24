context("Test handling of dendrogram labels")

test_that("Find labels of hclust object", {
      hc <- hclust(dist(USArrests), "ave")
      dend <- as.dendrogram(hc)
      expect_that(labels(hc), equals(labels(dend)))
    })

test_that("Assign labels of dendrogram object", {
      hc <- hclust(dist(USArrests), "ave")
      dend <- as.dendrogram(hc)
      new.labels=abbreviate(labels(dend))
      labels(dend)<-new.labels
      expect_that(labels(dend), equals(new.labels))
    })


context("Test cutting trees")

test_that("Ordered cut of an hclust object", {
      hc <- hclust(dist(USArrests), "ave")
      cut.ordered.k5 = structure(c(1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 
              3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 
              4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5), .Names = c("Florida", "North Carolina", 
              "California", "Maryland", "Arizona", "New Mexico", "Delaware", 
              "Alabama", "Louisiana", "Illinois", "New York", "Michigan", "Nevada", 
              "Alaska", "Mississippi", "South Carolina", "Washington", "Oregon", 
              "Wyoming", "Oklahoma", "Virginia", "Rhode Island", "Massachusetts", 
              "New Jersey", "Missouri", "Arkansas", "Tennessee", "Georgia", 
              "Colorado", "Texas", "Idaho", "Nebraska", "Kentucky", "Montana", 
              "Ohio", "Utah", "Indiana", "Kansas", "Connecticut", "Pennsylvania", 
              "Hawaii", "West Virginia", "Maine", "South Dakota", "North Dakota", 
              "Vermont", "Minnesota", "Wisconsin", "Iowa", "New Hampshire"))
      expect_that(slice(hc,k=5),equals(cut.ordered.k5),"Cut into k=5 groups")
    })

test_that("Compare group memberships with those returned by rect.hclust", {
      hc <- hclust(dist(USArrests), "ave")
      plot(hc)
      # find groups with rect.hclust (returns lists of named numeric ids of leaves)
      # list order does match dendrogram order for groups
      rhc=rect.hclust(hc,k=5)
      # Unlist and return as a single vector of named groups
      cutbyrect.hclust=structure(rep(seq(rhc),sapply(rhc,length)),
          .Names=names(unlist(rhc)))
      # reorder, first so that we are in original order 
      # and then in dendrogram order  
      cutbyrect.hclust=cutbyrect.hclust[order(unlist(rhc))][hc$order]
      expect_that(slice(hc,k=5),equals(cutbyrect.hclust),"Compare with rect.hclust")
    })

test_that("Compare slice.hclust and slice.dendrogram", {
      hc <- hclust(dist(USArrests), "ave")
      dend <- as.dendrogram(hc)
      slice_hc=slice(hc,k=5)
      expect_that(slice(dend,k=5),equals(slice_hc),"Compare slicing hclust and dendrogram objects for k=5")
    })

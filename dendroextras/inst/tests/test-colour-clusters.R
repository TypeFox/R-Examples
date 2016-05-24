context("Test colouring of dendrograms")

test_that("Can find leaf colours", {
      hc <- hclust(dist(USArrests), "ave")
      dend <- as.dendrogram(hc)
      colours5=c("red",'green','blue','cyan','magenta')
      cdendk5 <- colour_clusters(dend,k=5,col=colours5)
      leafcolours <- unlist(dendrapply(cdendk5,function(n) 
                if(is.leaf(n)) structure(attr(n,'edgePar')$col,
                      .Names=attr(n,'label')) else NULL))
      
      expect_that(leaf_colours(cdendk5), equals(leafcolours))
      no_colours <- structure(rep(NA_character_,length(leafcolours)),
          .Names=names(leafcolours))
      expect_that(leaf_colours(dend),equals(no_colours),
          'no edge colours in dendrogram')
      expect_that(leaf_colours(cdendk5,'label'),equals(no_colours),
          'no label colours in dendrogram')
      expect_that(leaf_colours(cdendk5,'node'),equals(no_colours),
          'no node colours ')
    })

test_that("Colour a dendrogram by cluster identity", {
      hc <- hclust(dist(USArrests), "ave")
      dend <- as.dendrogram(hc)
      cdend2 <- colour_clusters(dend,k=2,col=c("red",'green'))
      leafcolours <- leaf_colours(cdend2)
      expect_that(labels(hc), equals(names(leafcolours)))
      expect_equivalent(rep(c("red",'green'),c(16,34)),leafcolours)
      
      colours5=c("red",'green','blue','cyan','magenta')
      cdend5 <- colour_clusters(dend,k=5,col=colours5)
      leafcolours <- leaf_colours(cdend5)
      expect_that(labels(hc), equals(names(leafcolours)))
      expect_equivalent(rep(colours5,c(2L,14L,14L,10L,10L)),leafcolours)
    })

test_that("Colouring a dendrogram by cutting on height or number of groups", {
      hc <- hclust(dist(USArrests), "ave")
      dend <- as.dendrogram(hc)
      colours5=c("red",'green','blue','cyan','magenta')
      cdendk5 <- colour_clusters(dend,k=5,col=colours5)
      cdendh50 <- colour_clusters(dend,h=50,col=colours5)
      
      expect_that(cdendk5, equals(cdendh50))
    })

test_that("Can cut/colour an hclust object returning a dendrogram", {
      hc <- hclust(dist(USArrests), "ave")
      dend <- as.dendrogram(hc)
      colours5=c("red",'green','blue','cyan','magenta')
      cdendk5 <- colour_clusters(dend,k=5,col=colours5)
      chck5 <- colour_clusters(hc,k=5,col=colours5)
      
      expect_that(cdendk5, equals(chck5))
    })

test_that("Can colour a dendrogram and add group labels", {
      hc <- hclust(dist(USArrests), "ave")
      dend <- as.dendrogram(hc)
      colours5=c("red",'green','blue','cyan','magenta')
      cdendk51 <- colour_clusters(dend,k=5,col=colours5,groupLabels=TRUE)
      cdendk52 <- colour_clusters(dend,k=5,col=colours5,groupLabels=1:5)
      cdendk53 <- colour_clusters(dend,k=5,col=colours5,groupLabels=LETTERS[1:5])
      cdendk54 <- colour_clusters(dend,k=5,col=colours5,groupLabels=as.roman)
    })

test_that("Can colour a dendrogram with non-unique labels",{
  dmat=structure(c(4.66, 2.72, 2.75, 2.71, 4.21, 2.83, 3.45, 2.26, 3.4, 
                   2.68, 0.74, 3, 2.32, 5.09, 2.96), 
                 Labels = c("A", "B", "C", "D", "E", "F"), Size = 6L, 
                 class = "dist", Diag = FALSE, Upper = FALSE)
  
  dmat2=structure(c(4.66, 2.72, 2.75, 2.71, 4.21, 2.83, 3.45, 2.26, 3.4, 
              2.68, 0.74, 3, 2.32, 5.09, 2.96), 
            Labels = c("A", "A", "B", "A", "B", "A"), Size = 6L, 
            class = "dist", Diag = FALSE, Upper = FALSE)
  hc=hclust(dmat)
  hc2=hclust(dmat2)
  plot(hc)
  expect_is(hcd<-colour_clusters(hc, k=2), 'dendrogram')
  expect_is(hcd2<-colour_clusters(hc2, k=2), 'dendrogram')
  expect_equal(unname(leaf_colors(hcd)),unname(leaf_colors(hcd2)))  
})

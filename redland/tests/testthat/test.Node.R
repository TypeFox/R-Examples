context("Node tests")
test_that("redland library loads", {
  library(redland)
})
test_that("Node constructor", {
  library(redland)
  world <- new("World")
  expect_false(is.null(world))
  
  # Test creating a blank node
  node <- new("Node", world);
  expect_false(is.null(node))
  expect_that(class(node), matches("Node"))
  expect_that(getNodeType(node), matches("blank"))
  
  # Test creating a blank node with librdf generated identifier
  node <- new("Node", world, blank=NULL)
  expect_that(class(node@librdf_node), matches("_p_librdf_node_s")) 
  expect_that(getNodeType(node), matches("blank")) 

  # Test creating a blank node with librdf generated identifier
  node <- new("Node", world, blank="_:1234")
  expect_that(class(node@librdf_node), matches("_p_librdf_node_s")) 
  expect_that(getNodeType(node), matches("blank"))
  
  # Test creating a node with a literal value
  node <- new("Node", world, "Fee Fi Fo Fum")
  expect_that(class(node@librdf_node), matches("_p_librdf_node_s")) 
  expect_that(getNodeType(node), matches("literal"))
  #   err <- try(freeNode(node), silent=TRUE)
  #   expect_that(class(err), not(matches("try-error")))
  
  # Test creating a node with a literal value
  node <- new("Node", world, literal="Fee Fi Fo Fum")
  expect_that(class(node@librdf_node), matches("_p_librdf_node_s"))
  expect_that(getNodeType(node), matches("literal"))
  
  # Test creating a node with a URI value
  node <- new("Node", world, uri="http://www.example.com/")
  expect_that(class(node@librdf_node), matches("_p_librdf_node_s"))
  expect_that(getNodeType(node), matches("resource"))
  
  # Test that node creation fails if world is not provided or is null
  err <- try(node <- new("Node", literal="Fee Fi Fo Fum"), silent=TRUE)
  expect_that(class(err), matches("try-error"))
  err <- try(node <- new("Node", world=NULL, literal="Fee Fi Fo Fum"), silent=TRUE)
  expect_that(class(err), matches("try-error"))

})
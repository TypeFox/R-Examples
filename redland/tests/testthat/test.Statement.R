context("Statement tests")
test_that("redland library loads", {
  library(redland)
})
test_that("Statement constructor", {
  library(redland)
  world <- new("World")
  expect_false(is.null(world))
  
  # Test creating Subject, predicate, and object Nodes
  subject <- new("Node", world, literal="subject")
  expect_that(class(subject@librdf_node), matches("_p_librdf_node_s"))
  predicate <- new("Node", world, literal="subject")
  expect_that(class(predicate@librdf_node), matches("_p_librdf_node_s"))
  object <- new("Node", world, literal="subject")
  expect_that(class(object@librdf_node), matches("_p_librdf_node_s"))
  
  # Test creating the Statement
  stmt <- new("Statement", world, subject, predicate, object)
  expect_false(is.null(stmt))
  expect_that(class(stmt@librdf_statement), matches("_p_librdf_statement_s"))
  
  # Test that statement creation fails if world is not provided or is null
  err <- try(stmt <- new("Statement", world=NULL, subject=subject, predicate=predicate, object=object), silent=TRUE)
  expect_that(class(err), matches("try-error"))
  
  # Test that statement creation fails if subject, predicate, or object is not provided or is null
  err <- try(stmt <- new("Statement", world=world, subject=NULL, predicate=predicate, object=object), silent=TRUE)
  expect_that(class(err), matches("try-error"))
  err <- try(stmt <- new("Statement", world=world, subject=subject, predicate=NULL, object=object), silent=TRUE)
  expect_that(class(err), matches("try-error"))
  err <- try(stmt <- new("Statement", world=world, subject=subject, predicate=predicate, object=NULL), silent=TRUE)
  expect_that(class(err), matches("try-error"))
  
  # Test statement creation when subject, predicate, object are passed in as character and RDF type is not specified
  stmt <- new("Statement", world=world, 
              subject="https://cn.dataone.org/cn/v1/resolve/resourceMap_5cbcdecd-6b0e-4b24-a0be-20291b2e49a7#aggregation",
              predicate="http://purl.org/dc/terms/identifier",
              object="resourceMap_5cbcdecd-6b0e-4b24-a0be-20291b2e49a7^^xsd:string")
              
  expect_that(getTermType(stmt, "subject"), matches("resource"))
  expect_that(getTermType(stmt, "predicate"), matches("resource"))
  expect_that(getTermType(stmt, "object"), matches("literal"))
  
  # Test 
  stmt <- new("Statement", world=world, 
              subject="_:foo1",
              predicate="http://purl.org/dc/terms/identifier",
              object=NULL)
  
  expect_that(getTermType(stmt, "subject"), matches("blank"))
  expect_that(getTermType(stmt, "predicate"), matches("resource"))
  expect_that(getTermType(stmt, "object"), matches("blank"))
  
  err <- try(freeStatement(stmt), silent=TRUE)
  expect_false(class(err) == "try-error")
  
  stmt <- new("Statement", world=world, 
              subject=NULL,
              predicate="http://purl.org/dc/terms/identifier",
              object="id1234")
  
  expect_that(getTermType(stmt, "subject"), matches("blank"))
  expect_that(getTermType(stmt, "predicate"), matches("resource"))
  expect_that(getTermType(stmt, "object"), matches("literal"))
  
  err <- try(freeStatement(stmt), silent=TRUE)
  expect_false(class(err) == "try-error")
  
  # Test statement creation when subject, predicate and object are passed in as charater and RDF types are specified
  stmt <- new("Statement", world=world, 
              subject="http://www.exmaple.com/subject",
              predicate="http://purl.org/dc/terms/identifier",
              object="http://www.exmaple.com/object", subjectType="blank", objectType="literal")
  
  expect_that(getTermType(stmt, "subject"), matches("blank"))
  expect_that(getTermType(stmt, "predicate"), matches("resource"))
  expect_that(getTermType(stmt, "object"), matches("literal"))
  
  err <- try(freeStatement(stmt), silent=TRUE)
  expect_false(class(err) == "try-error")
  err <- try(freeWorld(world), silent=TRUE)
  expect_false(class(err) == "try-error")
})
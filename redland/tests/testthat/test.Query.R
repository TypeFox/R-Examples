context("Query tests")
test_that("redland library loads", {
  library(redland)
})
test_that("Query works", {
  library(redland)
  world <- new("World")
  expect_false(is.null(world))
  
  # Test creating the Storage system
  storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
  expect_false(is.null(storage))
  expect_that(class(storage@librdf_storage), matches("_p_librdf_storage_s"))
  
  # Test creating the Model
  model <- new("Model", world, storage, options="")
  expect_false(is.null(model))
  expect_that(class(model@librdf_model), matches("_p_librdf_model_s"))
  
  # Test that model creation fails if world is not provided or is null
  err <- try(model <- new("Model", world=NULL, storage, options=""), silent=TRUE)
  expect_that(class(err), matches("try-error"))
  
  expect_false(is.null(model))
  expect_that(class(model@librdf_model), matches("_p_librdf_model_s"))
  
  # Add typical statements to the model
  stmt <- new("Statement", world=world, 
              subject="https://cn.dataone.org/cn/v1/resolve/urn:uuid:1ef9ff61-c4d6-4314-9ffd-cb8b8fd71e5c",
              predicate="http://www.w3.org/ns/prov#wasGeneratedBy",
              object="https://cn.dataone.org/cn/v1/resolve/urn:uuid:10153c9c-0cf4-4d92-9544-8e14089df8d0")
  addStatement(model, stmt)
  
  stmt <- new("Statement", world=world,
              subject="https://cn.dataone.org/cn/v1/resolve/urn:uuid:1ef9ff61-c4d6-4314-9ffd-cb8b8fd71e5c",
              predicate="http://purl.org/dc/terms/identifier",
              object="urn:uuid:1ef9ff61-c4d6-4314-9ffd-cb8b8fd71e5c", objectType="literal", datatype_uri="http://www.w3.org/2001/XMLSchema#string")
  addStatement(model, stmt)
  
  stmt <- new("Statement", world=world, 
              subject="https://cn.dataone.org/cn/v1/resolve/urn:uuid:274a0c5c-3082-4562-bbd3-2b1288768cac",
              predicate="http://www.w3.org/ns/prov#hadPlan",
              object="https://cn.dataone.org/cn/v1/resolve/urn:uuid:01305f45-f22b-40c8-8d27-00357d01e4a5")
  addStatement(model, stmt)

  stmt <- new("Statement", world=world, 
              subject="https://orcid.org/0000-0002-2192-403X",
              predicate="http://www.w3.org/ns/prov#Agent",
              object="slaughter", 
              objectType="literal", datatype_uri="http://www.w3.org/2001/XMLSchema#string")
  addStatement(model, stmt)
  
  stmt <- new("Statement", world=world,
              subject="https://cn.dataone.org/cn/v1/resolve/urn:uuid:1ef9ff61-c4d6-4314-9ffd-cb8b8fd71e5c",
              predicate="http://www.w3.org/ns/prov#qualifiedAssociation",
              object="https://cn.dataone.org/cn/v1/resolve/urn:uuid:274a0c5c-3082-4562-bbd3-2b1288768cac")
  addStatement(model, stmt)

# Query the RDF model with a SPARQL query that should return all triples
queryString <- 'PREFIX orcid: <https://orcid.org/> PREFIX dataone: <https://cn.dataone.org/cn/v1/resolve/> PREFIX prov: <http://www.w3.org/ns/prov#> SELECT ?a ?b ?c WHERE { ?a ?b ?c . }'
query <- new("Query", world, queryString, base_uri=NULL, query_language="sparql", query_uri=NULL)
expect_that(class(query), matches("Query"))
queryResult <- executeQuery(query, model)
expect_that(class(queryResult), matches("QueryResult"))

# Retrieve query results and check the actual result count against the expected count
result <- getNextResult(queryResult)
i <- 0
while(!is.null(result)) {
  i <- i + 1
  # Something went wrong, break loop
  if(i > 5) {
    break
  }
  expect_that(class(result), matches("list"))
  result <- getNextResult(queryResult)
}
expect_that(i, equals(5))

freeQuery(query)
rm(query)
freeQueryResults(queryResult)
rm(queryResult)

# Query the RDF model with a new query that should only have one triple returned
queryString <- 'PREFIX orcid: <https://orcid.org/> PREFIX dataone: <https://cn.dataone.org/cn/v1/resolve/> PREFIX prov: <http://www.w3.org/ns/prov#> SELECT ?a ?c WHERE { ?a prov:Agent ?c . }'
query <- new("Query", world, queryString, base_uri=NULL, query_language="sparql", query_uri=NULL)
expect_that(class(query), matches("Query"))
queryResult <- executeQuery(query, model)
expect_that(class(queryResult), matches("QueryResult"))

# Retrieve query results and check the actual result count against the expected count
result <- getNextResult(queryResult)
i <- 0
while(! is.null(result)) {
  i <- i + 1
  # Something went wrong, break loop
  if(i > 5) {
    break
  }
  expect_that(class(result), matches("list"))
  result <- getNextResult(queryResult)
}
expect_that(i, equals(1))

freeQuery(query)
rm(query)
freeQueryResults(queryResult)
rm(queryResult)
  
})
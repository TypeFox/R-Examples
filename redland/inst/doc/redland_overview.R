## ------------------------------------------------------------------------
library(redland)
world <- new("World")
storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
model <- new("Model", world=world, storage, options="")
parser <- new("Parser", world)
parseFileIntoModel(parser, world, system.file("extdata", "dc.rdf", package="redland"), model)
  

## ------------------------------------------------------------------------
queryString <- 'PREFIX dc: <http://purl.org/dc/elements/1.1/> SELECT ?a ?c WHERE { ?a dc:creator ?c . }'
query <- new("Query", world, queryString, base_uri=NULL, query_language="sparql", query_uri=NULL)
queryResult <- executeQuery(query, model)
result <- getNextResult(queryResult)
cat(sprintf("Result from query: %s\n", result))

## ------------------------------------------------------------------------
stmt <- new("Statement", world=world, 
        subject="http://www.dajobe.org/",
        predicate="http://purl.org/dc/elements/1.1/language",
        object="en")
addStatement(model, stmt)

## ------------------------------------------------------------------------
serializer <- new("Serializer", world, mimeType="application/rdf+xml")
status <- setNameSpace(serializer, world, namespace="http://purl.org/dc/elements/1.1/", prefix="dc")  
filePath <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".rdf")
status <- serializeToFile(serializer, world, model, filePath)


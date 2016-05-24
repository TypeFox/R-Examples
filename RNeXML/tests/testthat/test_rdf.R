context("rdf")

test_that("we can extract rdf-xml", {

  if(require("Sxslt")){
    f <- system.file("examples", "meta_example.xml", package="RNeXML")
    rdf <- get_rdf(f)
    expect_is(rdf, "XMLInternalXSLTDocument")
  }
})

test_that("we can perform sparql queries with rrdf", {
  skip_on_travis()

  if(require("Sxslt")){
    f <- system.file("examples", "meta_example.xml", package="RNeXML")
    rdf <- get_rdf(f)

## Write to a file and read in with rrdf
    saveXML(rdf, "rdf_meta.xml")
    success <- require(rrdf)
    if(success){
      lib <- load.rdf("rdf_meta.xml")
## Perform a SPARQL query:
      out <- sparql.rdf(lib, "SELECT ?title WHERE { ?x <http://purl.org/dc/elements/1.1/title> ?title}")
    }
    unlink("rdf_meta.xml")
  }
})

context("librdf base API tests")
test_that("redland library loads", {
    library(redland)
})
test_that("librdf basic functions", {
    library(redland)
    world <- librdf_new_world();
    expect_that(class(world), matches("_p_librdf_world_s"))
    storage <- librdf_new_storage(world,'hashes','dummy',"new=yes,hash-type='memory'")
    expect_that(class(storage), matches("_p_librdf_storage_s"))

    model <- librdf_new_model(world,storage,'')
    expect_that(class(model), matches("_p_librdf_model_s"))

    parser <- librdf_new_parser(world,'rdfxml','application/rdf+xml',NULL)
    expect_that(class(parser), matches("_p_librdf_parser_s"))
    uri <- librdf_new_uri(world,paste0('file:',system.file('extdata/dc.rdf', package='redland')))
    expect_that(class(uri), matches("_p_librdf_uri_s"))
    rv <- librdf_parser_parse_into_model(parser,uri,uri,model)
    expect_that(rv, equals(0))
    librdf_free_uri(uri);
    librdf_free_parser(parser);

    query <- librdf_new_query(world, 'sparql', NULL, "PREFIX dc: <http://purl.org/dc/elements/1.1/> SELECT ?a ?c ?d WHERE { ?a dc:title ?c . OPTIONAL { ?a dc:related ?d } }", NULL)
    results <- librdf_model_query_execute(model, query);
    expect_false(is.null(results))
    expect_that(class(results), matches("_p_librdf_query_results"))

    # Convert the whole sparql result to a string and check its value
    qstr <- librdf_query_results_to_string(results, NULL, NULL)
    expect_that(qstr, matches("http://www.dajobe.org/"))
    expect_that(qstr, matches("Beckett"))

    # Re-execute the query, and gather each bound node and check its value
    results <- librdf_model_query_execute(model, query);
    while (!is.null(results) && librdf_query_results_finished(results) == 0) {
        num_nodes <- librdf_query_results_get_bindings_count(results)
        for (i in 1:num_nodes-1) {
            binding_name <- librdf_query_results_get_binding_name(results, i)
            val = librdf_query_results_get_binding_value(results, i)
            expect_that(class(val), matches("_p_librdf_node_s"))
            if (!is.null.externalptr(val@ref)) {
                nval <- librdf_node_to_string(val)
                expect_that(class(nval), matches("character"))
                expect_that(nval, matches("[a-z]"))
            } else {
                nval = "(unbound)"
            }
            if (i == 2) {
                expect_that(nval, matches("(unbound)"))
            }
        }
        rc <- librdf_query_results_next(results)
    }

    # Test adding a new Statement to the model
    about <- "http://matt.magisa.org/"
    subject <- librdf_new_node_from_uri_string(world, about)
    expect_false(is.null(subject))
    expect_that(class(subject), matches("_p_librdf_node_s"))
    pred <- librdf_new_node_from_uri_string(world, "http://purl.org/dc/elements/1.1/title")
    expect_false(is.null(pred))
    expect_that(class(pred), matches("_p_librdf_node_s"))
    object <- librdf_new_node_from_literal(world, "Matt Jones' Home Page", "", 0)
    expect_false(is.null(object))
    expect_that(class(object), matches("_p_librdf_node_s"))
    statement <- librdf_new_statement_from_nodes(world, subject, pred, object)
    expect_false(is.null(statement))
    expect_that(class(statement), matches("_p_librdf_statement_s"))

    rc <- librdf_model_add_statement(model, statement)
    expect_that(rc, equals(0))

    # Don't need to call librdf_free_node, just librdf_free_statement.
    librdf_free_statement(statement)

    # Test serialization of the model to a text file
    serializer <- librdf_new_serializer(world, "rdfxml", "", NULL);
    expect_false(is.null(serializer))
    expect_that(class(serializer), matches("_p_librdf_serializer"))
    base = librdf_new_uri(world, "http://example.org/base.rdf");
    filePath <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".rdf")

    librdf_serializer_serialize_model_to_file(serializer,filePath,base,model);
    expect_that(file.exists(filePath), is_true())
    unlink(filePath)

    # Free resources
    librdf_free_serializer(serializer);
    librdf_free_uri(base);
    librdf_free_model(model);
    librdf_free_storage(storage);
    expect_that("Reached end.", equals("Reached end."))
})

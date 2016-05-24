context("solr - methods")

# * d_solr_facet
# * d_solr_group
# * d_solr_highlight
# * d_solr_mlt
# * aa <- d_solr_search
# * d_solr_stats

test_that("d_solr_search works", {
  skip_on_cran()

  aa <- d_solr_search(q="Galliard", verbose = FALSE)

  # Basic search, restricting to certain fields
  bb <- d_solr_search(q="Galliard", fl='handle,dc.title_sort', verbose = FALSE)

  # Search all text for a string, but limits results to two specified fields:
  cc <- d_solr_search(q="dwc.ScientificName:drosophila", fl='handle,dc.title_sort', verbose = FALSE)

  # Dryad data based on an article DOI:
  dd <- d_solr_search(q="dc.relation.isreferencedby:10.1038/nature04863",
     fl="dc.identifier,dc.title_ac", verbose = FALSE)

  # All terms in the dc.subject facet, along with their frequencies:
  ee <- d_solr_facet(q="location:l2", facet.field="dc.subject_filter", facet.minCount=1,
     facet.limit=10, verbose = FALSE)

  # Article DOIs associated with all data published in Dryad over the past 90 days:
  ff <- d_solr_search(q="dc.date.available_dt:[NOW-90DAY/DAY TO NOW]",
     fl="dc.relation.isreferencedby", rows=10, verbose = FALSE)

  # Data DOIs published in Dryad during January 2011, with results returned in JSON format:
  query <- "location:l2 dc.date.available_dt:[2011-01-01T00:00:00Z TO 2011-01-31T23:59:59Z]"
  gg <- d_solr_search(q=query, fl="dc.identifier", rows=10, verbose = FALSE)

  expect_is(aa, "data.frame")
  expect_is(bb, "data.frame")
  expect_is(cc, "data.frame")
  expect_is(dd, "data.frame")
  expect_is(ee, "list")
  expect_is(ff, "data.frame")
  expect_is(gg, "data.frame")

  expect_true(any(grepl("Galliard", aa$dc.contributor.author, ignore.case = TRUE)))
  expect_more_than(NROW(aa), 0)

  expect_named(bb, c('handle', 'dc.title_sort'))

  expect_true(any(grepl("Drosophila", cc$dc.title_sort, ignore.case = TRUE)))
  expect_named(cc, c('handle', 'dc.title_sort'))

  expect_named(dd, c('dc.identifier', 'dc.title_ac'))

  expect_named(ee, c('facet_queries','facet_fields','facet_dates','facet_ranges'))
  expect_is(ee$facet_fields$dc.subject_filter, 'data.frame')

  expect_named(ff, 'dc.relation.isreferencedby')
  expect_true(all(grepl("doi", ff$dc.relation.isreferencedby)))

  expect_named(gg, 'dc.identifier')
  expect_true(all(grepl("doi", gg$dc.identifier)))
})

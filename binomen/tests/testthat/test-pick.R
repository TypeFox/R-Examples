context("pick")

out <- make_taxon(genus="Poa", epithet="annua", authority="L.",
                  family='Poaceae', clazz='Poales', kingdom='Plantae', variety='annua')

test_that("operating on taxon objects", {
  aa <- out %>% pick(family)
  bb <- out %>% pick(genus)
  cc <- out %>% pick(species, genus)

  expect_is(aa, "taxon")
  expect_is(aa$binomial, "binomial")
  expect_is(aa$grouping, "grouping")
  expect_is(aa$grouping$family, "taxonref")
  expect_equal(names(aa$grouping), "family")

  expect_is(bb, "taxon")
  expect_is(bb$binomial, "binomial")
  expect_is(bb$grouping, "grouping")
  expect_is(bb$grouping$genus, "taxonref")
  expect_equal(names(bb$grouping), "genus")

  expect_is(cc, "taxon")
  expect_is(cc$binomial, "binomial")
  expect_is(cc$grouping, "grouping")
  expect_is(cc$grouping$genus, "taxonref")
  expect_is(cc$grouping$species, "taxonref")
  expect_equal(names(cc$grouping), c('genus', 'species'))
})

test_that("parts can be picked out by name", {
  aa <- out %>% pick(species) %>% name()
  bb <- out %>% pick(species) %>% uri()

  expect_is(aa, "character")
  expect_equal(aa, "Poa annua")
  expect_is(bb, "character")
  expect_equal(bb, "none")
})

df <- data.frame(class=c('Magnoliopsida','Magnoliopsida','Magnoliopsida',
    'Magnoliopsida','Magnoliopsida','Magnoliopsida'),
  order=c('Asterales','Asterales','Fagales','Poales','Poales','Poales'),
  family=c('Asteraceae','Asteraceae','Fagaceae','Poaceae','Poaceae','Poaceae'),
  genus=c('Helianthus','Helianthus','Quercus','Poa','Festuca','Holodiscus'),
  stringsAsFactors = FALSE)
df2 <- taxon_df(df)

test_that("operating on taxonomic data.frames - select single or many taxonomic classes", {
  aa <- df2 %>% pick(order)
  bb <- df2 %>% pick(family, genus)

  expect_is(aa, "data.frame")
  expect_equal(NCOL(aa), 1)
  expect_named(aa, "order")

  expect_is(bb, "data.frame")
  expect_equal(NCOL(bb), 2)
  expect_named(bb, c('family', 'genus'))
})

test_that("from taxa objects, via scatter()", {
  aa <- df2 %>% scatter %>% pick(family)
  bb <- df2 %>% scatter %>% pick(family, species)

  expect_is(aa, "taxa")
  expect_equal(length(aa), NROW(df2))
  expect_is(bb, "taxa")
  expect_equal(length(bb), NROW(df2))

  expect_is(aa[[1]], "taxon")
  expect_is(aa[[1]]$binomial, "binomial")
  expect_is(aa[[1]]$binomial$authority, "character")
  expect_is(aa[[1]]$grouping, "grouping")
  expect_is(aa[[1]]$grouping$family, "taxonref")
  expect_is(aa[[1]]$grouping$family$rank, "character")
  expect_equal(length(aa[[1]]), 2)

  expect_true(all(vapply(aa, function(z) length(z$grouping), integer(1)) == 1))
  expect_true(all(vapply(bb, function(z) length(z$grouping), integer(1)) == 2))
})

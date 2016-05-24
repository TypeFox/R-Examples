context("pop")

out <- make_taxon(genus="Poa", epithet="annua", authority="L.",
                  family='Poaceae', clazz='Poales', kingdom='Plantae', variety='annua')

test_that("pop works with taxon class input - single tax group", {
  aa <- out %>% pop(family)
  bb <- out %>% pop(genus)
  cc <- out %>% pop(species)

  expect_is(out, "taxon")
  expect_is(aa, "taxon")
  expect_is(bb, "taxon")
  expect_is(cc, "taxon")

  expect_named(aa$grouping, c('kingdom', 'clazz', 'genus', 'species', 'variety'))
  expect_named(bb$grouping, c('kingdom', 'clazz', 'family', 'species', 'variety'))
  expect_named(cc$grouping, c('kingdom', 'clazz', 'family', 'genus', 'variety'))
})

test_that("pop works with taxon class input - many tax groups", {
  dd <- out %>% pop(genus, species)
  ee <- out %>% pop(family, genus, species)
  ff <- out %>% pop(family, genus, kingdom)

  expect_is(dd, "taxon")
  expect_is(ee, "taxon")
  expect_is(ff, "taxon")

  expect_named(dd$grouping, c('kingdom', 'clazz', 'family', 'variety'))
  expect_named(ee$grouping, c('kingdom', 'clazz', 'variety'))
  expect_named(ff$grouping, c('clazz', 'species', 'variety'))
})


df <- data.frame(class=c('Magnoliopsida','Magnoliopsida','Magnoliopsida',
                         'Magnoliopsida','Magnoliopsida','Magnoliopsida'),
         order=c('Asterales','Asterales','Fagales','Poales','Poales','Poales'),
         family=c('Asteraceae','Asteraceae','Fagaceae','Poaceae','Poaceae','Poaceae'),
         genus=c('Helianthus','Helianthus','Quercus','Poa','Festuca','Holodiscus'),
         stringsAsFactors = FALSE)
df2 <- taxon_df(df)

test_that("pop works with taxonomic DF's - single tax group", {
  expect_is(df2, "taxondf")
  expect_named(df2, c("class", "order", "family", "genus"))

  aa <- df2 %>% pop(order)
  bb <- df2 %>% pop(family)
  cc <- df2 %>% pop(genus)

  expect_is(aa, "taxondf")
  expect_is(bb, "taxondf")
  expect_is(cc, "taxondf")

  expect_named(aa, c("class", "family", "genus"))
  expect_named(bb, c("class", "order", "genus"))
  expect_named(cc, c("class", "order", "family"))
})

test_that("pop works with taxonomic DF's - many tax groups", {
  dd <- df2 %>% pop(order, family)
  ee <- df2 %>% pop(order, genus)

  expect_is(dd, "taxondf")
  expect_is(ee, "taxondf")

  expect_named(dd, c("class", "genus"))
  expect_named(ee, c("class", "family"))
})


test_that("pop works with taxa objects", {
  x <- df2 %>% scatter
  aa <- x %>% pop(family)
  bb <- x %>% pop(family, genus)
  cc <- x %>% pop(family, genus, species)

  expect_is(x, "taxa")

  expect_is(aa, "taxa")
  expect_is(bb, "taxa")
  expect_is(cc, "taxa")

  expect_null(names(aa))

  expect_is(aa[[1]], "taxon")
  expect_is(bb[[3]], "taxon")
  expect_is(cc[[5]], "taxon")

  expect_named(aa[[1]]$grouping, c("clazz", "order", "genus", "species"))
  expect_named(bb[[1]]$grouping, c("clazz", "order", "species"))
  expect_named(cc[[1]]$grouping, c("clazz", "order"))
})

test_that("pop check for class works", {
  expect_error(pop(), "no applicable method")
})

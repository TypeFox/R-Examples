context("sort")

df <- data.frame(class=c('Magnoliopsida','Magnoliopsida','Magnoliopsida',
                         'Magnoliopsida','Magnoliopsida','Magnoliopsida'),
         order=c('Asterales','Asterales','Fagales','Poales','Poales','Poales'),
         family=c('Asteraceae','Asteraceae','Fagaceae','Poaceae','Poaceae','Poaceae'),
         genus=c('Helianthus','Helianthus','Quercus','Poa','Festuca','Holodiscus'),
         stringsAsFactors = FALSE)
df2 <- taxon_df(df)

aa <- df2 %>% sort(desc(order))
bb <- df2 %>% sort(desc(family))
cc <- df2 %>% sort(genus)
dd <- df2 %>% sort(desc(genus))

test_that("input taxondf object is correct", {
  expect_is(df2, "taxondf")
  expect_equal(NROW(df2), 6)
})

test_that("sort gives back right objects", {
  expect_is(aa, "data.frame")
  expect_is(aa, "taxondf")
  expect_is(bb, "data.frame")
  expect_is(bb, "taxondf")
  expect_is(cc, "data.frame")
  expect_is(cc, "taxondf")
  expect_is(dd, "data.frame")
  expect_is(dd, "taxondf")
})

test_that("sort gives back right sort order", {
  ## aa
  expect_equal(df2$order[1], "Asterales")
  expect_equal(aa$order[1], "Poales")

  ## bb
  expect_equal(df2$family[1], "Asteraceae")
  expect_equal(bb$family[1], "Poaceae")

  ## cc
  expect_equal(df2$genus[1], "Helianthus")
  expect_equal(cc$genus[1], "Festuca")

  ## dd
  expect_equal(df2$genus[1], "Helianthus")
  expect_equal(dd$genus[1], "Quercus")
})

test_that("sort check for class works", {
  expect_error(sort(), "argument \"x\" is missing")
  expect_error(sort(), "argument \"x\" is missing")
})

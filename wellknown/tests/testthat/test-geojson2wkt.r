context("geojson2wkt")

test_that("convert point works", {
  point <- list('type' = 'Point', 'coordinates' = c(116.4, 45.2))
  a <- geojson2wkt(point)
  expect_is(a, "character")
  expect_match(a, "POINT")
  expect_equal(a, "POINT (116.4000000000000057  45.2000000000000028)")
})

test_that("convert multipoint works", {
  mp <- list(type = 'MultiPoint', coordinates=list(c(100.0, 3.101), c(101.0, 2.1), c(3.14, 2.18)))
  b <- geojson2wkt(mp)
  expect_is(b, "character")
  expect_match(b, "POINT")
  expect_match(b, "MULTIPOINT")
  expect_equal(b, "MULTIPOINT ((100.0000000000000000 3.1010000000000000), (101.0000000000000000 2.1000000000000001), (3.1400000000000001 2.1800000000000002))")
})

test_that("convert linestring works", {
  st <- list(type = 'LineString',
              coordinates = list(c(0.0, 0.0, 10.0), c(2.0, 1.0, 20.0),
                                c(4.0, 2.0, 30.0), c(5.0, 4.0, 40.0)))
  c <- geojson2wkt(st)
  expect_is(c, "character")
  expect_match(c, "LINESTRING")
  expect_equal(geojson2wkt(st, 0), "LINESTRING (0 0 10, 2 1 20, 4 2 30, 5 4 40)")
  expect_equal(length(strsplit(geojson2wkt(st, 1), "\\.")[[1]]), 13)
})

test_that("convert multilinestring works", {
  multist <- list(type = 'MultiLineString',
        coordinates = list(
          list(c(0.0, -1.0), c(-2.0, -3.0), c(-4.0, -5.0)),
          list(c(1.66, -31023.5), c(10000.9999, 3.0), c(100.9, 1.1), c(0.0, 0.0))
        ))
  d <- geojson2wkt(multist)
  expect_is(d, "character")
  expect_match(d, "LINESTRING")
  expect_match(d, "MULTILINESTRING")
  expect_equal(geojson2wkt(multist, 0), "MULTILINESTRING ((0 -1, -2 -3, -4 -5), (1.66 -31023.50, 10001 3, 100.9 1.1, 0 0))")
})

test_that("convert polygon works", {
  poly <- list(type = 'Polygon',
        coordinates=list(
          list(c(100.001, 0.001), c(101.12345, 0.001), c(101.001, 1.001), c(100.001, 0.001)),
          list(c(100.201, 0.201), c(100.801, 0.201), c(100.801, 0.801), c(100.201, 0.201))
  ))
  e <- geojson2wkt(poly)
  expect_is(e, "character")
  expect_match(e, "POLYGON")
  expect_equal(geojson2wkt(poly, 0), "POLYGON ((100.001 0.001, 101.1235 0.0010, 101.001 1.001, 100.001 0.001), (100.201 0.201, 100.801 0.201, 100.801 0.801, 100.201 0.201))")
})

test_that("convert multipolygon works", {
  mpoly2 <- list(type = "MultiPolygon",
               coordinates = list(list(list(c(30, 20), c(45, 40), c(10, 40), c(30, 20))),
                              list(list(c(15, 5), c(40, 10), c(10, 20), c(5 ,10), c(15, 5))))
  )
  f <- geojson2wkt(mpoly2, fmt=0)
  expect_is(f, "character")
  expect_match(f, "POLYGON")
  expect_match(f, "MULTIPOLYGON")
  expect_equal(f, "MULTIPOLYGON (((30 20, 45 40, 10 40, 30 20)), ((15 5, 40 10, 10 20, 5 10, 15 5)))")
})

test_that("convert geometrycollection works", {
  gmcoll <- list(type = 'GeometryCollection',
                 geometries = list(
                   list(type = "Point", coordinates = list(0.0, 1.0)),
                   list(type = 'LineString', coordinates = list(c(-100.0, 0.0), c(-101.0, -1.0))),
                   list(type = 'MultiPoint',
                        'coordinates' = list(c(100.0, 3.101), c(101.0, 2.1), c(3.14, 2.18))
                   )
                 )
  )
  g <- geojson2wkt(gmcoll, fmt=0)
  expect_is(g, "character")
  expect_match(g, "GEOMETRYCOLLECTION")
  expect_equal(g, "GEOMETRYCOLLECTION (POINT (0 1), LINESTRING (-100 0, -101 -1), MULTIPOINT ((100.000 3.101), (101.0 2.1), (3.14 2.18)))")
})

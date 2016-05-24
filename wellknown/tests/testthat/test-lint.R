context("lint")

# points --------------------
test_that("lint works for valid WKT strings - points", {
  # good
  expect_true(lint("POINT (1 2)"))
  expect_true(lint("POINT (1 2 3)"))
  expect_true(lint("POINT (1 2 3 4)"))
  expect_true(lint("POINT EMPTY"))
})

test_that("lint works for invalid WKT strings - points", {
  # bad
  expect_false(lint("POINT (1 2 3 4 5)"))
  expect_false(lint("POINT (1 a)"))
  expect_false(lint("POINT (1 )"))
  expect_false(lint("POINT (1)"))
  expect_false(lint("POINT ()"))
  expect_false(lint("POINT "))
  expect_false(lint("POINT empty"))
  expect_false(lint("POINT emp"))
  expect_false(lint("POINT12"))
  expect_false(lint("point (1 2)"))
})


# linestring --------------------
test_that("lint works for valid WKT strings - linestring", {
  # good
  expect_true(lint("LINESTRING EMPTY"))
  expect_true(lint("LINESTRING (100 0, 101 1)"))
  expect_true(lint("LINESTRING (100 0, 101 1, 4 6)"))
  expect_true(lint("LINESTRING (100 0,101 1,4 6)"))
  expect_true(lint("LINESTRING(100 0,101 1,4 6)"))
  expect_true(lint("LINESTRING (100 0, 101 1, 4 6, 4 5, 23434343 45454545)"))
  expect_true(lint("LINESTRING (100 3 4)"))
  expect_true(lint("LINESTRING (100 3 4.454545)"))
  expect_true(lint("LINESTRING (100 3.23434343 0)"))
  expect_true(lint("LINESTRING(7.127130925655365 44.872856971822685,50.545099675655365 42.5869135116598,31.560724675655365 33.02876067050816,42.986505925655365 24.408955129015624,24.002130925655365 22.311272019147477,34.724787175655365 14.308903161964137,15.037287175655365 14.649299083976166,19.431818425655365 6.535627729378449,3.084162175655365 13.113523650016484,4.666193425655365 0.7490445384366899,-7.638494074344635 11.223467666031517,-14.845525324344635 -2.5898954605716416,-20.470525324344635 11.912306151187774,-27.150212824344635 -1.1844435030433105,-29.435369074344635 18.68976902906233,-16.251775324344635 17.687796432467547,-22.579900324344635 26.472359363941862,-10.978337824344635 25.999343543885978,-19.240056574344635 31.3928609826332,-3.595525324344635 29.72795526109783,-10.978337824344635 36.21043323531069,3.963068425655365 33.76254102463172,-1.661931574344635 39.81006899259359)"))
})

test_that("lint works for valid WKT strings - linestring", {
  # bad
  expect_false(lint("LINESTRING (100 3 4 5)"))
  expect_false(lint("LINESTRING (100 3 4 5 7)"))
  expect_false(lint("LINESTRING (100 adf)"))
  expect_false(lint("LINESTRING (100)"))
  expect_false(lint("LINESTRING ()"))
  expect_false(lint("LINESTRING "))
  expect_false(lint("LINESTRING ("))
  expect_false(lint("LINESTRING )"))
  expect_false(lint("LINESTRING"))
  expect_false(lint("LINESTRING (100 4, 1)"))
  expect_false(lint("LINESTRING (100 4, 1 ad)"))
  expect_false(lint("LINESTRING (100 4, 1, 1)"))
  expect_false(lint("linestring (1 2)"))
  expect_false(lint("Linestring (1 2)"))
  expect_false(lint("LineString (1 2)"))
})

# polygon --------------------
test_that("lint works for valid WKT strings - polygon", {
  # good
  expect_true(lint("POLYGON EMPTY"))
  expect_true(lint("POLYGON ((1 2, 3 4, 0 5, 1 2))"))
  expect_true(lint("POLYGON ((1 2, 3 4, 0 5, 6 7, 1 2))"))
  expect_true(lint("POLYGON((1 2, 3 4, 0 5, 1 2))"))
  expect_true(lint("POLYGON ((1.23 2.23, 3.45 4.23, 0.5 5.12, 1.1 2.78))"))
  expect_true(lint("POLYGON((1.2784092128276825 30.31684905459481,12.792081087827682 31.373251008758597,12.704190462827682 18.167678453227037,-0.1278407871723175 27.392164134909823,1.2784092128276825 30.31684905459481))"))
})

test_that("lint works for valid WKT strings - polygon", {
  # bad
  expect_false(lint("POLYGON (100 3 4)"))
  expect_false(lint("POLYGON ((1 2, 3 4, 0 5, 1 a))"))
  expect_false(lint("POLYGON ((1 2, 3 4, 0 a.3, 1 3))"))
  expect_false(lint("POLYGON ((1 2, 3 4, 0 5, 1 a)"))
  expect_false(lint("POLYGON ((1 2, 3 4, 0 5, 1 6)))"))
  expect_false(lint("POLYGON ((1 2, 3 4, 05, 1 5))"))
  expect_false(lint("POLYGON (100)"))
  expect_false(lint("POLYGON ()"))
  expect_false(lint("POLYGON "))
  expect_false(lint("POLYGON ("))
  expect_false(lint("POLYGON )"))
  expect_false(lint("POLYGON"))
  expect_false(lint("POLYGON (100 4, 1)"))
  expect_false(lint("POLYGON (100 4, 1 ad)"))
  expect_false(lint("POLYGON (100 4, 1, 1)"))
  expect_false(lint("polygon (1 2)"))
  expect_false(lint("Polygon (1 2)"))
})

# multipolygon --------------------
test_that("lint works for valid WKT strings - multipolygon", {
  # good
  expect_true(lint("MULTIPOLYGON EMPTY"))
  expect_true(lint("MULTIPOLYGON (((30 20, 45 40, 10 40, 30 20)))"))
  expect_true(lint("MULTIPOLYGON (((30 20, 45 40, 10 40, 30 20)), ((15 5, 40 10, 10 20, 5 10, 15 5)))"))
  expect_true(lint("MULTIPOLYGON (((40 40, 20 45, 45 30, 40 40)), ((20 35, 45 20, 30 5, 10 10, 10 30, 20 35)))"))
})

test_that("lint works for valid WKT strings - multipolygon", {
  # should be good
  expect_false(lint("MULTIPOLYGON (((40 40, 20 45, 45 30, 40 40)), ((20 35, 45 20, 30 5, 10 10, 10 30, 20 35), (30 20, 20 25, 20 15, 30 20)))"))

  # bad
  expect_false(lint("MULTIPOLYGON (((30 20, 45 40, 10 40, 30)))"))
  expect_false(lint("MULTIPOLYGON (((30 20, 45 40, 10 40,)))"))
  expect_false(lint("MULTIPOLYGON (((30 20, 45 40, 10 40, 30 a)))"))
  expect_false(lint("MULTIPOLYGON (((30 20, 45 40, 10.adsfaf 40, 30 20)))"))
  expect_false(lint("MULTIPOLYGON (((30 20, 45 40, 10 40, 30 )"))
  expect_false(lint("MULTIPOLYGON (((30 20, 45 40, 10 40, 3020)))"))
  expect_false(lint("MULTIPOLYGON (((30)))"))
  expect_false(lint("MULTIPOLYGON ((()))"))
  expect_false(lint("MULTIPOLYGON "))
  expect_false(lint("MULTIPOLYGON ("))
  expect_false(lint("MULTIPOLYGON )"))
  expect_false(lint("MULTIPOLYGON"))
  expect_false(lint("MULTIPOLYGON (100 4, 1)"))
  expect_false(lint("MULTIPOLYGON (100 4, 1 ad)"))
  expect_false(lint("MULTIPOLYGON (100 4, 1, 1)"))
  expect_false(lint("MULTIpolygon (((30 20, 45 40, 10 40, 30 20)))"))
  expect_false(lint("MULTIPolygon (((30 20, 45 40, 10 40, 30 20)))"))
})


# triangle --------------------
test_that("lint works for valid WKT strings - triangle", {
  # good
  expect_true(lint("TRIANGLE EMPTY"))
  expect_true(lint("TRIANGLE ((0 0, 0 1, 1 1, 0 0))"))
  expect_true(lint("TRIANGLE ((0 0, 0 1, 1 1))"))
  expect_true(lint("TRIANGLE ((0 0, 0 1))"))
})

test_that("lint works for valid WKT strings - triangle", {
  # bad
  expect_false(lint("TRIANGLE ((30 20, 45 40, 10 40, 30))"))
  expect_false(lint("TRIANGLE ((30 20, 45 40, 10 40,))"))
  expect_false(lint("TRIANGLE ((30 20, 45 40, 10 40, 30 a))"))
  expect_false(lint("TRIANGLE ((30 20, 45 40, 10.adsfaf 40, 30 20))"))
  expect_false(lint("TRIANGLE ((30 20, 45 40, 10 40, 30 )"))
  expect_false(lint("TRIANGLE ((30 20, 45 40, 10 40, 3020))"))
  expect_false(lint("TRIANGLE ((30))"))
  expect_false(lint("TRIANGLE (())"))
  expect_false(lint("TRIANGLE "))
  expect_false(lint("TRIANGLE ("))
  expect_false(lint("TRIANGLE )"))
  expect_false(lint("TRIANGLE"))
  expect_false(lint("TRIANGLE (100 4, 1)"))
  expect_false(lint("TRIANGLE (100 4, 1 ad)"))
  expect_false(lint("TRIANGLE (100 4, 1, 1)"))
  expect_false(lint("TRIANGLE (((30 20, 45 40, 10 40, 30 20)))"))
  expect_false(lint("TRIANGLE (((30 20, 45 40, 10 40, 30 20)))"))
})



# circularstring --------------------
test_that("lint works for valid WKT strings - circularstring", {
  # good
  expect_true(lint("CIRCULARSTRING EMPTY"))
  expect_true(lint("CIRCULARSTRING(1 5, 6 2, 7 3)"))
})

test_that("lint works for valid WKT strings - circularstring", {
  # bad
  expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40, 30)"))
  expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40,)"))
  expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40, 30 a)"))
  expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10.adsfaf 40, 30 20)"))
  expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40, 30 )"))
  expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40, 3020)"))
  expect_false(lint("CIRCULARSTRING (30)"))
  expect_false(lint("CIRCULARSTRING ()"))
  expect_false(lint("CIRCULARSTRING "))
  expect_false(lint("CIRCULARSTRING ("))
  expect_false(lint("CIRCULARSTRING )"))
  expect_false(lint("CIRCULARSTRING"))
  expect_false(lint("CIRCULARSTRING (100 4, 1)"))
  expect_false(lint("CIRCULARSTRING (100 4, 1 ad)"))
  expect_false(lint("CIRCULARSTRING (100 4, 1, 1)"))
  ## fixme - these shouldn't pass linting
  # expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40, 30 20)"))
  # expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40, 30 20)"))
})


# circularstring --------------------
test_that("lint works for valid WKT strings - circularstring", {
  # good
  expect_true(lint("CIRCULARSTRING EMPTY"))
  expect_true(lint("CIRCULARSTRING(1 5, 6 2, 7 3)"))
})

test_that("lint works for valid WKT strings - circularstring", {
  # bad
  expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40, 30)"))
  expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40,)"))
  expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40, 30 a)"))
  expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10.adsfaf 40, 30 20)"))
  expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40, 30 )"))
  expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40, 3020)"))
  expect_false(lint("CIRCULARSTRING (30)"))
  expect_false(lint("CIRCULARSTRING ()"))
  expect_false(lint("CIRCULARSTRING "))
  expect_false(lint("CIRCULARSTRING ("))
  expect_false(lint("CIRCULARSTRING )"))
  expect_false(lint("CIRCULARSTRING"))
  expect_false(lint("CIRCULARSTRING (100 4, 1)"))
  expect_false(lint("CIRCULARSTRING (100 4, 1 ad)"))
  expect_false(lint("CIRCULARSTRING (100 4, 1, 1)"))
  ## fixme - these shouldn't pass linting
  # expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40, 30 20)"))
  # expect_false(lint("CIRCULARSTRING (30 20, 45 40, 10 40, 30 20)"))
})


#' lint("CIRCULARSTRING (1 5, 6 2, 7 3)")
#' lint("CIRCULARSTRING (1 5, 6 2, 7 3, 5 6, 4 3)")
#' lint('COMPOUNDCURVE (CIRCULARSTRING (1 0, 0 1, -1 0), (-1 0, 2 0))')

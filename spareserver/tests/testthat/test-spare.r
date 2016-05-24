
context("Redundant services")

my_services <- get("spare_services", envir = asNamespace("spareserver"))
on.exit(assign("spare_services", my_services,
               envir = asNamespace("spareserver")), add = TRUE)


test_that("Can create services", {
  clean_services()
  s1 <- server("http://server1.com", priority = 10)
  s2 <- server("http://server2.com", priority = 1)
  add_service("service1", s1, s2)
  expect_equal(names(services()), "service1")
})

test_that("Can add servers", {
  clean_services()

})

test_that("Can remove services", {
  clean_services()

})

context("Redundant queries")

test_that("First server is used if available", {
  clean_services()

  s1 <- server("http://google.com:80/", priority = 10)
  s2 <- server("http://192.0.2.1/foobar/", priority = 5)
  add_service("test", s1, s2)

  q <- spare_q("test", "", httr::GET)
  expect_equal(httr::status_code(q), 200)

})

test_that("Second server is used if first is not available", {
  clean_services()

  s1 <- server("http://google.com:80/", priority = 10)
  s2 <- server("http://192.0.2.1/foobar/", priority = 20)
  add_service("test", s1, s2)

  q <- spare_q("test", "", httr::GET)
  expect_equal(httr::status_code(q), 200)

})

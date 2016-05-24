context("core functionality")

test_that("bogons work", {
  expect_that(ipv4_bogons(force=TRUE), is_a("character"))
  expect_that(ipv6_bogons(force=TRUE), is_a("character"))
})

test_that("bulk_origin works", {
  expect_that(bulk_origin(c("68.22.187.5", "207.229.165.18", "198.6.1.65"),
                          timeout=1),
              is_a("data.frame"))
})

test_that("bulk_peer works", {
  expect_that(bulk_peer(c("68.22.187.5", "207.229.165.18", "198.6.1.65"),
                        timeout=1),
              is_a("data.frame"))
})

test_that("bulk_origin_asn works", {
  expect_that(bulk_origin_asn(c(22822, 1273, 2381, 2603, 2914, 3257, 3356, 11164,
                                174, 286, 1299, 2914, 3257, 3356, 3549, 22822),
                              timeout=1),
              is_a("data.frame"))
})

test_that("malware_hash works", {
  expect_that(malware_hash(c("1250ac278944a0737707cf40a0fbecd4b5a17c9d",
                             "7697561ccbbdd1661c25c86762117613",
                             "cbed16069043a0bf3c92fff9a99cccdc",
                             "e6dc4f4d5061299bc5e76f5cd8d16610",
                             "e1112134b6dcc8bed54e0e34d8ac272795e73d74"),
                           timeout=1),
              is_a("data.frame"))
})

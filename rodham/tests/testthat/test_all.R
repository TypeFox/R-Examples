library(rodham)

context("All tests")

test_that("test all", {
  # search
  emails <- search_emails()
  expect_equal(nrow(emails), 29444)
  expect_equal(ncol(emails), 9)
  expect_equal(names(emails), c("docID", "docDate", "to", "from",
                                "originalTo", "originalFrom", "subject",
                                "interesting", "not_interesting"))
  # search
  emails <- search_emails(internal = FALSE)
  expect_equal(nrow(emails), 29444)
  expect_equal(ncol(emails), 9)
  expect_equal(names(emails), c("docID", "docDate", "to", "from",
                                "originalTo", "originalFrom", "subject",
                                "interesting", "not_interesting"))
  # edges
  expect_error(edges_emails())
  expect_equal(nrow(edges_emails(emails)), 795)
  expect_equal(ncol(edges_emails(emails)), 3)
  expect_equal(ncol(edges_emails(emails, "subject")), 4)
  expect_equal(names(edges_emails(emails, "subject", "docDate")),
               c("from", "to", "subject", "docDate", "freq"))
  # errors get
  expect_error(get_emails())
})

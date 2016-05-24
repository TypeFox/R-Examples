context('read_blast_table')

b <- read_blast_table(
  system.file("testdata", "blast_table.txt", package="qiimer"))

test_that("Query IDs are parsed correctly", {
  expect_equal(b$query_id, rep("seq01", 6))
})

test_that("Subject IDs are parsed correctly", {
  expect_equal(b$subject_id[1], c("CF6VF"))
})

test_that("Subsequent columns are numeric or integer", {
  expect_is(b$pct_identity, "numeric")
  expect_is(b$alignment_len, "integer")
  expect_is(b$mismatches, "integer")
  expect_is(b$gap_openings, "integer")
  expect_is(b$query_start, "integer")
  expect_is(b$query_end, "integer")
  expect_is(b$subject_start, "integer")
  expect_is(b$subject_end, "integer")
  expect_is(b$e_value, "numeric")
  expect_is(b$bit_score, "numeric")
})

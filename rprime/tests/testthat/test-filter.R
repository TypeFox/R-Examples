# Test files
good_file <- "data/MINP_001L00XS1.txt"

context("Filtering functions")

test_that("filter in/out", {
  eprime_log <- read_eprime(good_file)
  chunked <- FrameList(eprime_log)

  just_header <- filter_in(chunked, "Running", "Header")
  header_frame <- just_header[[1]]
  header_df <- to_data_frame(just_header)

  # There should just be one frame when we filter to just the header
  expect_equal(length(just_header), 1)
  expect_equal(nrow(header_df), 1)
  expect_equal(header_frame$Running, "Header")

  # Filter down to just the Familiarization trials
  level_2 <- filter_in(chunked, "Eprime.Level", 2)
  just_practice <- filter_out(level_2, "Running", "Test")
  practice_df <- to_data_frame(just_practice)
  level_2_df <- to_data_frame(level_2)

  # There are six familiarization trials
  expect_equal(length(just_practice), 6)
  expect_equal(nrow(practice_df), 6)
  expect_equal(unique(practice_df$Running), "Familiarization")

  # Filtering with multiple values
  keys_peas <- filter_in(level_2, "CorrectResponse", c("keys", "peas"))
  keys_peas_df <- to_data_frame(keys_peas)
  # Test against subset
  keys_peas_rows <- nrow(subset(level_2_df, CorrectResponse %in% c("keys", "peas")))
  expect_equal(nrow(keys_peas_df), keys_peas_rows)
})

test_that("keep/drop levels", {
  eprime_log <- read_eprime(good_file)
  chunked <- FrameList(eprime_log)
  df_full <- to_data_frame(chunked)

  just_level_1 <- keep_levels(chunked, 1)
  df_level_1 <- to_data_frame(just_level_1)

  just_level_2 <- keep_levels(chunked, 2)
  df_level_2 <- to_data_frame(just_level_2)

  # Levels 1 and 2 partition the trials
  expect_equal(nrow(df_level_1) + nrow(df_level_2), nrow(df_full))

  # Dropping level 2 and keeping level 1 should return the same frames. Test by
  # comparing the dfs (need to omit NAs first to allow comparison.)
  not_level_2 <- drop_levels(chunked, 2)
  df_not_level_2 <- to_data_frame(not_level_2)

  remove_na_columns <- function(df) Filter(function(x) !any(is.na(x)), df)
  df_level_1 <- remove_na_columns(df_level_1)
  df_not_level_2 <- remove_na_columns(df_not_level_2)
  expect_equal(df_not_level_2, df_level_1)

})

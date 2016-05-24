## ------------------------------------------------------------------------
# Load the package
library(bpa)  # pipe operator %>% included for convenience

# Load the data
data(messy, package = "bpa")
head(messy)

## ------------------------------------------------------------------------
get_pattern("(123)  456-789")
get_pattern("(123)  456-789", show_ws = FALSE)  # same as ws_char = " "
get_pattern("(123)  456-789", ws_char = "<ws>")

## ------------------------------------------------------------------------
messy$Date %>%
  get_pattern %>%  # extract patterns
  table %>%        # tabulate frequencies
  as.data.frame    # display as a data frame

## ------------------------------------------------------------------------
messy$Date %>%
  unique %>%    # extract unique values
  head(50)      # look at first 50 observations

## ------------------------------------------------------------------------
# Standardize the entire data set (returns a data frame)
messy %>%
  basic_pattern_analysis %>%  # note: you can also use bpa for short
  head(10)                    # only look at first 10 observations

## ------------------------------------------------------------------------
# Only return unique patterns (returns a list)
bpa(messy, unique_only = TRUE)

## ------------------------------------------------------------------------
# Extract Gender values matching the pattern "Aaaa"
match_pattern(messy$Gender, pattern = "Aaaa", unique_only = TRUE)


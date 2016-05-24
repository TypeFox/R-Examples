patterns <- list(
  bracket       = "\\*{3} (.*) (Start|End) \\*{3}",
  bracket_start = "\\*{3} (.*) Start \\*{3}",
  bracket_end   = "\\*{3} (.*) End \\*{3}",
  header_start  = "^\\*{3} Header Start \\*{3}$",
  header_end    = "^\\*{3} Header End \\*{3}$",
  footer_start  = "^\\*{3} LogFrame Start \\*{3}$",
  footer_end    = "^\\*{3} LogFrame End \\*{3}$",
  row           = ".+: .*"
)

# Names of the metadata fields.
rprime_cols <- list(
  level      = "Eprime.Level",
  level_name = "Eprime.LevelName",
  basename   = "Eprime.Basename",
  frame      = "Eprime.FrameNumber",
  procedure  = "Procedure",
  running    = "Running"
)

# Make a list of the metadata fields initialized to NAs
default_metadata <- structure(
  .Data = as.list(rep(NA, length(rprime_cols))),
  names = unlist(unname(rprime_cols)))

str_which <- function(string, pattern) which(str_detect(string, pattern))
new_line <- function(key, value) sprintf("%s: %s", key, value)
is_bracket <- function(xs) str_detect(xs, patterns$bracket)
is_header <- function(xs) any(str_detect(xs, patterns$header_start))
is_row <- function(xs) any(str_detect(xs, patterns$row))
has_header <- function(xs) any(sapply(xs, is_header))


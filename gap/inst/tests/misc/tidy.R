options(keep.source = FALSE)
source("t.R")
dump(ls(all = TRUE), file = "tt.R")

# cp myfile.R t.R
# R --vanilla < tidy.R
# pico tt.R

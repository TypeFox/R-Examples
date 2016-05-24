library(dplyr)

src <- "https://www.kaggle.com/c/titanic/data/download/test.csv"
lcl <- "data-raw/test.csv"

# if (!file.exists(lcl)) {
#   download.file(src, lcl)
# }

raw <- read.csv(lcl, stringsAsFactors = FALSE)
titanic_test <- raw

save(airlines, file = "data/titanic_test.rda")

library(dplyr)

src <- "https://www.kaggle.com/c/titanic/data/download/train.csv"
lcl <- "data-raw/train.csv"

# if (!file.exists(lcl)) {
#   download.file(src, lcl)
# }

raw <- read.csv(lcl, stringsAsFactors = FALSE)
titanic_train <- raw

save(airlines, file = "data/titanic_train.rda")

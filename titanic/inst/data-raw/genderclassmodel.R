library(dplyr)

src <- "https://www.kaggle.com/c/titanic/data/download/genderclassmodel.csv"
lcl <- "data-raw/genderclassmodel.csv"

# if (!file.exists(lcl)) {
#   download.file(src, lcl)
# }

raw <- read.csv(lcl, stringsAsFactors = FALSE)
titanic_gender_class_model <- raw

save(airlines, file = "data/titanic_gender_class_model.rda")

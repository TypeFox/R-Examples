# docFile() will tell us where the file is on our computer
docFile("OneTrueLove.csv", package = "Lock5withR")
OneTrueLove3 <- read.file(docFile("OneTrueLove.csv", package = "Lock5withR"))
OneTrueLove4 <- read.file("OneTrueLove.csv", package = "Lock5withR")


## ---- echo = FALSE-------------------------------------------------------
pkgVersion <- packageDescription("seroincidence")$Version
pkgDate <- packageDescription("seroincidence")$Date
authorsString <- gsub("^ *|(?<= ) |\n| *$", "", packageDescription("seroincidence")$Authors, perl = TRUE)
authorList <- eval(parse(text = authorsString))
pkgAuthors <- paste(format(authorList, include = c("given", "family", "email", "comment"), braces = list(email = c("<", ">,<br />"), comment = c("", ""))), collapse = "<br /><br />")
pkgMaintainer <- packageDescription("seroincidence")$Maintainer
pkgBaseFileName <- paste("seroincidence", pkgVersion, sep = "_")
pkgUrl <- packageDescription("seroincidence")$URL


x <- scan("/Users/Stoffi/R_files/BEAST/alpinae.xml", what = "c", sep = "\n")

x <- gsub("\t", "", x, fixed = TRUE)
x <- gsub("[<]!--.+--[>]", "", x)
x <- x[x != ""]
x <- x[(grep("</alignment>", x) + 1) : length(x)]

y <- scan("/Users/Stoffi/R_files/BEAST/alpinaeR.xml", what = "c", sep = "\n")
y <- gsub("^[[:space:]]+", "", y)
y <- y[(grep("</alignment>", y) + 1) : length(y)]

maxlen <- 300
z <- vector(length = maxlen)
for (i in 1:maxlen)
z[i] <- identical(x[i], y[i])

id <- which(z == FALSE)[11]
x[id]
y[id]

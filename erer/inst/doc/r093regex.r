# Regular expressions not recognized
two.state <- tolower(paste(state.name[1:2], collapse = "a."))
da <- unlist(strsplit(x = two.state, split = "a.", fixed = TRUE))
two.state

# Regular expressions recognized
db <- unlist(strsplit(x = two.state, split = "a.", fixed = FALSE))
dc <- unlist(strsplit(x = two.state, split = "(a.){3}", fixed = FALSE))
da; db; dc

# Partial string matching
p1 <- match(x = "sum", table = c("mean", "median", "sum"))
p2 <- "sum" %in% c("mean", "median", "sum")
p3 <- pmatch(x = "med", table = c("mean", "median", "sum"))
p4 <- pmatch(x = "med", table = c("mean", "median", "medal"))
p5 <- pmatch(x = c("", "ab", "ab"), table = c("abc", "ab"), 
  duplicates.ok = FALSE)
p6 <- charmatch(x = c("a", "b", "c"), table = "c", nomatch = 0)
p1; p2; p3; p4; p5; p6

ma <- grep(pattern = "a", x = letters[1:5])
mb <- grep(pattern = "med", x = c("mean", "median", "sum"))
mc <- grep(pattern = "med", x = c("mean", "median", "medal"))
md <- grep(pattern = "[a-z]", x = letters[1:5])
sa <- sub(pattern = "a", replacement = "zzz", x = letters[1:5])
sb <- sub(pattern = "[ac]", replacement = "zzz", x = letters[1:5])
ma; mb; mc; md; sa; sb

# Formatting numbers as a string
gg1 <- c(5, 0, 18.34)
gg2 <- round(x = gg1, digits = 3)
gg3 <- format(x = gg1, trim = TRUE)
gg4 <- format(x = gg1, nsmall = 3, trim = FALSE)  # a space inserted
gg5 <- sprintf(fmt = "%.3f", gg1)  # no space before 5.000
gg1; gg2; gg3; gg4; gg5
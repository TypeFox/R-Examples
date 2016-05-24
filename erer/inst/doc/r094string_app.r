# A1. Create a string without recycling
dog2 <- rep(x = c("S", "J", "K", "A", "C"), each = 8)
pig2 <- rep(x = rep(x = c("V", "Q"), each = 4), times = 5)
cow2 <- rep(x = 1:4, times = 10)
all2 <- paste(dog2, pig2, cow2, sep = "")

# A2. Create a string with recycling
dog3 <- rep(x = c("S", "J", "K", "A", "C"), each = 8)
pig3 <- rep(x = c("V", "Q"), each = 4)
cow3 <- 1:4
all3 <- paste(dog3, pig3, cow3, sep = "")

# A3. Create a string with a concise command line
all4 <- paste(rep(x = c("S", "J", "K", "A", "C"), each = 8),
  rep(x = c("V", "Q"), each = 4), 1:4, sep = "")
noquote(all4)
  
# B1. Extract string elements by literal characters: substr(), pmatch()
begi <- substr(x = all4, start = 1, stop = 1)
loca <- pmatch(x = begi, table = "A", nomatch = 0, duplicates.ok = TRUE)
loc2 <- as.logical(loca); nam2 <- all4[loc2]
loc2; nam2

# B2. Extract string elements by regular expression: grep()
loc3 <- grep(x = all4, value = FALSE, pattern = "^A")
nam3 <- grep(x = all4, value = TRUE,  pattern = "^A")
nam3 <- grep(x = all4, value = TRUE,  pattern = "A..")
nam3 <- grep(x = all4, value = TRUE,  pattern = "A[VQ][1234]")
nam3 <- grep(x = all4, value = TRUE,  pattern = "A[VQ][1-4]")
nam3 <- grep(x = all4, value = TRUE,  pattern = "A[[:alpha:]][[:digit:]]")
loc3; nam3 
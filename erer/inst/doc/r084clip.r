# Manual data inputs in R session directly
v1 <- c(3, 5.1, 7)
v2 <- seq(from = 30, to = 40, by = 4)
v3 <- c(letters[1:2], LETTERS[4:5])
v4 <- c("red", "green")
v5 <- data.frame(state = c("FL", "CA"), size = c(25, 28))
v6 <- list(country = c("USA", "Mexico"), price = 1:8)
v1; v2; v3; v4; v5; v6

# Copying data from Excel to R interactively
# The following block is for illustration only.
# You need a sample Excel document and then use "Ctrl + c" to make a copy.
y <- readClipboard()  # A single character column
y <- as.numeric(readClipboard())  # A single numeric column
y <- read.table(file = "clipboard", header = TRUE, sep = "\t")
y <- read.delim(file = "clipboard", header = TRUE)
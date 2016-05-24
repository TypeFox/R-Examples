# A. Conditional if
# A1. Nesting if statements; negation on conditional object
a.default <- FALSE; b.default <- FALSE
tree <- 10
if (a.default) {
  leaf <- tree * 2
} else if (!b.default) {
  leaf <- tree * 30
}
leaf
  
# A2. A single logical value required
a <- 1:5; b <- 11:15; (one <- a > 2 && b < 14); (two <- a > 2 &  b < 14)
if (a > 2 && b < 14) {
  pig <- 50 
} else {
  pig <- 100
}
pig

# Multiple logical values: works but generates warning
if (a > 2 & b < 14) {dog <- 50} else {dog <- 100}; dog

# A3. Used in an assignment
x <- 87
y1 <- if (x != 2) {x + 10} else {x + 50}        # concise
if (x != 2) {y2 <- x + 10} else {y2 <- x + 50}  # clear
c(x, y1, y2)
               
# B. ifelse() function ----------------------------------------------------
# B1. length(test) = length(yes) = length(no) = 1
y3 <- ifelse(test = x != 2, yes = x + 10, no = x + 50)
c(x, y3)

# B2. length(test) > length(yes), length(no): recycle yes/no
w <- 1:6; mm <- 10; nn <- 31:33
y4 <- ifelse(test = w < 3, yes = mm + 1, no = nn + 50)
w; y4

# B3. length(test) < length(yes), length(no): truncate yes/no 
# The output length is always determined by "test".
z <- 1:10
y5 <- ifelse(test = x != 2, yes = z + 10, no = z + 50)
c(x, y5)

# B4. Nesting ifelse()
pvalue <- c(0.31, 0.04, 0.001, 0.08, 0.02)
svalue <- ifelse(test = pvalue >  0.10, yes = "", 
     no = ifelse(test = pvalue >= 0.05, yes = "*", 
     no = ifelse(test = pvalue >= 0.01, yes = "**", no = "***")))
(com <- data.frame(pvalue, svalue))

# C. switch() function ----------------------------------------------------
# C1. EXPR = a character string; return some results
ff <- c(10, 20, 30, 40)
choice1 <- "mmm"; choice2 <- "max"; choice3 <- "fix"; choice4 <- "other"
w1 <- switch(EXPR = choice1, mmm = mean(ff), max = max(ff), fix = 7, 99)
w2 <- switch(EXPR = choice2, mmm = mean(ff), max = max(ff), fix = 7, 99)
w3 <- switch(EXPR = choice3, mmm = mean(ff), max = max(ff), fix = 7, 99)
w4 <- switch(EXPR = choice4, mmm = mean(ff), max = max(ff), fix = 7, 99)
w5 <- switch(EXPR = choice4, mmm = mean(ff), max = max(ff))  # no default
c(w1, w2, w3, w4); w5

# C2. EXPR = a number: cannot set default value
sa <- switch(EXPR = 1, num = 1:5, LETTERS[1:3], 30, 40, 50)
sb <- switch(EXPR = 2, num = 1:5, LETTERS[1:3], 30, 40, 50)
sc <- switch(EXPR = 3, num = 1:5, LETTERS[1:3], 30, 40, 50)
sf <- switch(EXPR = 6, num = 1:5, LETTERS[1:3], 30, 40, 50) # empty slot 6
sa; sb; sc; sf 

# C3. EXPR = a character string; return a function
library(erer); data(daIns)
regress <- "lm"  # regress <- "glm"  
reg.model <- switch(EXPR = regress, glm = glm, lm = lm)
ra <- reg.model(formula = Y ~ Injury + HuntYrs, data = daIns)
str(reg.model); class(ra)

# A similar function selection by if statement
if (regress == "lm") {
  ra <- lm(formula = Y ~ Injury + HuntYrs, data = daIns)
}
if (regress == "glm") {
  ra <- glm(formula = Y ~ Injury + HuntYrs, data = daIns)
}
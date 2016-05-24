# Load two data frames: daIns, daInsNam
library(erer); data(daIns); data(daInsNam) 

# Task 1: Evaluate and display data frames quickly
dim(daIns); dim(daInsNam)
str(daIns); str(daInsNam) 
attributes(daInsNam); colnames(daIns); rownames(daInsNam) 
unique(daIns$Marital); range(daIns$Inc)
head(daIns); tail(daInsNam); daIns[1:5, 1:10]; fix(daIns); daInsNam 

# Task 2: Rename rows or columns of a data frame
(oldName <- rownames(daInsNam))
daIns.col <- daIns
colnames(daIns.col)[1:3] <- c("y", "injury", "hunt.yrs")
rownames(daIns.col)[1:4] <- paste("respondent", 1:4, sep = "")
head(daIns.col[, 1:5])

# Task 3: Change the mode of a column
daInsNam2 <- daInsNam
daInsNam2[, 1] <- as.character(daInsNam[, 1])
daInsNam2 <- apply(X = daInsNam, MARGIN = 2, FUN = as.character)
str(daInsNam); str(daInsNam2)

# Task 4: Reorder or sort a data frame by row or column
small <- daIns[1:20, 1:13]
small2 <- small[c(11:20, 1:10), c(1:11, 13, 12)]   # Manual reordering
age.sor <- sort(small$Age, decreasing = FALSE)     # sorted content
age.loc <- order(small$Age, decreasing = FALSE)    # location
age.sor2 <- small$Age[age.loc]
identical(age.sor, age.sor2)  # TRUE

sm2 <- small[order( small$Age, small$HuntYrs, decreasing = TRUE), ]
sm3 <- small[order(-small$Age, small$HuntYrs, decreasing = FALSE), ]
rownames(sm3) <- 1:nrow(sm3); sm2; sm3
sm4 <- small[, sort(colnames(small))]
sm5 <- small[, order(colnames(small))]
identical(sm4, sm5)  # TRUE

# Task 5: Extract or replace some elements
# Extract elements
daIns2 <- daIns[c(1:3, 5), c("Gender", "Age", "Race")]
daIns3 <- daIns[daIns$Gender == 0 & daIns$Inc > 100, ]
daIns4 <- daIns[daIns$Gender != 0 & daIns$Age >= 71, ]
daIns3 <- subset(daIns, subset = Gender == 0 & Inc > 100)
daIns4 <- subset(daIns, subset = Gender != 0 & Age >= 71)
head(daIns2); head(daIns3); head(daIns4)

# Replace elements
daIns[3, "Inc"] <- daIns[3, "Inc"] + 5
daIns[daIns$Gender == 0 & daIns$Inc > 100, "Inc"] <- 155
daInsNam[2, 1] <- "injury2"   # cannot change a factor content arbitrarily
daInsNam2[2, 1] <- "injury2"  # yes, change a character vector

# Task 6: Remove some rows or columns from a data frame
loc <- match(x = c("Edu", "FishYrs"), table = colnames(daIns))
daIns3 <- daIns[-c(1:1648), -c(11, 14)]
daIns3 <- daIns[-c(1:1648), -loc]

tiny <- data.frame(daIns[1:5, 1:6], city = c("a", "b", NA, "d", NA))
tiny[2, 3] <- NA
tiny2 <- na.omit(tiny); tiny3 <- na.exclude(tiny)
cc <- complete.cases(tiny); tiny4 <- tiny[cc, ]
str(tiny2); str(tiny3); str(tiny4)

# Task 7: Add or transform rows / columns
tx <- daIns[1:5, 1:3]
tx$abb <- LETTERS[1:nrow(tx)]
tx[, ncol(tx) + 1] <- paste("person", 1:nrow(tx), sep = "")
colnames(tx)[ncol(tx)] <- "person.name"
tx[nrow(tx) + 1, 1:3] <- c(1, 1, 30)
tx[6, 4:5] <- c("F", "person6")
tx; str(tx)

wa <- with(tx, expr = log(HuntYrs))  # return the expression value
wb <- within(tx, expr = log(HuntYrs))  # return the unmodified data
wc <- within(tx, expr = {logHuntYrs <- log(HuntYrs)}) #return modified data
wd <- transform(tx, logHuntYrs = log(HuntYrs)) # similar to within()

# Task 8: Generate a random sample
ny <- tx
set.seed(2); sa <- sample(x = 1:nrow(ny), size = 3,  replace = FALSE)
set.seed(8); sb <- sample(x = 1:nrow(ny), size = 6,  replace = FALSE)
set.seed(5); sc <- sample(x = 1:nrow(ny), size = 20, replace = TRUE)
set.seed(4); sf <- sample(x = 1:nrow(ny), size = 20, 
  replace = FALSE)  # Error: conflict between size and replace
ny.a <- ny[sa, ]; ny.b <- ny[sb, ]; ny.c <- ny[sc, ]
ny.c2 <- ny.c[order(ny.c$person.name), ]

# Task 9: Combine two or more data frames
ak1 <- merge(x = ny.a, y = ny.b, all = FALSE)
ak2 <- cbind(tx, ny.b); ak3 <- rbind(tx, ny.b) 
rownames(ak3) <- 1:nrow(ak3)

# Task 10: Reshape between narrow and wide formats
dog <- daIns[1:3, 1:3]; rownames(dog) <- paste("person", 1:3, sep = ".")
ua <- stack(x = dog[, 1:3])
rownames(ua) <- paste(rownames(ua), rownames(dog), sep = '.'); ua
unstack(ua)

ub <- reshape(data = dog, idvar = "person", ids = row.names(dog),
  times = names(dog), timevar = "characteristic", 
  varying = list(names(dog)), direction = "long") 
reshape(ub, direction = "wide")
uc <- reshape(ub, direction = "wide", new.row.names = unique(ub$person))
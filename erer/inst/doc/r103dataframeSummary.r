# 0. Create a small data set
library(erer); data(daIns)
sam <- daIns[, c("Y",  "Race", "Marital", "Nonres", "Injury", 
  "HuntYrs", "Age", "Inc")]
head(sam, n = 2)

# A. Quick summarization by row or column
# A1. Summary()
md <- summary(sam); md; str(md); is.matrix(md)
md2 <- substr(x = md[c(1, 4, 6), 1:6], start = 9, stop = 14)
md3 <- matrix(data = as.numeric(md2), nrow = 3)
rownames(md3) <- c("Minimum", "Mean", "Maximum")
colnames(md3) <- colnames(md2); md3

# A2. Single columns or rows by mean()
s1 <- sam$HuntYrs 
mean(s1); sd(s1); max(s1); min(s1); quantile(s1)
range(s1); unique(sam[, "Marital"]); length(s1)
min(5:1, 3.14); pmin(5:1, 3.14)  # Output: 1 versus 5 numbers
# determine if a column contains a dummy variable
identical(as.numeric(sort(unique(sam[, "Marital"]))), c(0, 1))

cumsum(1:4); cumprod(1:4); cummax(c(1, 5, 0, 8)); cummin(c(1, 5, 0, 8))
vec <- c(10, 50, 10, 50) 
which.min(vec); which.max(vec)  # The first extreme value
which(vec == min(vec))          # The minimum of 10 occurs twice.

# A3. Multiple columns or rows by colMeans()
ca <- colMeans(sam[, 1:6]); ca 
cb <- colSums(sam[, 1:6]); cb
rowMeans(sam[1:5, 2:4]); rowSums(sam[1:5, 2:4])

# B. Pivot table
# B1. Contingency table by table()
wa <- table(daIns[, c("Race", "Marital", "Nonres")]); wa
as.data.frame(wa)  # A nice change here
wb <- xtabs( ~ Race + Marital + Nonres, daIns); wb
as.data.frame(wb)
xtabs( ~ Race + Inc, daIns)
ftable(daIns[, c("Race", "Marital", "Nonres")])
ftable(daIns[, c("Race", "Marital", "Nonres")], row.vars = c(1, 2, 3))

# B2. aggregate() on data frame
ka <- aggregate(x = daIns$Y, 
  by = listn(daIns$Race, daIns$Marital), FUN = length)
str(ka); ka; as.data.frame(ka)
kb <- aggregate(x = daIns[, c("Y", "HuntYrs", "Inc")], 
  by = listn(daIns$Race, daIns$Marital), FUN = mean)
kc <- aggregate(x = daIns[, c("Y", "HuntYrs", "Inc")], 
  by = listn(daIns$Race, daIns$Marital), FUN = range)

# B3. aggregate() on time series
data(daBedRaw); head(daBedRaw)
kd <- aggregate(x = daBedRaw[, 1:3], nfrequency = 1, FUN = sum); head(kd)
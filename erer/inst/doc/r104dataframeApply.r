library(erer); data(daIns)
sam <- daIns[, c("Y",  "Race", "Marital", "Nonres", "Injury", 
  "HuntYrs", "Age", "Inc")]

# 1. apply(): over array margins (columns or rows)
qa <- apply(X = sam[, 3:6], MARGIN = 2, FUN = unique); qa
qb <- apply(X = sam[, 1:6], MARGIN = 2, FUN = mean); round(qb, digits = 2)
qc <- NULL
for (i in 1:6) {
  qc <- c(qc, mean(sam[, i]))
}
names(qc) <- colnames(sam[, 1:6])

# 2. lapply(), sapply(), mapply(): over a list or vector
# On a numeric vector; FUN has two arguments
price <- c(35, 10); names(price) <- c("Atlanta", "Boston")
lapply(X = price, FUN = rep, times = 4)
sapply(X = price, FUN = rep, times = 4)
mapply(FUN = rep, x = price, times = 4)
mapply(FUN = rep, x = price, times = 1:2)
mapply(FUN = rep, times = 1:2, MoreArgs = list(x = price))

# On a list from a data frame; FUN has one argument
lapply(X = listn(sam$Y, sam$Race, sam$HuntYrs), FUN = mean)
sapply(X = listn(sam$Y, sam$Race, sam$HuntYrs), FUN = mean)
mapply(FUN = mean, x = listn(sam$Y, sam$Race, sam$HuntYrs))

# On a data frame directly; FUN has two arguments
sam3 <- sam[, c("Y", "Race", "HuntYrs")]
lapply(X = sam3, FUN = quantile, probs = 0.5)
sapply(X = sam3, FUN = quantile, probs = 0.5)
mapply(FUN = quantile, x = sam3, probs = 0.5)
mapply(FUN = quantile, x = sam3, probs = c(0, 0.5, 1))
mapply(FUN = quantile, x = sam3, MoreArgs = list(probs = c(0, 0.5, 1)))

# 3. tapply(): over a ragged array with a factor index
# Create three factor vectors with labels
f.race <- factor(x = sam$Race, levels = c(0, 1), 
  labels = c("Other", "Caucasian"))
f.marital <- factor(x = sam$Marital, levels = c(0, 1), 
  labels = c("Unmarried", "Married"))
f.nonres <- factor(x = sam$Nonres, levels = c(0, 1), 
  labels = c("Resident", "Nonresident"))

# INDEX can be coerced into a factor internally
ta <- tapply(X = sam$HuntYrs, FUN = mean, 
  INDEX = listn(f.race, f.marital, f.nonres))
ta; as.data.frame(as.table(ta))  
tb <- tapply(X = sam$HuntYrs, FUN = mean, 
  INDEX = listn(sam$Race, sam$Marital, sam$Nonres))
tb; as.data.frame(as.table(tb))  
tc <- mean(sam[sam$Race == 1 & sam$Marital == 1 & sam$Nonres == 1, 
  "HuntYrs"])
tc
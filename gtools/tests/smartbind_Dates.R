library(gtools)

today <- Sys.Date()
tenweeks <- seq(today, length.out=10, by="1 week")

df1 <- data.frame(dates=tenweeks, chars=letters[1:10], ints=1:10, numeric=1.1:10.1)
df2 <- data.frame(chars=letters[11:20], ints=11:20, numeric=11.1:20.1)

smartbind(df1, df2)

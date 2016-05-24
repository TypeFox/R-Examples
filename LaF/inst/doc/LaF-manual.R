### R code from vignette source 'LaF-manual.Rnw'

###################################################
### code chunk number 1: LaF-manual.Rnw:47-72
###################################################
options(width=54)
library(LaF)

# Generate data
n <- 10000
data <- data.frame(
        id = trunc(runif(n, 1, 1E6)),
        gender = sample(c("M", "F"), n, replace=TRUE),
        postcode = paste(
            trunc(runif(n, 1000, 9999)),
            sample(LETTERS, n, replace=TRUE),
            sample(LETTERS, n, replace=TRUE), sep=""),
        age = round(runif(n, 0, 109)),
        income = round(rexp(n, 1/1000), 2),
    stringsAsFactors=FALSE)

# Generate fwf file
lines <- sprintf("%6.0f%1s%6s%3d%8.2f", data$id, data$gender, data$postcode,
    data$age, data$income)
writeLines(lines, con="file.fwf")

# Generate CSV file
lines <- sprintf("%.0f,%s,%s,%d,%f", data$id, data$gender, data$postcode,
    data$age, data$income)
writeLines(lines, con="file.csv")


###################################################
### code chunk number 2: LaF-manual.Rnw:127-129
###################################################
lines <- readLines("file.fwf", n=5)
cat(paste(lines,collapse="\n"), "\n")


###################################################
### code chunk number 3: LaF-manual.Rnw:132-137
###################################################
dat <- laf_open_fwf(filename="file.fwf", 
    column_types=c("integer", "categorical", 
        "string", "integer", "double"), 
    column_names=c("id", "gender", "postcode", "age", "income"),
    column_widths=c(6, 1, 6, 3, 8))


###################################################
### code chunk number 4: LaF-manual.Rnw:143-144
###################################################
alldata <- dat[ , ]


###################################################
### code chunk number 5: LaF-manual.Rnw:211-213
###################################################
lines <- readLines("file.csv", n=5)
cat(paste(lines,collapse="\n"), "\n")


###################################################
### code chunk number 6: LaF-manual.Rnw:216-220
###################################################
dat <- laf_open_csv(filename="file.csv", 
    column_types=c("integer", "categorical", 
        "string", "integer", "double"), 
    column_names=c("id", "gender", "postcode", "age", "income"))


###################################################
### code chunk number 7: LaF-manual.Rnw:226-227
###################################################
alldata <- dat[ , ]


###################################################
### code chunk number 8: LaF-manual.Rnw:240-241
###################################################
write_dm(dat, "model.yaml")


###################################################
### code chunk number 9: LaF-manual.Rnw:244-246
###################################################
lines <- readLines("model.yaml")
cat(paste(lines,collapse="\n"), "\n")


###################################################
### code chunk number 10: LaF-manual.Rnw:251-252
###################################################
dat <- laf_open(read_dm("model.yaml"))


###################################################
### code chunk number 11: LaF-manual.Rnw:277-278
###################################################
begin(dat)


###################################################
### code chunk number 12: LaF-manual.Rnw:284-285
###################################################
goto(dat, 1000)


###################################################
### code chunk number 13: LaF-manual.Rnw:291-293
###################################################
d <- next_block(dat)
nrow(d)


###################################################
### code chunk number 14: LaF-manual.Rnw:298-300
###################################################
d <- next_block(dat, columns=c(1,3), nrows=100)
dim(d)


###################################################
### code chunk number 15: LaF-manual.Rnw:311-319
###################################################
n <- 0
begin(dat)
while (TRUE) {
    d <- next_block(dat, 2)
    n <- n + sum(d$gender == 'M')
    if (nrow(d) == 0) break;
}
print(n)


###################################################
### code chunk number 16: LaF-manual.Rnw:333-338
###################################################
count <- function(d, prev) {
  if (is.null(prev)) prev <- 0
  return(prev + sum(d$gender == 'M'))
}
(n <- process_blocks(dat, count))


###################################################
### code chunk number 17: LaF-manual.Rnw:350-363
###################################################
ave <- function(d, prev) {
  # initialisation
  if (is.null(prev)) {
    prev <- c(sum=0, n=0)
  }
  # finilisation
  if (nrow(d) == 0) {
    return(as.numeric(prev[1]/prev[2]))
  }
  result <- prev + c(sum(d$income), nrow(d))
  return(result)
}
(n <- process_blocks(dat, ave, columns=5))


###################################################
### code chunk number 18: LaF-manual.Rnw:374-380
###################################################
# select the first 10 rows
result <- dat[1:10, ]
# select the second column
result <- dat[ , 2]
# select the first 10 rows and the second column
result <- dat[1:10, 2]


###################################################
### code chunk number 19: LaF-manual.Rnw:393-394
###################################################
result <- dat[dat$age[] > 65, ]


###################################################
### code chunk number 20: LaF-manual.Rnw:397-398
###################################################
result <- dat[dat[[4]][] > 65, ]


###################################################
### code chunk number 21: LaF-manual.Rnw:401-402
###################################################
result <- dat[dat[ , 4] > 65, ]


###################################################
### code chunk number 22: LaF-manual.Rnw:405-406
###################################################
result <- dat[dat[ , "age"] > 65, ]


###################################################
### code chunk number 23: LaF-manual.Rnw:417-419
###################################################
levels(dat)[["age"]] <- data.frame(levels=0:100, labels=paste(0:100, "years"))
dat$age[1:10]


###################################################
### code chunk number 24: LaF-manual.Rnw:424-427
###################################################
write_dm(dat, "model.yaml")
lines <- readLines("model.yaml")
cat(paste(c(lines[1:29], "..."),collapse="\n"), "\n")


###################################################
### code chunk number 25: LaF-manual.Rnw:452-454
###################################################
(m1 <- colmean(dat, columns=4))
(m1 <- colmean(dat$age))



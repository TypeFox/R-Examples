## ----, include=FALSE-----------------------------------------------------
library(rmarkdown)
library(knitr)
opts_chunk$set(results='markup', echo=TRUE)
knitr_options(opts_chunk=list(echo=TRUE, results='markup'))

## ------------------------------------------------------------------------
library(refset)
employees <- data.frame(
      id=1:4, 
      name=c("James", "Sylvia", "Meng Qi", "Luis"), 
      age=c(28,44,38, 23), 
      gender=factor(c("M", "F", "F", "M")), 
      stringsAsFactors=FALSE)
      
refset(rs, employees[1:2,])

## ------------------------------------------------------------------------
rs
employees$name[1] <- "Jimmy"
rs

## ------------------------------------------------------------------------
rs$age <- c(29, 45)
employees$age

## ------------------------------------------------------------------------
ss <- rs
employees$name[2] <- "Silvia"
rs$name[2]
ss$name[2]

## ------------------------------------------------------------------------
refset(rs2, rs$id)
rs2 
rs$id <- rs$id + 1000
rs2
rs2 <- 101:102
employees$id

## ------------------------------------------------------------------------
# the multi-argument form. Note the empty argument, to select all columns:
refset(rsd, employees, age < 30, , drop=FALSE)
rsd
employees$age <- employees$age + 1
rsd

## ------------------------------------------------------------------------
vec <- 1:10
refset(rs, vec, 4:6)
rs <- rs*10
vec

## ------------------------------------------------------------------------
lst <- list(a="text", b=42, NA)
refset(rsl, lst$b)
rsl <- "more text"
lst$b

## ------------------------------------------------------------------------
rs %r% employees[1:3,] # equivalent to refset(rs, employees[1:3,])

## ----wrap----------------------------------------------------------------
f <- function(x) {
  cx <- contents(x)
  contents(x)$name <- paste(cx$name, "the", sample(c("Kid", "Terrible", "Silent", 
        "Fair"), nrow(cx), replace=TRUE))
}

parcel <- wrap(rs)
f(parcel)
employees


## ------------------------------------------------------------------------
f <- function(x) {x <- x*2}
a <- 4
f(a)
a

## ------------------------------------------------------------------------
dfr <- data.frame(x1=1:5, x2=rnorm(5), alpha=letters[1:5])
refset(rs, dfr[dfr$x1 <= 3, c("x1", "alpha")])

## ------------------------------------------------------------------------
ss <- dfr[dfr$x1 <= 3, c("x1", "alpha")]
rs
ss

## ------------------------------------------------------------------------
c(class(rs), class(ss))
c(mean(rs$x1), mean(ss$x1))

## ------------------------------------------------------------------------
dfr$alpha <- c(NA, letters[23:26])
rs
ss

## ------------------------------------------------------------------------
rs$alpha <- LETTERS[1:3]
rs
dfr

## ------------------------------------------------------------------------
vec <- 1:10
refset(rvec, vec[2:3])
mylist <- list(a="some", b="more", c="data")
refset(rls, mylist$b)
refset(rls2, mylist[["c"]])
rvec
c(rls, rls2)

## ----, error=TRUE--------------------------------------------------------
myss <- subset(dfr, x1>1)
refset(rs, myss)

## ------------------------------------------------------------------------
top4 %r% dfr[1:4,]
exists("top4")

## ------------------------------------------------------------------------
refset(large, dfr, x2 > 0,)
large

## ------------------------------------------------------------------------
employees <- data.frame(
      id=1:4, 
      name=c("James", "Sylvia", "Meng Qi", "Luis"), 
      age=c(28,44,38, 23), 
      gender=factor(c("M", "F", "F", "M")),
      hours=c(160, 130, 185, 145),
      pay=c(60000, 50000, 70000, 60000),
      stringsAsFactors=FALSE)

## ------------------------------------------------------------------------
overtimers %r% employees[employees$hours > 140,]
overtimers

## ------------------------------------------------------------------------
employees$hours <- c(135, 150, 70, 145)
overtimers

## ------------------------------------------------------------------------
# people who worked long hours last month:
refset(overtimers_static, employees, hours > 140, , dyn.idx=FALSE)
# give them a holiday...
overtimers_static$hours <- 0
# ... and a pay rise
overtimers_static$pay <- overtimers_static$pay * 1.1 
overtimers_static

## ------------------------------------------------------------------------
copy <- overtimers
copy$pay <- copy$pay * 2
employees$pay # still the same :/

## ------------------------------------------------------------------------
rs %r% employees[1:3,]

## ----, ref.label='wrap'--------------------------------------------------
f <- function(x) {
  cx <- contents(x)
  contents(x)$name <- paste(cx$name, "the", sample(c("Kid", "Terrible", "Silent", 
        "Fair"), nrow(cx), replace=TRUE))
}

parcel <- wrap(rs)
f(parcel)
employees


## ------------------------------------------------------------------------
f <- function(parcel) {
  unwrap_as(emps, parcel)
  emps$name <- paste(emps$name, "the", sample(c("Kid", "Terrible", "Silent", 
        "Fair"), nrow(emps), replace=TRUE))
}
f(parcel)
employees

## ------------------------------------------------------------------------
parcel <- wrapset(employees, grepl("Terrible", employees$name), , drop=FALSE)
# note: drop=FALSE works just as for standard subsetting, see ?Extract
contents(parcel)


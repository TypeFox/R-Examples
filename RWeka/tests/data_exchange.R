
## test the data exchange between R and Weka

library(RWeka)

## low-level

set.seed(123)
x <- data.frame(R = runif(10),
                L = sample(c(FALSE, TRUE), 10, rep = TRUE),
                C = sample(LETTERS[1:3], 10, rep = TRUE),
                F = factor(sample(c("good","bad"), 10, rep = TRUE)),
                stringsAsFactors = FALSE)
x

jx <- RWeka:::read_data_into_Weka(x, length(x))
rJava::.jcall(jx, "I", "classIndex")

xj <- RWeka:::read_instances_from_Weka(jx)
all.equal(x, xj)    # TRUE

## mixed data

f <- system.file("arff", "test1.arff", package = "RWeka")
x <- read.arff(f)
x

temp <- tempfile()
write.arff(x, temp)

xx <- read.arff(temp)
unlink(temp)

all.equal(x, xx)    # TRUE

# test the R parser

xx <- RWeka:::read.arff.R(f)
xx

all.equal(x, xx)    # TRUE

## sparse data

f <- system.file("arff", "test_sparse.arff", package = "RWeka")
x <- read.arff(f)
x

temp <- tempfile()
write.arff(x, temp)

xx <- read.arff(temp)
unlink(temp)

all.equal(x, xx)    # TRUE

# test the R parser

all.equal(x, RWeka:::read.arff.R(f))	# TRUE

## uppercase and quoting

x <- read.arff(system.file("arff", "test2.arff", package = "RWeka"))
x

## connections

f <- file("")
write.arff(x, f)
xx <- read.arff(f)
close(f)
print(all.equal(x, xx))                 # TRUE

## date normalization
x <- data.frame(date1 = as.POSIXct("2012-12-12 12:12:12", tz = ""),
		date2 = as.POSIXct("2012-12-12 12:12:12", tz = "GMT"))
x

temp <- tempfile()
write.arff(x, temp)

xx <- read.arff(temp)
unlink(temp)
xx

# differ in representation
all.equal(x, xx)


## test the R parser and writer
temp <- tempfile()
RWeka:::write.arff.R(x, temp)

zz <- RWeka:::read.arff.R(temp)
unlink(temp)
zz

all.equal(xx, zz)   # TRUE 

###

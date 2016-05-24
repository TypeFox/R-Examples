context("convert_raw")

# load data with NA as missing
data(mouse_msats)
data <- mouse_msats

## constructing different data sets

# data with NA missing
data_NA <- data
# data with -1 missing
data[is.na(data)] <- -1
data_NUM <- data

# character data frame (for SNP data)
bp <- c('A','T','G','C')
m <- as.data.frame(matrix(data = bp, nrow = 100, ncol = 100))
m[] <- lapply(m, as.character)
m[] <- lapply(m, sample)
m[sample(1:100, 10), sample(1:100, 10)] <- NA
data_CHAR <- m

# convert raw
# different missing data
data_conv_NA <- convert_raw(data_NA)
# Characters as genotypes
data_conv_CHAR <- convert_raw(data_CHAR)

test_that("Converted data contains half the number columns", {
    expect_equal(ncol(convert_raw(data_NA)), ncol(data_NA)/2)
    expect_equal(ncol(convert_raw(data_CHAR)), ncol(data_CHAR)/2)
})

test_that("Converted data contains just 1, 0 and -1", {
    expect_equal(sum(!((data_conv_NA == 1) | (data_conv_NA == 0) | (is.na(data_conv_NA)))), 0)
    expect_equal(sum(!((data_conv_CHAR == 1) | (data_conv_CHAR == 0) | (is.na(data_conv_CHAR)))), 0)
})








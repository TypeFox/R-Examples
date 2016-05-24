## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library(brotli)
library(ggplot2)

# Example data
myfile <- file.path(R.home(), "COPYING")
x <- readBin(myfile, raw(), file.info(myfile)$size)

# The usual suspects
y1 <- memCompress(x, "gzip")
y2 <- memCompress(x, "bzip2")
y3 <- memCompress(x, "xz")
y4 <- brotli_compress(x)

## ------------------------------------------------------------------------
stopifnot(identical(x, memDecompress(y1, "gzip")))
stopifnot(identical(x, memDecompress(y2, "bzip2")))
stopifnot(identical(x, memDecompress(y3, "xz")))
stopifnot(identical(x, brotli_decompress(y4)))

## ------------------------------------------------------------------------
# Combine data
alldata <- data.frame (
  algo = c("gzip", "bzip2", "xz (lzma2)", "brotli"),
  ratio = c(length(y1), length(y2), length(y3), length(y4)) / length(x)
)

ggplot(alldata, aes(x = algo, fill = algo, y = ratio)) + 
  geom_bar(color = "white", stat = "identity") +
  xlab("") + ylab("Compressed ratio (less is better)")

## ------------------------------------------------------------------------
library(microbenchmark)
bm <- microbenchmark(
  memDecompress(y1, "gzip"),
  memDecompress(y2, "bzip2"),
  memDecompress(y3, "xz"),
  brotli_decompress(y4),
  times = 1000
)

alldata$decompression <- summary(bm)$median
ggplot(alldata, aes(x = algo, fill = algo, y = decompression)) + 
  geom_bar(color = "white", stat = "identity") +
  xlab("") + ylab("Decompression time (less is better)")

## ------------------------------------------------------------------------
library(microbenchmark)
bm <- microbenchmark(
  memCompress(x, "gzip"),
  memCompress(x, "bzip2"),
  memCompress(x, "xz"),
  brotli_compress(x),
  times = 20
)

alldata$compression <- summary(bm)$median
ggplot(alldata, aes(x = algo, fill = algo, y = compression)) + 
  geom_bar(color = "white", stat = "identity") +
  xlab("") + ylab("Compression time (less is better)")


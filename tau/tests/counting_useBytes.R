
### ceeboo 2008

library("tau")

.format.count <- function(x)
    data.frame(counts = x, bytes = nchar(names(x), type = "bytes", 
					 allowNA = TRUE),
               Encoding = Encoding(names(x)),
	       ## otherwise the code checks for illegal multibyte
               ## sequences and fails. great!
	       row.names = names(x),
	       stringsAsFactors = FALSE)

## latin capital letter a with diaresis

t1 <- c("abc", "a\xc4", "", NA)
Encoding(t1) <- c("unknow", "latin1", "unknown", "unknown") 
t1[2] <- iconv(t1[2], from = "latin1", to = "UTF-8")
t1
Encoding(t1)

t2 <- c(paste("_", t1[1:2], "_", sep = ""), t1[3:4])
t2

## count n-grams
r <- .Call("R_utf8CountNgram", list(t2), 3L, 0L, TRUE, FALSE, TRUE)
.format.count(r)

## incremental
.Call("R_utf8CountNgram", list(t2), 3L, 0L, TRUE, TRUE, TRUE)
r <- .Call("R_utf8CountNgram", list(t2), 3L, 0L, TRUE, FALSE, TRUE)
.format.count(r)

## count strings
r <- .Call("R_utf8CountString", list(t1), 3L, 0L, 0L, TRUE, FALSE, TRUE)
.format.count(r)

## count prefixes
r <- .Call("R_utf8CountString", list(t1), 3L, 0L, 1L, TRUE, FALSE, TRUE)
.format.count(r)

## count suffixes
r <- .Call("R_utf8CountString", list(t1), 3L, 0L, 2L, TRUE, FALSE, TRUE)
.format.count(r)

## FIXME add to interface
r <- .Call("R_utf8CountString", list(t1), 3L, 0L, 3L, TRUE, FALSE, TRUE)
.format.count(r)

## incremental
.Call("R_utf8CountString", list(t1), 3L, 0L, 0L, TRUE, TRUE, TRUE)
r <- .Call("R_utf8CountString", list(t1), 3L, 0L, 0L, TRUE, FALSE, TRUE)
.format.count(r)

###

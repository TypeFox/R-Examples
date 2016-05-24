
### ceeboo 2015

library("arulesSequences")

## use example
data(zaki)
s7 <- cspade(zaki, parameter = list(support = .25), 
	           control   = list(verbose = TRUE, tidLists = TRUE))

k <- support(s7, zaki, control = list(verbose = TRUE))
table(size(s7), sign(quality(s7)$support - k))

## random atomic sequences
set.seed(20150921)

f <- tempfile()
local({
n <- 1000
l <- 30
k <- 5

s <- lapply(sample(1:l, n, TRUE), function(x)
       sample(1:k, x, TRUE)
)

m <- sapply(s, length)
s <- cbind(
    rep(1:n, m),
    unlist(lapply(m, seq)),
    unlist(s)
)

write.table(s, file = f, row.names = FALSE, col.names = FALSE)
})

## use generated data
t <- read_baskets(con  = f,
                  info = c("sequenceID", "eventID"))
unlink(f)

## use low support
s1 <- cspade(t, parameter = list(support = .17), 
                control   = list(verbose = TRUE))
summary(s1)

##
k <- support(s1, t, control = list(verbose = TRUE))
table(size(s1), sign(quality(s1)$support - k))

## internal
all.equal(
    .Call(arulesSequences:::R_pnscount, s1@data, s1@data, 
	  s1@elements@items@data, FALSE),
    .Call(arulesSequences:::R_pnscount, s1@data, s1@data, NULL, TRUE)
)

all.equal(
    .Call(arulesSequences:::R_pnsclosed, s1@data,
	  s1@elements@items@data, rep(1L, length(s1)), TRUE),
    is.maximal(s1)
)


all.equal(
    .Call(arulesSequences:::R_pnsredundant, s1@data,
	  s1@elements@items@data, rep(1L, length(s1)), TRUE),
    size(s1) > 1L
)

warnings()
###


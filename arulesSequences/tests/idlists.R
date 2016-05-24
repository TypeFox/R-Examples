
### ceeboo 2015

library("arulesSequences")

.string2sequences <-
function(x) {
    x <- lapply(strsplit(x, split = " "), strsplit, split = "")
    x <- unlist(x, recursive = FALSE)
    as(x, "sequences")
}

## repeating
t <- paste(
    1, 
    c(10, 15, 18, 19), 
    c("A", "A", "A", "A"),
    collapse = "\n"
)
t <- textConnection(t)
t <- read_baskets(t, info = c("sequenceID", "eventID"))
as(t, "data.frame")

## The second pattern has two matches starting at 10 and 15
## extending over 8 and 4 time periods. 
s <- c("A AAA")
s <- .string2sequences(s)

k <- 
list(
    list(parameter = list()),
    ## Test if cummulation of gaps works.
    list(parameter = list(maxwin = 5))
)

k <- mapply(function(x)
    support(s, t, control = x),
    k
)

cbind(as(s, "data.frame"), support = k)

## Test if optimization works.
s <- c("A AB")
s <- .string2sequences(s)
k <- support(s, t, control = list(verbose = TRUE, parameter = list()))
cbind(as(s, "data.frame"), support = k)

## Test if conversion works.
all.equal(k, support(s, as(t, "timedsequences"),
		     control = list(parameter = list())))

k <- support(s, s, type = "absolute", control = list(parameter = list()))
all.equal(k, as.integer(rowSums(is.subset(s))))

## Test internal
k <- arulesSequences:::support.idlists(s, t, type = "tidLists")
all.equal(k, supportingTransactions(s, t))

##

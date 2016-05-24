require("simsalapar")

sessionInfo()# see it interactively

### A dummy simulation, illustrating how the .Random.seeds are set

vlist <- varlist(
    ## replications
    n.sim = list(value = 7), ## <- small, just for illustration
    ## sample size
    n = list(type="grid", value = c(20, 40, 100)),
    method = list(type = "grid", ## complete dummy; methods *equal on purpose*:
        value = list( ## alternatively, return runif(1) or (4) ..
            f = function() .Random.seed,
            g = function() .Random.seed)))

mkGrid(vlist)# simplistic visual check

do1 <- function(n, method) if(is.function(method)) method() else if(is.list(method)) method[[1]]()

options(mc.cores=2)
res <- doClusterApply(vlist, doOne=do1, monitor=TRUE)

str(a <- getArray(res))
 ## int [1:626, 1:3, 1:2, 1:7] 403 624 -169270483 -442010614 ..........
 ## - attr(*, "dimnames")=List of 4
 ##  ..$ D.1   : NULL
 ##  ..$ n     : chr [1:3] "20" "40" "100"
 ##  ..$ method: chr [1:2] "f" "g"
 ##  ..$ n.sim : NULL
all.eq.1st <- function(x) all(x == x[1])

str(ae <- apply(a, c("D.1", "n.sim"), all.eq.1st))
## 626 x 7 =  length(.Random.seed) x n.sim
stopifnot(all(ae))
## --> so indeed, the .Random.seeds are all equal within one simulation, i.e.,
## above for the different 'n' and 'method'

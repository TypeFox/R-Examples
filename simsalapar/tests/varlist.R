require("simsalapar")

## as long as R <= 3.0.1 is possible:
source(system.file("xtraR/assertErr-etc.R", package="simsalapar", mustWork=TRUE))

vli <- varlist(
    ## replications
    ## suposed to be 500, choosen 2 here for performance purposes
    n.sim=list(expr=quote(N[sim]),type="N", value=1),
    ## sample size
    n = list(type="grid", value = c(50,100)),
    ## Beta
    beta = list(type="grid",value=list(beta1=0:1, beta2=c(3,5), beta3=c(1,10))),
    # Method
    method = list(type="grid",value=list(g1= function(x,y) x+y)),
    # formula
    formula = list(type="frozen", value= y~x.1+x.2+x.3+x.4+x.5),
    ## outlier locations
    kstep = list(type="inner",value=seq(0,100,by=20))
    )
vl <- vli
vl $ kstep $type <- "grid"

ws <- warnings()
##-> Warning should end in " :..method"
stopifnot(inherits(ws, "warnings"),
	  grepl("method", sub(".*:", "",
			      names(ws)[[1]])))## wrongly had "beta"

## now works without error:
toLatex(vl, label = "tab:var", caption = "varlist 'vl' - regression")
## TODO MM: 1) formula  2) extra "," in  "function (x, y) , x + y"


doNil <- function(...) rpois(3, lambda=10)

doNamed <- function(...) setNames(rpois(4, lambda=10), LETTERS[1:4])

doMatr	 <- function(...) matrix(rpois(6, lambda=10), 2,3)
doMat.wrong <- function(n, method, beta, kstep, formula)
    matrix(rpois(6, lambda=10), length(kstep), 11)
doMat.i	 <- function(n, method, beta, kstep, formula)
    matrix(rpois(6, lambda=10), 11, length(kstep))

doMat.n <- function(...)
    matrix(rpois(6, lambda=10), 2,3,
	   dimnames=list(paste0("r", 1:2), paste0("C.",1:3)))
doMat.in <- function(n, method, beta, kstep, formula) {
    nr <- length(kstep)
    matrix(rpois(6, lambda=10), 2, nr,
	   dimnames=list(paste0("r", 1:2), paste0("C.",1:nr)))
}

doArr1 <- function(...) array(rpois(7, lambda=10), 7)
doArr1n <- function(...) array(rpois(7, lambda=10), 7, dimnames=list(letters[1:7]))
doArr1.i <- function(n, method, beta, kstep, formula) {
    n <- length(kstep); array(rpois(n, lambda=10), n) }
doArr3n <- function(...)
  array(rpois(24, lambda=10), c(3,4,2), dimnames=
	list(letters[1:3],paste0("d",1:4), paste0("D",1:2)))

doArr3 <- function(...) array(rpois(24, lambda=10), c(3,4,2))
doArr3.n <- function(...)
    array(rpois(30, lambda=10), c(3,5,2), dimnames=
	  list(letters[1:3],NULL, paste0("D",1:2)))
doArr3.i <- function(n, method, beta, kstep, formula)  {
    n <- length(kstep)
    array(rpois(2*4*n, lambda=10), c(2,4,n)) }

##' basic structural check or do*Apply() result 'array-list' {as made by mkAL() }
chkArray <- function(res) {
    stopifnot(is.list(res), !is.null(d <- dim(res)), prod(d) == length(res))
    v1 <- res[[1]]$value

    val   <- getArray(res) # array of values
    err   <- getArray(res, "error") # array of error indicators
    pwarn <- getArray(res, "warning") # array of warning indicators
    time  <- getArray(res, "time") # array of user times in ms
    Dim <- function(x) if(is.null(d <- dim(x))) length(x) else d
    stopifnot(is.numeric(val),
	      all.equal(c(Dim(v1), dim(err)), dim(val), tol=0),
	      identical(dim(err), dim(pwarn)),
	      identical(dim(err), dim(time)))
}

###--------- n.sim missing ( = 1 ) : -----------------------------------------

## 1. "Nil": unnamed vector ------------------------------
chkArray(ra1   <- doLapply(vl, seed="seq", sfile=NULL, doOne= doNil))
## 2. "Named": *named* vector ------------------------------ had bug
chkArray(raN   <- doLapply(vl, seed="seq", sfile=NULL, doOne= doNamed))
## 3. Arrays (including matrices)
##   1D
chkArray(raA1  <- doLapply(vl, seed="seq", sfile=NULL, doOne= doArr1))
chkArray(raA1n <- doLapply(vl, seed="seq", sfile=NULL, doOne= doArr1n))
chkArray(raA1.i<- doLapply(vli,seed="seq", sfile=NULL, doOne= doArr1.i))
##   2D
chkArray(raM   <- doLapply(vl, seed="seq", sfile=NULL, doOne= doMatr))
chkArray(raM.  <- doLapply(vl, seed="seq", sfile=NULL, doOne= doMat.n))
assertError(doLapply(vli,seed="seq", sfile=NULL, doOne= doMat.wrong))
chkArray(raM.i <- doLapply(vli,seed="seq", sfile=NULL, doOne= doMat.i))
chkArray(raMin <- doLapply(vli,seed="seq", sfile=NULL, doOne= doMat.in))

##   3D
chkArray(raA3  <- doLapply(vl, seed="seq", sfile=NULL, doOne= doArr3))
chkArray(raA3n <- doLapply(vl, seed="seq", sfile=NULL, doOne= doArr3n))
chkArray(raA3n.<- doLapply(vl, seed="seq", sfile=NULL, doOne= doArr3.n))
chkArray(raA3.i<- doLapply(vli,seed="seq", sfile=NULL, doOne= doArr3.i))

###--------- n.sim = 3 : -----------------------------------------
vl3 <- set.n.sim(vl, 3)

## 1. "Nil": unnamed vector ------------------------------
chkArray(sa1   <- doLapply(vl3, seed="seq", sfile=NULL, doOne= doNil))
## 2. "Named": *named* vector ------------------------------ had bug
chkArray(saN   <- doLapply(vl3, seed="seq", sfile=NULL, doOne= doNamed))
## 3. Arrays (including matrices)
##   1D
chkArray(saA1  <- doLapply(vl3, seed="seq", sfile=NULL, doOne= doArr1))
chkArray(saA1n <- doLapply(vl3, seed="seq", sfile=NULL, doOne= doArr1n))
##   2D
chkArray(saM   <- doLapply(vl3, seed="seq", sfile=NULL, doOne= doMatr))
chkArray(saM.  <- doLapply(vl3, seed="seq", sfile=NULL, doOne= doMat.n))
##   3D
chkArray(saA3  <- doLapply(vl3, seed="seq", sfile=NULL, doOne= doArr3))
chkArray(saA3n <- doLapply(vl3, seed="seq", sfile=NULL, doOne= doArr3n))
chkArray(saA3n.<- doLapply(vl3, seed="seq", sfile=NULL, doOne= doArr3.n))

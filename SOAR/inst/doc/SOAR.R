### R code from vignette source 'SOAR.Rnw'

###################################################
### code chunk number 1: SOAR.Rnw:78-88
###################################################
## attach the package, checking that it is installed
stopifnot(require("SOAR"))

## create some dummy data
X <- matrix(rnorm(1000*50), 1000, 50)
S <- var(X)
Xb <- colMeans(X)

## and store part of it in the default cache
Store(X)


###################################################
### code chunk number 2: SOAR.Rnw:96-98
###################################################
objects()   ## or ls()
find("X")


###################################################
### code chunk number 3: SOAR.Rnw:125-128
###################################################
Store(Xb, S)
Search()
Ls()        # short alias for Objects()


###################################################
### code chunk number 4: SOAR.Rnw:153-154
###################################################
vc <- gc(); vc ; v0 <- vc["Vcells", "used"]


###################################################
### code chunk number 5: SOAR.Rnw:164-172
###################################################
Vcells <- function()
  c(Vcells = gc()["Vcells", "used"])            ; Vcells()-v0
bigX <- matrix(rnorm(1000^2), 1000, 1000)       ; Vcells()-v0
Store(bigX)                                     ; Vcells()-v0
d <- dim(bigX)                                  ; Vcells()-v0
bigX[1,1] <- 0                                  ; Vcells()-v0
Attach()                                        ; Vcells()-v0
Store(bigX)                                     ; Vcells()-v0


###################################################
### code chunk number 6: SOAR.Rnw:197-198
###################################################
Remove(bigX)                                   ; (Vcells()-v0)


###################################################
### code chunk number 7: SOAR.Rnw:225-226 (eval = FALSE)
###################################################
## Store(objects())


###################################################
### code chunk number 8: SOAR.Rnw:232-234 (eval = FALSE)
###################################################
## objs <- objects()
## Store(list = objs)


###################################################
### code chunk number 9: SOAR.Rnw:240-242 (eval = FALSE)
###################################################
## objs <- ls()
## Store(objs, list = objs)


###################################################
### code chunk number 10: SOAR.Rnw:247-248 (eval = FALSE)
###################################################
## Remove(Objects())


###################################################
### code chunk number 11: SOAR.Rnw:273-275 (eval = FALSE)
###################################################
## cat("\nR_LOCAL_CACHE=.R_Store\n",
##     file = "~/.Renviron", append=TRUE)


###################################################
### code chunk number 12: SOAR.Rnw:288-290 (eval = FALSE)
###################################################
## Attach(lib = ".R_Store")
## Attach(lib = .R_Store)


###################################################
### code chunk number 13: SOAR.Rnw:294-296 (eval = FALSE)
###################################################
## lib <- ".R_Store"
## Store(X, Y, Z, lib = lib)


###################################################
### code chunk number 14: SOAR.Rnw:319-320 (eval = FALSE)
###################################################
## Attach(lib.loc = "..")


###################################################
### code chunk number 15: SOAR.Rnw:375-376 (eval = FALSE)
###################################################
## `%ni%` <- Negate(`%in%`)


###################################################
### code chunk number 16: SOAR.Rnw:384-385 (eval = FALSE)
###################################################
## lsCache <- Objects


###################################################
### code chunk number 17: SOAR.Rnw:390-391 (eval = FALSE)
###################################################
## StoreUtils(Vcells, `%ni%`, lsCache)


###################################################
### code chunk number 18: SOAR.Rnw:395-396 (eval = FALSE)
###################################################
## AttachUtils()


###################################################
### code chunk number 19: SOAR.Rnw:405-406 (eval = FALSE)
###################################################
## AttachData()


###################################################
### code chunk number 20: SOAR.Rnw:488-493 (eval = FALSE)
###################################################
## if(require("SOAR")) {
##     lst <- paste('autoload("', objects("package:SOAR"),
##                  '", "SOAR")\n', sep="")
##     cat("\n", lst, sep="", file = "~/.Rprofile", append = TRUE)
## }


###################################################
### code chunk number 21: SOAR.Rnw:547-551 (eval = FALSE)
###################################################
## x <- 1
## Store(x, lib = ".A")
## x <- 2
## Store(x, lib = ".a")


###################################################
### code chunk number 22: SOAR.Rnw:671-672 (eval = FALSE)
###################################################
## Remove(Objects())


###################################################
### code chunk number 23: SOAR.Rnw:692-699
###################################################
bad <- c(" ", "<", ">", ":", "\"", "/", "\\", "|", "?", "*")
rpl <- paste("@", 0:9, sep = "")

out <- rbind("Code:" = rpl, "Character:" = bad)
colnames(out) <- rep("    ", 10)
noquote(format(out, justify = "right"))
rm(bad, rpl, out)



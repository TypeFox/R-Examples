library(VLMC)

set.seed(21001)
c0 <- proc.time()

f1 <- c(1,0,0,0)
f2 <- rep(1:0,2)
(dt2 <- rep(c(f1,f1,f2,f1,f2,f2,f1),2))
dt3 <- c(unlist(sapply(1:5, function(k) c(rep(0,k),1))),
         rep(1,6))# patterns NOI in dt2

(vlmc.dt2c15  <- vlmc(dt2, cutoff = 1.5))
draw(vlmc.dt2c15)
## Predict "myself" (in sample)
noquote(cbind(fit   = predict(vlmc.dt2c15, type = "response"),
              probs = pr <- predict(vlmc.dt2c15),
              lev   = predict(vlmc.dt2c15, type = "depth"),
              ID    = (pid <- predict(vlmc.dt2c15, type = "id")),
              flags = {ff <- attr(pr,"flags"); ff[ff == "0"] <- ""; ff},
              context= id2ctxt(pid, m=2, alpha=TRUE)
              ))

(vlmc.dt2c.1 <- vlmc(dt2, cutoff = .1))
(vlmc.dt2c0  <- vlmc(dt2, cutoff = 0, thresh = 1))# need thresh=1 for "perfect fit"

draw(vlmc.dt2c.1)

## Predict "myself" (in sample)
pp.1  <- predict(vlmc.dt2c.1)# pred.prob
p.1   <- predict(vlmc.dt2c.1, type="response")
pid.1 <- predict(vlmc.dt2c.1, type="id.node") ## now give "-" for initial ones!
(pl.1 <- predict(vlmc.dt2c.1, type="depth"))

noquote(cbind(pp.1,fit=p.1, "="=c("<>","")[1+(dt2==p.1)], lev= pl.1,
              id= pid.1, flags = attr(p.1,"flags"),
              context = id2ctxt(pid.1, m=2, alpha=TRUE)
))
## Almost the same :
predict(vlmc.dt2c.1, type = "ALL")


draw(vlmc.dt2c0)
draw(vlmc.dt2c0, cumulative = FALSE)

## Predict "myself" (in sample)
(p0 <- predict(vlmc.dt2c0, type="response"))
all(dt2[-1] == abs(p0[-1])) # TRUE !
pp     <- predict(vlmc.dt2c0)# pred.prob
!any(is.na(match(pp[-1,],0:1)))#> TRUE
## since all pred.probabilties are 0 or 1 !

p.id   <- predict(vlmc.dt2c0, type="id.node")
p.lev <- predict(vlmc.dt2c0, type="depth")
noquote(cbind(pp, flags = attr(p.id,"flags"),
              ID= p.id, lev = p.lev, ctxt = id2ctxt(p.id,m=2, a=TRUE)))
## or almost the same:
predict(vlmc.dt2c0, type="ALL")

## This once gave a seg.fault!
p.dt2c0dt3 <- predict(vlmc.dt2c0, dt3) ## quite many NA (but not at end!)
## which(is.na(p.dt2c0dt3[,1]))
which(is.na(unname(p.dt2c0dt3[-1,1])))
##  1  3  6 10 15 ---- SAME as ../VLMC/dt2_c0.dt3.pred !!
## NONE!

predict(vlmc.dt2c0, dt3, type = "ALL")

## NA NA's :
prm <- cbind(fit   = predict(vlmc.dt2c0, dt3, type = "response"),
              probs = (pr  <- predict(vlmc.dt2c0, dt3)),
              lev   =         predict(vlmc.dt2c0, dt3, type = "depth"),
              ID    = (pid <- predict(vlmc.dt2c0, dt3, type = "id")),
              flags = attr(pr,"flags"))
noquote(cbind(prm, context = id2ctxt(pid, m=2, alpha=TRUE)))

## Constant data prediction:  |alphabet| = 1 (!!)
mod <- vlmc(rep(1, 99))
pr <- predict(mod, c(1, 0, 1))  # crash did occur here (in C code)
stopifnot(pr[-1] == 1)


cat("Time elapsed:", format(proc.time() - c0),"\n")

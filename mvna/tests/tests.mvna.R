require(mvna)

### First test with fake data

time <- 1:10
from <- rep(0, 10)
to <- c(1, 1, "cens", 1, "cens", rep(1, 5))
id <- 1:10
d1 <- data.frame(id, from, to, time)

tra1 <- matrix(FALSE, 2, 2)
tra1[1, 2] <- TRUE

test1 <- mvna(d1, c("0", "1"), tra1, "cens")

tna <- cumsum(as.numeric(to != "cens") / 10:1)
tvarna <- cumsum(as.numeric(to != "cens") / (10:1)^2)

all.equal(c(0, tna), test1$"0 1"$na)
all.equal(c(0, tvarna), test1$"0 1"$var1)

## same test but with left-truncation

from <- rep(0, 10)
to <- c(1, 1, "cens", 1, "cens", rep(1, 5))
entry <- c(rep(1, 3), rep(2, 3), rep(4, 4))
exit <- c(2, 7, 7, 4, 6, 9, 7, 9, 6, 7)
id <- 1:10
d2 <- data.frame(id, from, to, entry, exit)

test2 <- mvna(d2, c("0", "1"), tra1, "cens")

tna <- cumsum(c(1, 1, 1, 3, 2) /
              c(3, 5, 8, 6, 2))

all.equal(test2$"0 1"$na, c(0, tna))

## with a fake illness-death model with recovery

id <- c(1:20, 11:20)
from <- c(rep(0, 20), rep(1, 10))
to <- c(rep(2, 5), rep("cens", 5), rep(1, 10), rep(2, 5), rep("cens", 5))
entry <- c(rep(0, 20), rep(10, 10))
exit <- c(rep(9, 5), rep(8, 5), rep(10, 10), rep(15, 5), rep(13, 5))
d3 <- data.frame(id, from, to, entry, exit)

tra2 <- matrix(FALSE, 3, 3)
tra2[1, 2:3] <- TRUE
tra2[2, 3] <- TRUE

test3 <- mvna(d3, c("0", "1", "2"), tra2, "cens")

tna01 <- c(0, cumsum(c(0, 0, 10) / c(20, 15, 10)))
tna02 <- c(0, cumsum(c(0, 5, 0) / c(20, 15, 10)))
tna12 <- cumsum(c(0, 5) / c(10, 5))

all.equal(tna01, test3$"0 1"$na)
all.equal(tna02, test3$"0 2"$na)
all.equal(tna12, test3$"1 2"$na)

### tests with the included data sets

## sir.adm
data(sir.adm)

## data set transformation
data(sir.adm) 
id <- sir.adm$id
from <- sir.adm$pneu
to <- ifelse(sir.adm$status==0,"cens",sir.adm$status+1)
times <- sir.adm$time
dat.sir <- data.frame(id,from,to,time=times)

## Possible transitions
tra <- matrix(ncol=4,nrow=4,FALSE)
tra[1:2,3:4] <- TRUE

na.pneu <- mvna(dat.sir,c("0","1","2","3"),
                tra,"cens")

na.pneu

na.pneu$"0 2"$na
na.pneu$"0 2"$var.aalen
na.pneu$"0 2"$var.greenwood

na.pneu$"1 2"$na
na.pneu$"1 2"$var.aalen
na.pneu$"1 2"$var.greenwood

na.pneu$n.risk


## sir.cont
data(sir.cont)

## Matrix of possible transitions
tra <- matrix(ncol=3,nrow=3,FALSE)
tra[1, 2:3] <- TRUE
tra[2, c(1, 3)] <- TRUE

## Modification for patients entering and leaving a state
## at the same date
sir.cont <- sir.cont[order(sir.cont$id, sir.cont$time), ]
for (i in 2:nrow(sir.cont)) {
    if (sir.cont$id[i]==sir.cont$id[i-1]) {
        if (sir.cont$time[i]==sir.cont$time[i-1]) {
            sir.cont$time[i-1] <- sir.cont$time[i-1] - 0.5
        }
    }
}

## Computation of the Nelson-Aalen estimates
na.cont <- mvna(sir.cont,c("0","1","2"),tra,"cens")

na.cont

na.cont$"0 1"$na
na.cont$"0 1"$var.aalen
na.cont$"0 1"$var.greenwood

na.cont$"0 2"$na
na.cont$"0 2"$var.aalen
na.cont$"0 2"$var.greenwood

na.cont$"1 2"$na
na.cont$"1 2"$var.aalen
na.cont$"1 2"$var.greenwood

na.cont$n.risk
na.cont$n.event

summ.na.cont <- summary(na.cont)

all.equal(summ.na.cont$"0 1"$na, na.cont$"0 1"$na)
all.equal(summ.na.cont$"0 1"$var.aalen, na.cont$"0 1"$var.aalen)

aa <- predict(na.cont, tr.choice = "0 1", times = 9.5)

aa

## abortion data
data(abortion)

## Data set modification in order to be used by mvna
names(abortion) <- c("id", "entry", "exit", "from", "to")
abortion$to <- abortion$to + 1

## computation of the matrix giving the possible transitions
tra <- matrix(FALSE, nrow = 5, ncol = 5)
tra[1:2, 3:5] <- TRUE

na.abortion <- mvna(abortion, as.character(0:4), tra, NULL)

na.abortion$"0 2"$na
na.abortion$"0 2"$var.aalen
na.abortion$"0 2"$var.greenwood

na.abortion$"1 2"$na
na.abortion$"1 2"$var.aalen
na.abortion$"1 2"$var.greenwood

na.abortion$"1 3"$na
na.abortion$"1 3"$var.aalen
na.abortion$"1 3"$var.greenwood

na.abortion$"0 3"$na
na.abortion$"0 3"$var.aalen
na.abortion$"0 3"$var.greenwood

na.abortion$n.risk

na.abortion$n.event

na.abortion


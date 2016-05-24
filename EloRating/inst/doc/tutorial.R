### R code from vignette source 'tutorial.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: tutorial.Rnw:23-25
###################################################
options(width=80)
options(continue=" ")


###################################################
### code chunk number 2: "eloinstall2" (eval = FALSE)
###################################################
## install.packages("EloRating")


###################################################
### code chunk number 3: reading1
###################################################
library(EloRating)
xdata <- read.table(system.file("ex-sequence.txt", package = 'EloRating'),
                    header=T)


###################################################
### code chunk number 4: "reading1b" (eval = FALSE)
###################################################
## xdata <- read.table("c:\\temp\\ex-sequence.txt", 
##                     header=TRUE, sep="\t")


###################################################
### code chunk number 5: "seqcheck"
###################################################
seqcheck(winner=xdata$winner, loser=xdata$loser, Date=xdata$Date)


###################################################
### code chunk number 6: "eloseq" (eval = FALSE)
###################################################
## res <- elo.seq(winner=xdata$winner, loser=xdata$loser, Date=xdata$Date,
##                runcheck=TRUE)


###################################################
### code chunk number 7: "eloseqhidden"
###################################################
res <- elo.seq(winner=xdata$winner, loser=xdata$loser, Date=xdata$Date,
               runcheck=TRUE)


###################################################
### code chunk number 8: summary1
###################################################
summary(res)


###################################################
### code chunk number 9: extract.elo
###################################################
extract.elo(res, "2000-05-28")
extract.elo(res, "2000-05-28", IDs=c("s", "a", "c", "k"))


###################################################
### code chunk number 10: fig1plot
###################################################
eloplot(res)


###################################################
### code chunk number 11: fig1
###################################################
eloplot(res)


###################################################
### code chunk number 12: fig2plot
###################################################
eloplot(res,ids=c("s", "a", "w", "k", "c"), 
        from="2000-06-05", to="2000-07-04")


###################################################
### code chunk number 13: fig2
###################################################
eloplot(res,ids=c("s", "a", "w", "k", "c"), 
        from="2000-06-05", to="2000-07-04")


###################################################
### code chunk number 14: "reading2"
###################################################
#xpres <- read.table("ex-presence.txt", 
                    #header=T, sep="\t")
xpres <- read.table(system.file("ex-presence.txt", package = 'EloRating'), header=T)

xpres[,1] <- as.Date(as.character(xpres[,1]))


###################################################
### code chunk number 15: "seqcheck2"
###################################################
seqcheck(winner=xdata$winner, loser=xdata$loser, 
         Date=xdata$Date, presence=xpres)


###################################################
### code chunk number 16: "eloseqhidden2"
###################################################
res2 <- elo.seq(winner=xdata$winner, loser=xdata$loser, 
               Date=xdata$Date, presence=xpres, draw=xdata$Draw)


###################################################
### code chunk number 17: "eloseq2" (eval = FALSE)
###################################################
## res2 <- elo.seq(winner=xdata$winner, loser=xdata$loser, 
##                Date=xdata$Date, presence=xpres, draw=xdata$Draw)


###################################################
### code chunk number 18: extract.elo
###################################################
extract.elo(res2, "2000-05-28")
# note that "s" is absent and omitted
extract.elo(res2, "2000-05-28", IDs=c("s", "a", "c", "k"))
# note that "s" is absent and returned as NA


###################################################
### code chunk number 19: fig3plot
###################################################
eloplot(res2)


###################################################
### code chunk number 20: fig3
###################################################
eloplot(res2)


###################################################
### code chunk number 21: fig4plot
###################################################
eloplot(res2, ids=c("s", "a", "w", "k", "c"), 
        from="2000-06-05", to="2000-07-04")


###################################################
### code chunk number 22: fig4
###################################################
eloplot(res2, ids=c("s", "a", "w", "k", "c"), 
        from="2000-06-05", to="2000-07-04")


###################################################
### code chunk number 23: stabil
###################################################
stab.elo(res2, from="2000-05-05", to="2000-06-05")


###################################################
### code chunk number 24: traj
###################################################
traj.elo(res2, ID=c("f", "n"),
         from="2000-05-05", to="2000-06-05")


###################################################
### code chunk number 25: indiv
###################################################
individuals(res2, from="2000-05-05", to="2000-05-05", outp="N")
individuals(res2, from="2000-05-05", to="2000-06-05", outp="N")
individuals(res2, from="2000-05-05", to="2000-06-05", outp="CV")
individuals(res2, from="2000-05-05", to="2000-06-05", outp="IDs")


###################################################
### code chunk number 26: winprob
###################################################
winprob(1000,1200)
winprob(1200,1000)
winprob(1200,1200)


###################################################
### code chunk number 27: mat
###################################################
creatematrix(res2)
sum(creatematrix(res2))
creatematrix(res2, drawmethod="0.5")
sum(creatematrix(res2, drawmethod="0.5"))
# "c" and "n" are omitted
creatematrix(res2, c("2000-06-10", "2000-06-16"))
creatematrix(res2, c("2000-06-10", "2000-06-16"), 
             onlyinteracting=TRUE)


###################################################
### code chunk number 28: simu1
###################################################
rdata <- randomsequence()


###################################################
### code chunk number 29: simu2 (eval = FALSE)
###################################################
## xres <- elo.seq(rdata$seqdat$winner, rdata$seqdat$loser, 
##                 rdata$seqdat$Date, presence=rdata$pres)


###################################################
### code chunk number 30: simu3
###################################################
xres <- elo.seq(rdata$seqdat$winner, rdata$seqdat$loser, 
                rdata$seqdat$Date, presence=rdata$pres)


###################################################
### code chunk number 31: simu4
###################################################
summary(xres)


###################################################
### code chunk number 32: DS
###################################################
data(bonobos)
DS(bonobos)


###################################################
### code chunk number 33: tutorial.Rnw:306-308
###################################################
data(bonobos)
xdata <- randomelo(bonobos, 5)


###################################################
### code chunk number 34: tutorial.Rnw:310-312 (eval = FALSE)
###################################################
## data(bonobos)
## xdata <- randomelo(bonobos, 100)


###################################################
### code chunk number 35: tutorial.Rnw:314-315
###################################################
res <- data.frame(ID=colnames(xdata[[1]]), avg=round(colMeans(xdata[[1]]),1))


###################################################
### code chunk number 36: tutorial.Rnw:319-321
###################################################
ds <- DS(bonobos)
ds <- ds[order(ds$ID), ]


###################################################
### code chunk number 37: fig5plot
###################################################
plot(ds$normDS, res$avg, xlab="David's score", 
     ylab="randomized average Elo rating", las=1)


###################################################
### code chunk number 38: fig5
###################################################
plot(ds$normDS, res$avg, xlab="David's score", 
     ylab="randomized average Elo rating", las=1)


###################################################
### code chunk number 39: tutorial.Rnw:342-345
###################################################
data(adv); data(advpres)
x <- elo.seq(winner=adv$winner, loser=adv$loser, Date=adv$Date, 
             presence=advpres)


###################################################
### code chunk number 40: tutorial.Rnw:347-350
###################################################
prunk(x, c("2010-01-01", "2010-01-15"))
mat <- creatematrix(x, c("2010-01-01", "2010-01-15"))
prunk(mat)



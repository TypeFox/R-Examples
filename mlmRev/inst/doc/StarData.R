### R code from vignette source 'StarData.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(show.signif.stars = FALSE,width=68)
library(mlmRev)
library(lme4)


###################################################
### code chunk number 2: datainput
###################################################
str(orig <- read.delim(unz(system.file("original/text-star.zip",
                                       package = "mlmRev"),
                           "webstar.txt"),
                       colCl = rep(c("factor", "integer"), c(5, 48))))


###################################################
### code chunk number 3: strorig
###################################################
str(orig)


###################################################
### code chunk number 4: missval
###################################################
mv <- rep("9", 53)
mv[c(4,17,26,34,45)] <- "99"
mv[c(19,20,27,28,35,36,39,40,46:53)] <- "999"
mv[5] <- "9999"
mv[1] <- "999999"
for (i in seq(a = orig)) orig[[i]][orig[[i]] == mv[i]] <- NA
summary(orig[1:5])


###################################################
### code chunk number 5: dropunused
###################################################
for (i in seq(a = orig)) orig[[i]] <- orig[[i]][drop = TRUE]
summary(orig[1:5])


###################################################
### code chunk number 6: lcase
###################################################
names(orig) <- tolower(names(orig))


###################################################
### code chunk number 7: factorlevels
###################################################
orig$hdegk <- ordered(orig$hdegk, levels = 1:6,
                      labels = c("ASSOC","BS/BA","MS/MA/MEd","MA+","Ed.S","Ed.D/Ph.D"))
orig$hdeg1 <- ordered(orig$hdeg1, levels = 1:4,
                      labels = c("BS/BA","MS/MA/MEd","Ed.S","Ed.D/Ph.D"))
orig$hdeg2 <- ordered(orig$hdeg2, levels = 1:4,
                     labels = c("BS/BA","MS/MA/MEd","Ed.S","Ed.D/Ph.D"))
orig$hdeg3 <- ordered(orig$hdeg3, levels = 1:4,
                     labels = c("BS/BA","MS/MA/MEd","Ed.S","Ed.D/Ph.D"))
orig$cladk <- factor(orig$cladk, levels = c(1:3,5:8),
                     labels = c("1","2","3","APPR","PROB","NOT","PEND"))
orig$clad1 <- factor(orig$clad1, levels = 1:6,
                     labels = c("NOT","APPR","PROB","1","2","3"))
orig$clad2 <- factor(orig$clad2, levels = 1:6,
                     labels = c("NOT","APPR","PROB","1","2","3"))
orig$clad3 <- factor(orig$clad3, levels = 1:6,
                     labels = c("NOT","APPR","PROB","1","2","3"))


###################################################
### code chunk number 8: student
###################################################
student <- orig[1:5]
names(student) <- c("id", "sx", "eth", "birthq", "birthy")
levels(student$sx) <- c("M", "F")
levels(student$eth) <- c("W", "B", "A", "H", "I", "O")
student$birthy <- ordered(student$birthy)
student$birthq <- ordered(paste(student$birthy,student$birthq,sep=":"))
summary(student)


###################################################
### code chunk number 9: long
###################################################
long <- data.frame(id = rep(orig$newid, 4),
                   gr = ordered(rep(c("K", 1:3), each = nrow(orig)),
                                levels = c("K", 1:3)),
                   star = factor(unlist(orig[6:9])),
                   cltype = factor(unlist(orig[10:13])),
                   schtype = factor(unlist(orig[c(14,22,30,38)])),
                   hdeg = ordered(unlist(lapply(orig[c(15,24,32,43)],as.character)),
                                  levels = c("ASSOC","BS/BA","MS/MA/MEd","MA+","Ed.S","Ed.D/Ph.D")),
                   clad = factor(unlist(lapply(orig[c(16,25,33,44)],as.character)),
                                  levels = c("NOT","APPR","PROB","PEND","1","2","3")),
                   exp = unlist(orig[c(17,26,34,45)]),
                   trace = factor(unlist(orig[c(18,23,31,42)]), levels=1:6,
                                  labels=c("W", "B", "A", "H", "I", "O")),
                   read = unlist(orig[c(19,27,35,39)]),
                   math = unlist(orig[c(20,28,36,40)]),
                   ses = factor(unlist(orig[c(21,29,37,41)]),labels=c("F","N")),
                   sch = factor(unlist(orig[50:53])))


###################################################
### code chunk number 10: summarylong
###################################################
summary(long)


###################################################
### code chunk number 11: isnacheck
###################################################
with(long, all.equal(is.na(schtype), is.na(sch)))
with(long, all.equal(is.na(cltype), is.na(sch)))


###################################################
### code chunk number 12: longsubset
###################################################
long <- long[!is.na(long$sch),]


###################################################
### code chunk number 13: summlong
###################################################
summary(long[1:5])


###################################################
### code chunk number 14: dropstar
###################################################
long$star <- NULL


###################################################
### code chunk number 15: rownames
###################################################
rownames(long) <- paste(long$id, long$gr, sep = '/')


###################################################
### code chunk number 16: schooltable
###################################################
school <- unique(long[, c("sch", "schtype")])
length(levels(school$sch)) == nrow(school)
row.names(school) <- school$sch
school <- school[order(as.integer(as.character(school$sch))),]
long$schtype <- NULL
levels(school$schtype) <- c("inner","suburb","rural","urban")
levels(long$cltype) <- c("small", "reg", "reg+A")


###################################################
### code chunk number 17: merging
###################################################
star <- merge(merge(long, school, by = "sch"), student, by = "id")
star$time <- as.integer(star$gr) - 1


###################################################
### code chunk number 18: teacher
###################################################
teacher <- unique(star[, c("cltype", "trace", "exp", "clad", "gr",
                           "hdeg", "sch")])
teacher <- teacher[with(teacher, order(sch, gr, cltype, exp, hdeg,
                                       clad, trace)), ]


###################################################
### code chunk number 19: teacherlabels
###################################################
row.names(teacher) <- tch <- teacher$tch <- seq(nrow(teacher))
names(tch) <- do.call("paste", c(teacher[,1:7], list(sep=":")))
star$tch <- tch[do.call("paste", c(star[c("cltype", "trace", "exp",
                                          "clad", "gr", "hdeg", "sch")],
                                   list(sep = ":")))]


###################################################
### code chunk number 20: classsizes
###################################################
table(table(star$tch))
table(table(subset(star, cltype == "small")$tch))
table(table(subset(star, cltype == "reg")$tch))
table(table(subset(star, cltype == "reg+A")$tch))


###################################################
### code chunk number 21: initialModel
###################################################
library(lme4)
(mm1 <- lmer(math ~ gr + sx + eth + cltype + (1|id) + (1|sch), star))
(rm1 <- lmer(read ~ gr + sx + eth + cltype + (1|id) + (1|sch), star))


###################################################
### code chunk number 22: sessionInfo
###################################################
toLatex(sessionInfo())



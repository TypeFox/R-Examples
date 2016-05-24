### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/data.tex'

###################################################
### code chunk number 1: data.tex:10-11
###################################################
library(HH)


###################################################
### code chunk number 2: data.tex:14-19
###################################################
## the standard lattice color 2 is difficult for people with color deficient vision
data(col3x2)
## These colors look like a 3x2 color array when run through
## the vischeck simulator to see how they look for the three most
## common color vision deficiencies: Protanope, Deuteranope, Tritanope.


###################################################
### code chunk number 3: data.tex:422-431
###################################################
## hhcapture("meltcast.Rout", '
library(reshape2)
wide <- data.frame(Names=LETTERS[1:5], x=1:5, y=6:10)
wide
long <- melt(wide, id="Names")
long
wideagain <- dcast(Names ~ variable, value="value", data=long)
wideagain
## ')


###################################################
### code chunk number 4: data.tex:608-618
###################################################
## hhcapture("NA-0.Rout", '
AA <- read.table(text="
x y
1 2
3 NA
5 6
", header=TRUE)
AA
sapply(AA, class)
## ')


###################################################
### code chunk number 5: data.tex:620-632
###################################################
## hhcapture("NA-0BB.Rout", '
BB <- read.table(text="
x y
1 2
3 999
5 6
7 .
9 10
", header=TRUE, na.strings=c("999", "."))
BB
sapply(BB, class)
## ')


###################################################
### code chunk number 6: data.tex:634-647
###################################################
## hhcapture("NA-0CC.Rout", '
CC <- read.table(text="
x y
1 2
3 999
5 6
7 .
9 10
", header=TRUE)
CC
sapply(CC, class)
CC$y
## ')


###################################################
### code chunk number 7: data.tex:699-707
###################################################
## hhcapture("NA-1.Rout", '
abcd <- data.frame(x=c(1, 2, NA, 4, 5, 6, 7, 8),
                   y=c(6, 5, 8, NA, 10, 9, 12, 11),
                   ch=c(NA, "N", "O", "P", "Q", "R", "S", "T"),
                   stringsAsFactors=FALSE)
abcd
sapply(abcd, class)
## ')


###################################################
### code chunk number 8: data.tex:726-732
###################################################
## hhpdf("NA-2.pdf")
## hhcapture("NA-2.Rout", '
xyplot(y ~ x, data=abcd, labels=abcd$ch, panel=panel.text,
       col=col3x2[2:3], cex=2)
## ')
## hhdev.off()


###################################################
### code chunk number 9: data.tex:748-756
###################################################
## hhcapture("NA-3.Rout", '
3 + NA
sum(3, NA)
sum(3, NA, na.rm=TRUE)
abcd$x
mean(abcd$x)
mean(abcd$x, na.rm=TRUE)
## ')


###################################################
### code chunk number 10: data.tex:777-782
###################################################
## hhcapture("NA-4.Rout", '
a.lm <- lm(y ~ x, data=abcd)
summary(a.lm)
predict(a.lm)
## ')


###################################################
### code chunk number 11: data.tex:802-807
###################################################
## hhcapture("NA-5.Rout", '
b.lm <- lm(y ~ x, data=abcd, na.action=na.exclude)
summary(b.lm)
predict(b.lm)
## ')



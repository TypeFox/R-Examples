### R code from vignette source 'maSAE.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: maSAE.Rnw:58-59
###################################################
library('maSAE')


###################################################
### code chunk number 2: maSAE.Rnw:63-65
###################################################
data('s2')
data('s1')


###################################################
### code chunk number 3: maSAE.Rnw:70-77
###################################################
s2$phase1 <- s2$phase2 <- TRUE
s1$phase1 <- TRUE
s1$phase2 <- FALSE
eval(parse(text=(paste('s1$'
		       , setdiff(names(s2), names(s1))
		       , ' <- NA' , sep = ''))))
s12 <- rbind(s1, s2)


###################################################
### code chunk number 4: maSAE.Rnw:80-82
###################################################
saeO <- saObj(data = s12, f = y ~x1 + x2 + x3 | g
	      , s2 = 'phase2')


###################################################
### code chunk number 5: maSAE.Rnw:85-86
###################################################
predict(saeO)


###################################################
### code chunk number 6: maSAE.Rnw:94-101
###################################################
data('s0')
s0$x1 <- s0$x3 <- NULL
s0$phase1 <- s0$phase2 <- FALSE
eval(parse(text=(paste('s0$'
		       , setdiff(names(s12), names(s0))
		       , ' <- NA' , sep = ''))))
s012 <- rbind(s0, s12)


###################################################
### code chunk number 7: maSAE.Rnw:104-106
###################################################
predict(saObj(data = s012,  f = y ~x1 + x2 + x3 | g
	      , s2 = 'phase2', s1 = 'phase1'))


###################################################
### code chunk number 8: maSAE.Rnw:113-117
###################################################
tm1 <- as.data.frame(tapply(s012$x2 , s012$g, mean))
 names(tm1)[1] <- c('x2'); tm1$g <- row.names(tm1)
predict(saObj(data = s12, f = y ~x1 + x2 + x3 | g
	      , s2 = 'phase2', smallAreaMeans = tm1))


###################################################
### code chunk number 9: maSAE.Rnw:125-132
###################################################
preds <- paste('x',1:3, sep='')
tm <- as.data.frame(
    rbind(
        colMeans(subset(s12, g =='a')[, preds])
        , colMeans(subset(s12, g =='b')[, preds])
        )
    ); tm$g=c('a', 'b')


###################################################
### code chunk number 10: maSAE.Rnw:135-137
###################################################
predict(saObj(data = s12, f = y ~x1 + x2 + x3 | g
	      , s2 = 'phase2', smallAreaMeans = tm))


###################################################
### code chunk number 11: maSAE.Rnw:143-152
###################################################
source('Rao.R')
library(nlme)
dat <- subset(s2, ! is.na(s2$g))
dat <- dat[with(dat, order(g)), ]
aLmeObj <- lme(y ~x1 + x2 + x3
               , data = dat
               , random =  ~1|g)
foo <- new(Class = "sae", lmeObj = aLmeObj, domain.df = tm)
sae(foo)


###################################################
### code chunk number 12: maSAE.Rnw:164-165
###################################################
 grep('.*clust', capture.output(str(s1)), value = TRUE)


###################################################
### code chunk number 13: maSAE.Rnw:173-176 (eval = FALSE)
###################################################
## predict(saObj(data = s12, f = y ~x1 + x2 + x3 | g
## 	      , s2 = 'phase2'
## 	      , cluster = 'clustid'))


###################################################
### code chunk number 14: maSAE.Rnw:178-181 (eval = FALSE)
###################################################
## predict(saObj(data = s012,  f = y ~x1 + x2 + x3 | g
## 	      , s2 = 'phase2', s1 = 'phase1'
## 	      , cluster = 'clustid'))


###################################################
### code chunk number 15: maSAE.Rnw:183-186 (eval = FALSE)
###################################################
## predict(saObj(data = s12, f = y ~x1 + x2 + x3 | g
## 	      , s2 = 'phase2', smallAreaMeans = tm1
## 	      , cluster = 'clustid'))


###################################################
### code chunk number 16: maSAE.Rnw:188-191 (eval = FALSE)
###################################################
## predict(saObj(data = s12, f = y ~x1 + x2 + x3 | g
## 	      , s2 = 'phase2', smallAreaMeans = tm
## 	      , cluster = 'clustid'))



### R code from vignette source 'partykit.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width = 70)
library("partykit")
set.seed(290875)


###################################################
### code chunk number 2: weather-data
###################################################
data("WeatherPlay", package = "partykit")
WeatherPlay


###################################################
### code chunk number 3: weather-plot0
###################################################
py <- party(
  partynode(1L,
    split = partysplit(1L, index = 1:3),
    kids = list(
      partynode(2L,
        split = partysplit(3L, breaks = 75),
        kids = list(
          partynode(3L, info = "yes"),
          partynode(4L, info = "no"))),
      partynode(5L, info = "yes"),
      partynode(6L,
        split = partysplit(4L, index = 1:2),
        kids = list(
          partynode(7L, info = "yes"),
          partynode(8L, info = "no"))))),
  WeatherPlay)
plot(py)


###################################################
### code chunk number 4: weather-splits
###################################################
sp_o <- partysplit(1L, index = 1:3)
sp_h <- partysplit(3L, breaks = 75)
sp_w <- partysplit(4L, index = 1:2)


###################################################
### code chunk number 5: weather-nodes
###################################################
pn <- partynode(1L, split = sp_o, kids = list(
  partynode(2L, split = sp_h, kids = list(
    partynode(3L, info = "yes"),
    partynode(4L, info = "no"))),
  partynode(5L, info = "yes"),
  partynode(6L, split = sp_w, kids = list(
    partynode(7L, info = "yes"),
    partynode(8L, info = "no")))))


###################################################
### code chunk number 6: weather-nodes-print
###################################################
pn


###################################################
### code chunk number 7: weather-party
###################################################
py <- party(pn, WeatherPlay)
print(py)


###################################################
### code chunk number 8: weather-plot (eval = FALSE)
###################################################
## plot(py)


###################################################
### code chunk number 9: weather-predict
###################################################
predict(py, head(WeatherPlay))


###################################################
### code chunk number 10: weather-methods-dim
###################################################
length(py)
width(py)
depth(py)


###################################################
### code chunk number 11: weather-methods-subset
###################################################
py[6]


###################################################
### code chunk number 12: weather-methods-names
###################################################
py2 <- py
names(py2)
names(py2) <- LETTERS[1:8]
py2


###################################################
### code chunk number 13: weather-methods-nodeids
###################################################
nodeids(py)
nodeids(py, terminal = TRUE)


###################################################
### code chunk number 14: weather-methods-nodeapply
###################################################
nodeapply(py, ids = c(1, 7), FUN = function(n) n$info)
nodeapply(py, ids = nodeids(py, terminal = TRUE),
  FUN = function(n) paste("Play decision:", n$info))


###################################################
### code chunk number 15: weather-methods-predict
###################################################
predict(py, FUN = function(n) paste("Play decision:", n$info))


###################################################
### code chunk number 16: weather-methods-print
###################################################
print(py, terminal_panel = function(n)
  c(", then the play decision is:", toupper(n$info)))


###################################################
### code chunk number 17: weather-methods-plot (eval = FALSE)
###################################################
## plot(py, tp_args = list(FUN = function(i) 
##   c("Play decision:", toupper(i))))


###################################################
### code chunk number 18: weather-methods-plot1
###################################################
plot(py[6])


###################################################
### code chunk number 19: weather-methods-plot2
###################################################
plot(py, tp_args = list(FUN = function(i) 
  c("Play decision:", toupper(i))))


###################################################
### code chunk number 20: weather-prune
###################################################
nodeprune(py, 2)
nodeprune(py, c(2, 6))


###################################################
### code chunk number 21: partysplit-1
###################################################
sp_h <- partysplit(3L, breaks = 75)
class(sp_h)


###################################################
### code chunk number 22: partysplit-2
###################################################
unclass(sp_h)


###################################################
### code chunk number 23: partysplit-3
###################################################
character_split(sp_h, data = WeatherPlay)


###################################################
### code chunk number 24: partysplit-4
###################################################
kidids_split(sp_h, data = WeatherPlay)


###################################################
### code chunk number 25: partysplit-5
###################################################
as.numeric(!(WeatherPlay$humidity <= 75)) + 1


###################################################
### code chunk number 26: partysplit-6
###################################################
sp_o2 <- partysplit(1L, index = c(1L, 1L, 2L))
character_split(sp_o2, data = WeatherPlay)
table(kidids_split(sp_o2, data = WeatherPlay), WeatherPlay$outlook)


###################################################
### code chunk number 27: partysplit-6
###################################################
unclass(sp_o2)


###################################################
### code chunk number 28: partysplit-7
###################################################
sp_o <- partysplit(1L, index = 1L:3L)
character_split(sp_o, data = WeatherPlay)


###################################################
### code chunk number 29: partysplit-8
###################################################
sp_t <- partysplit(2L, breaks = c(69.5, 78.8), index = c(1L, 2L, 1L))
character_split(sp_t, data = WeatherPlay)
table(kidids_split(sp_t, data = WeatherPlay),
  cut(WeatherPlay$temperature, breaks = c(-Inf, 69.5, 78.8, Inf)))


###################################################
### code chunk number 30: partynode-1
###################################################
n1 <- partynode(id = 1L)
is.terminal(n1)
print(n1)


###################################################
### code chunk number 31: partynode-2
###################################################
n1 <- partynode(id = 1L, split = sp_o, kids = lapply(2L:4L, partynode))
print(n1, data = WeatherPlay)


###################################################
### code chunk number 32: partynode-3
###################################################
fitted_node(n1, data = WeatherPlay)


###################################################
### code chunk number 33: partynode-4
###################################################
kidids_node(n1, data = WeatherPlay)


###################################################
### code chunk number 34: party-1a
###################################################
t1 <- party(n1, data = WeatherPlay)
t1


###################################################
### code chunk number 35: party-1b
###################################################
party(n1, data = WeatherPlay[0, ])


###################################################
### code chunk number 36: party-2
###################################################
t2 <- party(n1, 
  data = WeatherPlay,
  fitted = data.frame(
    "(fitted)" = fitted_node(n1, data = WeatherPlay),
    "(response)" = WeatherPlay$play,
    check.names = FALSE),
  terms = terms(play ~ ., data = WeatherPlay),
)


###################################################
### code chunk number 37: party-3
###################################################
t2 <- as.constparty(t2)
t2


###################################################
### code chunk number 38: constparty-plot
###################################################
plot(t2, tnex = 1.5)


###################################################
### code chunk number 39: party-4
###################################################
nd <- data.frame(outlook = factor(c("overcast", "sunny"),
  levels = levels(WeatherPlay$outlook)))
predict(t2, newdata = nd, type = "response")
predict(t2, newdata = nd, type = "prob")
predict(t2, newdata = nd, type = "node")



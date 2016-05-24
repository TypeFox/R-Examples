library("mlogit")
data("Electricity", package = "mlogit")
Electr <- mlogit.data(Electricity, id="id", choice="choice", 
                      varying=3:26, shape="wide", sep="")

Elec.mxl <- mlogit(choice~pf+cl+loc+wk+tod+seas|0, Electr, 
              rpar=c(pf='n', cl='n', loc='n', wk='n', tod='n', seas='n'), 
              R=100, halton=NA, print.level=0, panel=TRUE)

Elec.mxl2 <- mlogit(choice~pf+cl+loc+wk+tod+seas|0, Electr, 
               rpar=c(cl='n', loc='n', wk='n', tod='n', seas='n'), 
               R=100, halton=NA, print.level=0, panel=TRUE)

Elec.mxl3 <- update(Elec.mxl, rpar=c(cl='n', loc='n', wk='u', tod='n', seas='n'))

Electr <- mlogit.data(Electricity, id="id", choice="choice", 
                      varying=3:26, shape="wide", sep="",
                      opposite=c('tod', 'seas'))
Elec.mxl4 <- mlogit(choice~pf+cl+loc+wk+tod+seas|0, Electr, 
              rpar=c(cl='n', loc='n', wk='u', tod='ln', seas='ln'), 
              R=100, halton=NA, print.level=0, panel=TRUE)

Elec.mxl5 <- update(Elec.mxl4, correlation = TRUE)


data("Train", package = "mlogit")
Tr <- mlogit.data(Train, shape = "wide", varying = 4:11, 
                  choice = "choice", sep = "", 
                  opposite = c("price", "time", "change", "comfort"),
                  alt.levels=c("choice1", "choice2"), id="id")
Train.ml <- mlogit(choice ~ price + time + change + comfort, Tr)
Train.mxlc <- mlogit(choice ~ price + time + change + comfort, Tr,
               panel = TRUE, rpar = c(time = "n", change = "n", comfort = "n"),
               correlation = TRUE, R = 100, halton = NA)
Train.mxlu <- update(Train.mxlc, correlation = FALSE)


data("Fishing", package = "mlogit")
Fish <- mlogit.data(Fishing, shape="wide", varying=2:9, choice="mode")
Fish.mprobit <- mlogit(mode~price | income | catch, Fish, probit = TRUE, alt.subset=c('beach', 'boat','pier'))

data("Mode", package="mlogit")
Mo <- mlogit.data(Mode, choice='choice', shape='wide', varying=c(2:9))
p1 <- mlogit(choice~cost+time, Mo, seed = 20, R = 100, probit = TRUE)
p2 <- mlogit(choice~cost+time, Mo, seed = 21, R = 100, probit = TRUE)


save(Elec.mxl, Elec.mxl2, Elec.mxl3, Elec.mxl4, Elec.mxl5,
     Train.ml, Train.mxlc, Train.mxlu,
     Fish.mprobit, p1, p2,
     file='../../data/Examples.rda')



## ----set-options,echo=FALSE---------------------------------------------------
options(width=80)

## ---- warning=FALSE-----------------------------------------------------------
library(EMbC)

## ---- fig.width=5, fig.height=3.5, fig.align='center'-------------------------
plot(x2d@D,col=x2d@L)

## -----------------------------------------------------------------------------
mybc <- embc(x2d@D)

## -----------------------------------------------------------------------------
slotNames(mybc)

## ---- fig.width=3.8, fig.height=3, fig.show='hold'----------------------------
# the lkhp() function allows an offset parameter;
lkhp(mybc)      # left panel
lkhp(mybc,10)   # right panel

## -----------------------------------------------------------------------------
stts(mybc)

## -----------------------------------------------------------------------------
mybc@P[[1]]

## -----------------------------------------------------------------------------
mybc@R

## ---- fig.width=5, fig.height=3.5, fig.align='center'-------------------------
sctr(mybc)

## ---- fig.width=5, fig.height=6, fig.align='center'---------------------------
sctr(mybc,x2d@L)

## -----------------------------------------------------------------------------
cnfm(mybc,x2d@L)

## -----------------------------------------------------------------------------
head(expth)

## -----------------------------------------------------------------------------
# info=-1 supresses any step wise output information
mybcp <- stbc(expth, info=-1)

## -----------------------------------------------------------------------------
slotNames(mybcp)

## -----------------------------------------------------------------------------
stts(mybcp)

## ---- fig.width=5, fig.height=3.5, fig.align='center'-------------------------
sctr(mybcp)

## -----------------------------------------------------------------------------
# the expert labelling given in expth$lbl is used by default
cnfm(mybcp)

## ---- fig.width=5, fig.height=5, fig.align='center'---------------------------
# lims=c(a,b) limits the plot to a chunk of the trajectory
lblp(mybcp,lims=c(100,500))

## ---- fig.width=6, fig.height=4.5, fig.align='center'-------------------------
# this function allows a parameter lims=c(a,b) as well
view(mybcp,lims=c(100,500))

## ---- eval=FALSE--------------------------------------------------------------
#  # point-wise kml doc generation
#  # by setting display=TRUE we can automatically launch google-earth from within R
#  pkml(bc,display=TRUE)

## ---- fig.width=5, fig.height=4.2, fig.align='center'-------------------------
# time-spans, distances and heading directions;
varp(mybcp)

## ---- fig.width=5, fig.height=3.5, fig.align='center'-------------------------
# clustering input data (estimated local values of velocity and turn)
varp(mybcp@X)

## ---- fig.width=5, fig.height=3.5, fig.align='center'-------------------------
# certainties associated to each data-point
varp(mybcp@U)

## -----------------------------------------------------------------------------
# create a move object from a Movebank csv file
library(move)
moveObj <- move(system.file("extdata","leroy.csv.gz",package="move"))

## -----------------------------------------------------------------------------
# the moveObj is passed directly to the constructor
mybc3 <- stbc(moveObj,scv='height',info=-1)

## ---- fig.width=3.8, fig.height=3, fig.show='hold'----------------------------
# the lkhp() function allows an offset parameter;
lkhp(mybc3)      # left panel
lkhp(mybc3,50)   # right panel

## -----------------------------------------------------------------------------
stts(mybc3)

## ---- fig.width=6, fig.height=3.5, fig.align='center'-------------------------
# showVars=c(1,2,3) is the default option and it is only shown for illustrative purposes
sctr(mybc3,showVars=c(1,2,3))

## ---- eval=FALSE--------------------------------------------------------------
#  # with showClst=c() we can specify a particular subset of clusters to be shown
#  sct3(mybc3,showClst=c(5,6,7,8))

## -----------------------------------------------------------------------------
# dlta is the maximum likelihood difference to accept a relabelling
# dlta=1 (accept all changes) is the default behaviour
postbc3 <- smth(mybc3,dlta=0.9)

## -----------------------------------------------------------------------------
# smth sets the smoothing time window length in hours
prebc3 <- stbc(moveObj,smth=1,scv='height',info=-1)

## ---- fig.width=6, fig.height=3.5, fig.align='center'-------------------------
# regardless of a pre-smoothing, we can still aply a post-smoothing;
# there is no real need to instantiate the smoothed copy of prebc3;
# this is useful for saving memory in case of long trajectories;
lblp(postbc3,smth(prebc3),lims=c(200,600))

## ---- eval=FALSE--------------------------------------------------------------
#  pkml(smth(prebc3),showClst=6,display=TRUE)

## -----------------------------------------------------------------------------
# the relabelling can be reversed with reset=TRUE
rlbl(prebc3,6,5)

## ---- fig.width=6, fig.height=3.5, fig.align='center'-------------------------
# the solar height is the control variable used by default
# note the previous relabelling in *prebc3*, changing HLH labels to HLL.
chkp(smth(prebc3),lims=c(200,600))

## ---- fig.width=3.2, fig.height=3.3, fig.show='hold'--------------------------
tmp <- runif(nrow(expth))
expth1 <- expth[which(tmp<=0.5),]
expth2 <- expth[which(tmp>=0.5),]
view(expth1)
view(expth2)

## -----------------------------------------------------------------------------
# we can combine data.fame trajectories and move objects
# only for illustrative purposes !!!
mystckbc <- stbc(list(expth1,expth2,moveObj),info=-1)

## -----------------------------------------------------------------------------
stts(mystckbc)

## ---- fig.width=5, fig.height=3.5, fig.align='center'-------------------------
sctr(mystckbc)

## -----------------------------------------------------------------------------
# this will only work when expert labelling is given for all trajectories in the stack
cnfm(mystckbc)

## -----------------------------------------------------------------------------
slotNames(mystckbc)

## -----------------------------------------------------------------------------
class(mystckbc@bC)

## -----------------------------------------------------------------------------
class(mystckbc@bCS)

## -----------------------------------------------------------------------------
lapply(mystckbc@bCS,class)

## -----------------------------------------------------------------------------
bcInd1 <- slct(mystckbc,1)

## ---- fig.width=3.7, fig.height=3.2, fig.show='hold'--------------------------
sctr(slct(mystckbc,1))  # left panel
sctr(slct(mystckbc,3))  # right panel

## -----------------------------------------------------------------------------
# comparing individual 1 with its correspondant out of the population
cnfm(stbc(expth1,info=-1),slct(mystckbc,1))

## ---- fig.width=6, fig.height=3.5, fig.align='center'-------------------------
# comparing inividuals 1 and 2 within the population
lblp(slct(mystckbc,1),slct(mystckbc,2))


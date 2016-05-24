

source("main.R") 
load("../simData_example.Rdata") #loads example data set called simData


D <- tsdata$event
T <- tsdata$trt
Y1 <- tsdata$Y1
Y2 <- tsdata$Y2

#trtsel objects
trtsel.Y1 <- TrtSel(disease = D, trt = T, marker = Y1, cohort.type="randomized cohort")
trtsel.Y1


trtsel.Y2 <- TrtSel(disease = D, trt = T, marker = Y2, study.design="randomized cohort")
trtsel.Y2

#plot
tmp <- plot(trtsel.Y1, plot.type = "cdf", bootstraps = 50)
head(tmp)

plot(trtsel.Y1, bootstraps = 200, ci = "vertical", plot.type = "treatment effect")

plot(trtsel.Y1, plot.type = "cdf", conf.bands = FALSE)#, fixed.values = seq(from=0.01, to=.4, by=.01))
plot(trtsel.Y1, plot.type = "risk", ylim = c(0, .8), main = "NEW MAIN HERE", bootstraps = 100 )
plot(trtsel.Y1, plot.type = "risk", ci = "horizontal" , fixed.values = c(.2, .25, .3, .35, .4))
plot(trtsel.Y2, bootstraps = 500)


plot(trtsel.Y2)


#eval
eval.Y1 <- evalTrtSel(trtsel.Y1, bootstraps= 100)
eval.Y1


eval.Y2 <- evalTrtSel(trtsel.Y2, bootstraps = 0)
eval.Y2

#compare
mycompare <- compare(trtsel1 = trtsel.Y1, trtsel2 = trtsel.Y2, bootstraps = 100)
mycompare

tmp <- plot(mycompare, bootstraps = 100)


#calibrate
cali.coh.Y1 <- calibrate(trtsel.Y1, plot.type = "risk.t0")



cali.coh.Y2 <- calibrate(trtsel.Y2)




##### BELOW is not functional anymore, I have been using it to check the code

### different sample designs

source("../trtsel_Aug2013/sim_functions.R")
load("../trtsel_Aug2013/my_sim_FY.Rdata")
alpha.strong <- c( -1.2402598, -0.6910426,  0.6, -2.25) #need to provide this for the bounded marker
#y.strong <- seq( -15, 15, by = .01)

n = 50000
simData <- sim.data(n=n, d.vec = d.vec,
                grid.y = grid.y, FY.11 = FY.11, FY.10 = FY.10, FY.01 = FY.01, FY.00 = FY.00)



nmatch = 1

D <- simData$D
T <- simData$T
Y1 <- simData$Y1
Y2 <- simData$Y2

# generate case-control subset (sample based on D only)
S <- NULL
S[D==1] <- 1 #select all cases
numcontrols <- length(D[D==1])*nmatch
S[D==0] <- sample(c(rep(1,numcontrols),rep(0,length(D[D==0])-numcontrols)))

myD<-D[S==1]; myT<-T[S==1]; myY<-Y2[S==1]

my.trtsel<-trtsel(event="D",trt="T",marker="Y2", data = simData,
                  default.trt = "trt none")

cc.trtsel<-trtsel(event="D",trt="T",marker="Y2", data = simData[S==1,],
                  cohort.attributes = c(n, mean(T), mean(D), 1),
                  study.design="nested case control", 
                  default.trt = "trt none")
 

#rho = c(mean(D), 1000000, mean(D[T==0]), mean(D[T==1]), nmatch, sum(T==1),0)

plot(cc.trtsel, bootstraps=500, plot.type = "risk", trt.names = c("marshall", "brownsworth"))


mean(1-myT[myD==1])

##STRATIFIED CASE CONTROL

nmatch = 1
# generate case-control subset (sample based on R and T)
S <- NULL

S[D==1] <- 1 #select all cases
numcontrols <- length(D[D==1 & T==0])*nmatch
#numcontrols <- sum(myconts.t0)*nmatch
S[D==0 & T==0] <- sample(c(rep(1,numcontrols),rep(0,length(D[D==0 & T==0])-numcontrols)))

#numcontrols <- sum(myconts.t0)*nmatch
numcontrols <- length(D[D==1 & T==1])*nmatch
S[D==0 & T==1] <- sample(c(rep(1,numcontrols),rep(0,length(D[D==0 & T==1])-numcontrols)))

# fit risk model

myD<-D[S==1]; myT<-T[S==1]; myY<-Y2[S==1]
#rho[1] = Pr(D = 1 | T = 0)
#rho[2] = Pr(D = 1 | T = 1)
#  N.t0.r0 <- rho[3]
#  N.t1.r0 <- rho[4]
#  N.t1    <- rho[5]
#  N       <- rho[6]   
scc.trtsel<-trtsel(event="D",trt="T",marker="Y2", data = simData[S==1,],
                   cohort.attributes = c(n, mean(D==0 & T==0), mean(D==1 & T==0), mean(D==0 & T==1), 1,1),
                  study.design="stratified nested case control", 
                  default.trt = "trt none")

coh <- eval.trtsel(my.trtsel, bootstraps = 0)#500)

cc <- eval.trtsel(cc.trtsel, bootstraps = 0)#500)

scc <- eval.trtsel(scc.trtsel, bootstraps = 0)#500)


rbind(coh$estimates,  
      cc$estimates, 
      scc$estimates)

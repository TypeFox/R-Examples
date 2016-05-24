## ----options,echo=FALSE--------------------------------------------------
# colorblind friendly palette
palette(c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
library('knitr')
knitr::opts_chunk$set(fig.width=7,fig.height=4.7,fig.show='hold')

## ----loading-------------------------------------------------------------
suppressMessages(library('nonlinearTseries'))
library('plot3D')
# by default, the simulation creates a RGL plot of the system's phase space
lor = lorenz(do.plot = F)
# let's plot the phase space of the simulated lorenz system
scatter3D(lor$x, lor$y, lor$z,
          main = "Lorenz's system phase space",
          col = 1, type="o",cex = 0.3)

## ----tauEstimation-------------------------------------------------------
# suppose that we have only measured the x-component of the Lorenz system
lor.x = lor$x

old.par = par(mfrow = c(1, 2))
# tau-delay estimation based on the autocorrelation function
tau.acf = timeLag(lor.x, technique = "acf",
                  lag.max = 100, do.plot = T)
# tau-delay estimation based on the mutual information function
tau.ami = timeLag(lor.x, technique = "ami", 
                  lag.max = 100, do.plot = T)
par(old.par)

## ----mEstimation---------------------------------------------------------
emb.dim = estimateEmbeddingDim(lor.x, time.lag = tau.ami,
                               max.embedding.dim = 15)


## ----buildTakens---------------------------------------------------------
tak = buildTakens(lor.x,embedding.dim = emb.dim, time.lag = tau.ami)
scatter3D(tak[,1], tak[,2], tak[,3],
          main = "Lorenz's system reconstructed phase space",
          col = 1, type="o",cex = 0.3)

## ----corrDim-------------------------------------------------------------
cd = corrDim(lor.x,
             min.embedding.dim = emb.dim,
             max.embedding.dim = emb.dim + 5,
             time.lag = tau.ami, 
             min.radius = 0.001, max.radius = 50,
             n.points.radius = 40,
             do.plot=FALSE)
plot(cd)
cd.est = estimate(cd, regression.range=c(0.75,3),
                  use.embeddings = 5:7)
cat("expected: 2.05  --- estimate: ",cd.est,"\n")

## ----sampEnt-------------------------------------------------------------
se = sampleEntropy(cd, do.plot = F)
se.est = estimate(se, do.plot = F,
                  regression.range = c(8,15))
cat("Sample entropy estimate: ", mean(se.est), "\n")

## ----maxLyap,fig.show='hide'---------------------------------------------
# get the sampling period of the lorenz simulation
# computing the differences of time (all differences should be equal)
sampling.period = diff(lor$time)[1]
ml = maxLyapunov(lor.x, 
                 sampling.period=0.01,
                 min.embedding.dim = emb.dim,
                 max.embedding.dim = emb.dim + 3,
                 time.lag = tau.ami, 
                 radius=1,
                 max.time.steps=1000,
                 do.plot=FALSE)
plot(ml,type="l", xlim = c(0,8))
ml.est = estimate(ml, regression.range = c(0,3),
                  do.plot = T,type="l")
cat("expected: 0.906  --- estimate: ", ml.est,"\n")

## ----surrogateTest-------------------------------------------------------
st = surrogateTest(lor.x,significance = 0.05,one.sided = F,
                   FUN = timeAsymmetry, do.plot=F)
plot(st)


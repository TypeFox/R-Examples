## ----init, message=FALSE-------------------------------------------------

require(peakPick)


## ----construct spikes----------------------------------------------------

set.seed(123)
spikes <- matrix(runif(200), ncol=2)
spikes[40, 1] <- 300
spikes[50, 2] <- 40


## ----detect spikes-------------------------------------------------------

spikehits <- detect.spikes(spikes, c(11, 90), 10, spike.min.sd=12)

## Spikes are correctly recovered
which(spikehits[, 1])
which(spikehits[, 2])


## ----construct peaks-----------------------------------------------------

set.seed(123)
peaks <- 1 + 0.5*sin((1:1000)*pi/100) + 0.1*rnorm(1000)
smoothpeaks <- filter(peaks, rep(1/20, 20), sides=2)


## ----detect peaks, fig.show='hold'---------------------------------------

peakhits <- peakpick(matrix(smoothpeaks, ncol=1), 100, peak.npos=50)

## Plot results
plot(smoothpeaks)
points((1:1000)[peakhits], smoothpeaks[peakhits], col="red")



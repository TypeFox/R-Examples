# sample analysis of eye movement data using emov in R
# by Simon Schwab
library(calibrate) # textxy

setwd("~/Work/code/emov/pkg/R")
source("emov.R")
# install.packages("circular")
#library("circular")

# read raw data file
data = emov.read_iviewsamples(
  "~/Data/nscenes/natural_scenes_samples.txt", 46)

# handle missing data: Iview has 0 for missing data
data$L.Raw.X..px.[data$L.Raw.X..px. == 0] = NA
data$L.Raw.Y..px.[data$L.Raw.Y..px. == 0] = NA
data$L.POR.X..mm.[data$L.POR.X..mm. == 0] = NA
data$L.POR.Y..mm.[data$L.POR.Y..mm. == 0] = NA
data$L.GVEC.X[data$L.GVEC.X == 0] = NA
data$L.GVEC.Y[data$L.GVEC.Y == 0] = NA
data$L.GVEC.Z[data$L.GVEC.Z == 0] = NA

# select channels to use
data = data.frame(t = data$Time, x = data$L.POR.X..mm, y = -data$L.POR.Y..mm)

# filter data, 1 deg is 14.4 mm, filtering > 750 deg/s (10800 mm/s)
flt = emov.filter(data$x, data$y, 10500/200)
data$x = flt$x
data$y = flt$y

# cart2sphere
#data = emov.cart2sphere(data$L.GVEC.X, data$L.GVEC.Y, data$L.GVEC.Z)
#data = data.frame(x=deg(data$az), y=-deg(data$elev))

# trial segmentation
n = 12 # number of trials
idx = c() # index
start = 1
for (i in 1:n) {
  idx = c(idx, start, start - 1 + 2000)
  start = start + 2000
}
idx <- matrix(idx, nrow=n, ncol=2, byrow=TRUE)
idx <- data.frame(start=idx[,1], end=idx[,2]) # easy to access

# fixation detecton for each trial
max_disp  = 19.0 # in cm, 28.8 cm (2 deg)
min_dur =  80/1000*200
fix = emov.idt(data$t, data$x, data$y, max_disp, min_dur)

# fixation segmentation for easy ploting
fixseg = list()
for (i in 1:n) {  
  start = data$t[idx[i,1]]
  end   = data$t[idx[i,2]]
  
  fixseg[[i]] = fix[fix$start >= start & fix$end <= end, ]
}

# Plot all trials, raw data and fixations
#my_xlim = c(-35, 25)
#my_ylim = c(-20, 15)
my_xlim = c(0, 770)
my_ylim = c(-680, -250)
c = sqrt(80/pi) # constant, r=1 corresponds to fixation duration of 50 ms.
par(mfcol=c(4,3))
for (i in 1:n) {
  plot(fixseg[[i]]$x, fixseg[[i]]$y,
       xlim=my_xlim, ylim=my_ylim,
       xlab=NA, ylab=NA, pch=19,
       cex=sqrt(fixseg[[i]][, 3] * 10^-3 * pi^-1) * c^-1, col='gray')
  textxy(fixseg[[i]]$x, fixseg[[i]]$y, 1:length(fixseg[[i]]$x), cx=1)
  par(new=TRUE)
  plot(fixseg[[i]]$x, fixseg[[i]]$y,
       xlim=my_xlim, ylim=my_ylim,
       xlab=NA, ylab=NA, cex=1)
  par(new=TRUE)
  plot(data$x[idx$start[i]:idx$end[i]],
       data$y[idx$start[i]:idx$end[i]],
       type="l", xlim=my_xlim,  ylim=my_ylim,
       xlab="Horizontal (px)", ylab="Vertical (px)")
}

# Plot single trial
par(mfcol=c(1,1))
nr = 1
plot(fixseg[[nr]]$x, fixseg[[nr]]$y, xlim=my_xlim,  ylim=my_ylim, 
     xlab=NA, ylab=NA, pch=19,
     cex=sqrt(fixseg[[nr]][, 3] * 10^-3 * pi^-1) * c^-1,
     col='gray')
textxy(fixseg[[nr]]$x, fixseg[[nr]]$y, 1:length(fixseg[[nr]]$x), cx=1)
par(new=TRUE)
plot(fixseg[[nr]]$x, fixseg[[nr]]$y, xlim=my_xlim,  ylim=my_ylim, 
     xlab=NA, ylab=NA, cex=1)
par(new=TRUE)
plot(data$x[idx$start[nr]:idx$end[nr]],
     data$y[idx$start[nr]:idx$end[nr]],  
     type="l", xlim=my_xlim,  ylim=my_ylim,
     xlab="Horizontal (px)", ylab="Vertical (px)")

# Plot stimuli
# install.packages("jpeg")
# library(jpeg)
# 
# img = list()
# img[[1]] <- readJPEG("/home/simon/Data/nscenes/stimuli/000.jpg")
# img[[2]] <- readJPEG("/home/simon/Data/nscenes/stimuli/001.jpg")
# img[[3]] <- readJPEG("/home/simon/Data/nscenes/stimuli/002.jpg")
# img[[4]] <- readJPEG("/home/simon/Data/nscenes/stimuli/003.jpg")
# img[[5]] <- readJPEG("/home/simon/Data/nscenes/stimuli/004.jpg")
# img[[6]] <- readJPEG("/home/simon/Data/nscenes/stimuli/005.jpg")
# img[[7]] <- readJPEG("/home/simon/Data/nscenes/stimuli/006.jpg")
# img[[8]] <- readJPEG("/home/simon/Data/nscenes/stimuli/007.jpg")
# img[[9]] <- readJPEG("/home/simon/Data/nscenes/stimuli/008.jpg")
# img[[10]] <- readJPEG("/home/simon/Data/nscenes/stimuli/009.jpg")
# img[[11]] <- readJPEG("/home/simon/Data/nscenes/stimuli/010.jpg")
# img[[12]] <- readJPEG("/home/simon/Data/nscenes/stimuli/011.jpg")
# 
# par(mfcol=c(4,3))
# 
# for (i in 1:n) {
#   plot(c(0,1), c(0,1))
#   rasterImage(img[[i]], 0, 0, 1, 1)
# }

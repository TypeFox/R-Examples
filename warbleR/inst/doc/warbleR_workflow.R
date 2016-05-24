## ---- eval=FALSE---------------------------------------------------------
#  
#  install.packages("warbleR")
#  library(warbleR)
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  # Create a new directory
#  dir.create(file.path(getwd(),"warbleR_example"))
#  setwd(file.path(getwd(),"warbleR_example"))
#  
#  # Check the location of the directory
#  getwd()
#  

## ---- eval=TRUE, echo=FALSE, message=FALSE-------------------------------

# this sets my working directory
library(warbleR)
library(knitr)


## ---- eval=TRUE----------------------------------------------------------

# Query Xeno-Canto for all recordings of the hummingbird genus Phaethornis
Phae <- querxc(qword = "Phaethornis", download = FALSE) 

# Find out what kind of metadata we have
names(Phae) 


## ---- eval=FALSE---------------------------------------------------------
#  
#  View(Phae)
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  # Query Xeno-Canto for all recordings of the species Phaethornis longirostris
#  Phae.lon <- querxc(qword = "Phaethornis longirostris", download = FALSE)
#  View(Phae.lon)
#  

## ---- eval=TRUE, echo=FALSE, message=FALSE-------------------------------

Phae.lon <- querxc(qword = "Phaethornis longirostris", download = FALSE) 


## ---- eval=FALSE---------------------------------------------------------
#  
#  # Image type default is jpeg, but tiff files have better resolution
#  xcmaps(X = Phae, img = TRUE, it = "tiff")
#  xcmaps(X = Phae.lon, img = FALSE)
#  

## ---- eval=TRUE, echo=FALSE, message=FALSE-------------------------------

xcmaps(X = Phae.lon, img = FALSE) 


## ---- eval=TRUE----------------------------------------------------------

# Find out number of available recordings
nrow(Phae.lon) 

# Find out how many types of signal descriptions exist in the Xeno-Canto metadata
levels(Phae.lon$Vocalization_type)

# How many recordings per signal type?
table(Phae.lon$Vocalization_type)


## ---- eval=TRUE----------------------------------------------------------

# There are many levels to the Vocalization_type variable. 
# Some are biologically relevant signals, but most just 
# reflect variation in data entry.

# Luckily, it's very easy to filter the signals we want 
Phae.lon.song <- droplevels(Phae.lon[grep("song", Phae.lon$Vocalization_type, 
                                ignore.case = TRUE), ])

# Check resulting data frame
str(Phae.lon.song) 


## ---- eval=FALSE---------------------------------------------------------
#  
#  # Now, how many recordings per locatity?
#  table(Phae.lon.song$Locality)
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  # In case you want more than one signal type you can try something like this:
#  Phae.lon.sc <- Phae.lon[grep("song|call", Phae.lon$Vocalization_type,ignore.case = TRUE), ]
#  str(Phae.lon.sc)
#  

## ---- eval=TRUE----------------------------------------------------------

#first filter by location
Phae.lon.LS <- Phae.lon.song[grep("La Selva Biological Station, Sarapiqui, Heredia", Phae.lon.song$Locality,
                              ignore.case = FALSE),]

# And only those of the highest quality
Phae.lon.LS <- Phae.lon.LS[Phae.lon.LS$Quality == "A", ]

## ---- eval=TRUE----------------------------------------------------------
# map in the graphic device (img = FALSE)
xcmaps(Phae.lon.LS, img = FALSE)


## ---- eval=FALSE, echo=FALSE---------------------------------------------
#  
#  # This copies the selected sound files to a dropbox folder so they can be shared
#  # do not show this code
#  fn <- with(Phae.lon.LS,paste(paste(Genus, Specific_epithet, Recording_ID, sep = "-"), ".wav", sep = " "))
#  file.copy(from = file.path("/home/m/Documents/Biblioteca de cantos/Trochilidae/XC/wavs",fn), to = file.path("/home/m/Dropbox/Projects/warbleR package/vignette files", fn), overwrite = TRUE)
#  
#  wlist <- lapply(fn,function(x) downsample(readWave(file.path("/home/m/Dropbox/Projects/warbleR package/vignette files", x)), samp.rate = 22050))
#  
#  names(wlist) <- fn
#  
#  saveRDS(wlist, file = "/home/m/Dropbox/Sharing/warbleR/recs.RDS")
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  # Download sound files
#  querxc(X = Phae.lon.LS)
#  
#  # Save each data frame object as a .csv file
#  write.csv(Phae.lon.LS, "Phae_lon.LS.csv", row.names = FALSE)
#  

## ---- eval=FALSE---------------------------------------------------------
#  # Neither of these functions requires arguments
#  # Always check you're in the right directory beforehand
#  # getwd()
#  mp32wav()
#  
#  # You can use checkwavs to see if wav files can be read
#  checkwavs()
#  
#  # Let's create a list of all the recordings in the directory
#  wavs <- list.files(pattern="wav$")
#  
#  # We will use this list to downsample the wav files so the following analyses go a bit faster
#  lapply(wavs, function(x) writeWave(downsample(readWave(x), samp.rate = 22050),
#                                    filename = x))

## ---- eval=FALSE, echo=FALSE---------------------------------------------
#  
#  ### If you were unable to convert _mp3_ to _wav_ format:
#    + download the file in [this link](https://www.dropbox.com/s/htpbxbdw8s4i23k/recs.RDS?dl=0) and put it in your working directory
#    + then run the following code:
#  
#  
#  recs <- readRDS(file = "recs.RDS")
#  
#  for(i in 1:length(recs))
#    writeWave(recs[[i]], filename = names(recs)[i])
#  
#  *Note: In case you have your own recordings in _wav_ format and have skipped previous sections, you must specify the location of your sound files prior to running downstream functions.*
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  # Let's first create a subset for playing with arguments
#  # This subset is based on the list of wav files we created above
#  sub <- wavs[c(1,3)]
#  
#  # ovlp = 10 speeds up process a bit
#  # tiff image files are better quality and are faster to produce
#  lspec(flist = sub, ovlp = 10, it = "tiff")
#  
#  # We can zoom in on the frequency axis by changing flim,
#  # the number of seconds per row, and number of rows
#  lspec(flist = sub, flim = c(1.5, 11), sxrow = 6, rows = 15, ovlp = 10, it = "tiff")
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  lspec(flim = c(1.5, 11), ovlp = 10, sxrow = 6, rows = 15, it = "tiff")
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  # List the image files in the directory
#  # Change the pattern to "jpeg" if you used that image type
#  imgs <- list.files(pattern = ".tiff")
#  
#  # If the maps we created previously are still there, you can remove them from this list easily
#  imgs <- imgs[grep("Map", imgs, invert = TRUE)]
#  
#  # Extract the recording IDs of the files for which image files remain
#  kept <- unique(sapply(imgs, function(x){
#    strsplit(x, split = "-", fixed = TRUE)[[1]][3]
#    }, USE.NAMES = FALSE))
#  
#  # Now we can get rid of sound files that do not have image files
#  snds <- list.files(pattern = ".wav", ignore.case = TRUE)
#  file.remove(snds[grep(paste(kept, collapse = "|"), snds, invert = TRUE)])
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  # Select a subset of the recordings
#  wavs <- list.files(pattern = ".wav", ignore.case = TRUE)
#  
#  # Set a seed so we all have the same results
#  set.seed(1)
#  sub <- wavs[sample(1:length(wavs), 3)]
#  
#  # Run autodetec() on subset of recordings
#  autodetec(flist = sub, bp = c(1, 10), threshold = 10, mindur = 0.05, maxdur = 0.5, envt="abs",
#            ssmooth = 300, ls = TRUE, res = 100,
#            flim = c(1, 12), wl = 300, set =TRUE, sxrow = 6, rows = 15,
#            redo = FALSE, it = "tiff")
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  autodetec(flist = sub, bp = c(2, 9), threshold = 20, mindur = 0.09, maxdur = 0.22,
#                       envt = "abs", ssmooth = 900, ls = TRUE, res = 100,
#                       flim= c(1, 12), wl = 300, set =TRUE, sxrow = 6, rows = 15,
#                       redo = TRUE, it = "tiff", img = TRUE, smadj = "end")

## ---- eval=FALSE---------------------------------------------------------
#  
#  Phae.ad <- autodetec(bp = c(2, 9), threshold = 20, mindur = 0.09, maxdur = 0.22,
#                       envt = "abs", ssmooth = 900, ls = TRUE, res = 100,
#                       flim= c(1, 12), wl = 300, set =TRUE, sxrow = 6, rows = 15,
#                       redo = TRUE, it = "tiff", img = TRUE, smadj = "end")
#  
#  str(Phae.ad)

## ---- eval=FALSE---------------------------------------------------------
#  
#  table(Phae.ad$sound.files)
#  

## ---- eval=TRUE, echo=FALSE----------------------------------------------

Phae.snr <- Phae.ad <- read.csv("Phae.snr.csv")


# Look at the number of selections per sound file 
table(Phae.ad$sound.files)


## ---- eval=FALSE---------------------------------------------------------
#  
#  # A margin that's too large causes other signals to be included in the noise measurement
#  # Re-initialize X as needed, for either autodetec or manualoc output
#  
#  # Let's try it on 10% of the selections so it goes a faster
#  # Set a seed first, so we all have the same results
#  set.seed(5)
#  
#  X <- Phae.ad[sample(1:nrow(Phae.ad),(nrow(Phae.ad)*0.1)), ]
#  
#  snrspecs(X = X, flim = c(2, 110), snrmar = 0.5, mar = 0.7, it = "tiff")
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  # This smaller margin is better
#  snrspecs(X = X, flim = c(2, 11), snrmar = 0.2, mar = 0.7, it = "tiff")
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  snrspecs(X = Phae.ad, flim = c(2, 11), snrmar = 0.2, mar = 0.7, it = "tiff")
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  Phae.snr <- sig2noise(X = Phae.ad[seq(1, nrow(Phae.ad), 2), ], mar = 0.04)
#  

## ---- eval=TRUE----------------------------------------------------------

Phae.hisnr <- Phae.snr[ave(-Phae.snr$SNR, Phae.snr$sound.files, FUN = rank) <= 5, ]

# Double check the number of selection per sound files 
table(Phae.hisnr$sound.files)


## ---- eval=FALSE---------------------------------------------------------
#  
#  write.csv(Phae.hisnr, "Phae_lon_autodetec_selecs.csv", row.names = FALSE)
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  # Run manualoc() with frequency range set for Phaethornis longirostris
#  # Recording comments are enabled to mark recording quality
#  # Selection comments enabled to include visual classifications
#  manualoc(flim = c(1, 11), reccomm = TRUE, selcomm = TRUE, osci = TRUE)
#  
#  # Read manualoc() output back into RStudio as an object
#  # This data frame object can be used as input for the functions that follow
#  manualoc_out <- read.csv("manualoc_output.csv", header = TRUE)
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  # Create a subset of 5 recordings analyzed by autodetec() or manualoc()
#  # Speeds up process of playing around with arguments
#  # Run either line below to reinitialize X with either autodetec
#  # or manualoc subset as desired
#  
#  set.seed(50)
#  X <- Phae.hisnr[sample(1:nrow(Phae.hisnr), 5), ]
#  
#  # Plot selection lines from manualoc() or autodetec()
#  specreator(X, osci = FALSE, line = TRUE, wl = 300, flim = c(1, 11), it = "tiff")
#  
#  # Change frequency limits of y-axis
#  specreator(X, flim = c(1, 11), osci = TRUE, line = TRUE, wl = 300, it = "tiff")
#  
#  # Change width of spectrogram to be proportional to signal duration
#  specreator(X, flim = c(1, 11), osci = TRUE, line = TRUE, propwidth = TRUE, wl = 300, it = "tiff")
#  
#  # Change spectrogram size
#  # Changing inner.mar and outer.mar arguments improves picsize results
#  specreator(X, flim = c(1, 11), osci = TRUE, line = TRUE, picsize = 1.5, wl = 300,
#             ovlp = 90, inner.mar = c(4, 4.5, 2, 1), outer.mar = c(4, 2, 2, 1), it = "tiff")
#  
#  # Run function for all recordings, with final argument settings
#  specreator(Phae.hisnr, flim = c(1, 11), osci = TRUE, line = TRUE, wl = 300,
#             ovlp = 90, it = "tiff", res = 300)

## ---- eval=FALSE---------------------------------------------------------
#  
#  # Note that the dominant frequency measurements are almost always more accurate
#  trackfreqs(Phae.hisnr, flim = c(1, 11), bp = c(1, 12), it = "tiff")
#  
#  # We can change the lower end of bandpass to make the frequency measurements more precise
#  trackfreqs(Phae.hisnr, flim = c(1, 11), bp = c(2, 12), col = c("purple", "orange"),
#             pch = c(17, 3), res = 300, it = "tiff")
#  
#  # If the frequency measurements look acceptable with this bandpass setting,
#  # that's the setting we should use when running specan()
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  # Use the bandpass filter to your advantage, to filter out low or high background
#  # noise before performing measurements
#  # The amplitude threshold will change the amplitude at which noises are
#  # detected for measurements
#  params <- specan(Phae.hisnr, bp = c(1, 11), threshold = 15)
#  
#  View(params)
#  
#  str(params)
#  
#  # As always, it's a good idea to write .csv files to your working directory

## ---- eval=TRUE, echo=FALSE----------------------------------------------
params <- read.csv("acoustic_parameters.csv")

str(params)


## ---- eval=TRUE----------------------------------------------------------

params <- params[, grep("fun|peakf", colnames(params), invert = TRUE)]


## ---- eval=TRUE, dpi=220-------------------------------------------------

# Run the PCA with only numeric variables of params
pca <- prcomp(x = params[, sapply(params, is.numeric)], scale. = TRUE)

# Check loadings
summary(pca)

# Extract PCA scores
pcascor <- as.data.frame(pca[[5]])

# Plot the 2 first PCs
plot(pcascor[, 1], pcascor[, 2], col = as.numeric(params$sound.files), pch = 20, 
     cex = 1, xlab = "PC1", ylab = "PC2")

# Add recordings/individuals labels 
x <- tapply(pcascor[, 1], params$sound.files, mean)
y <- tapply(pcascor[, 2], params$sound.files, mean)

labs <- gsub(".wav", "", unique(sapply(as.character(params$sound.files), function(x){
  strsplit(x, split = "-", fixed = TRUE)[[1]][3]
  }, USE.NAMES = FALSE)))

text(x, y, labs, cex=0.5)


## ---- eval=TRUE, dpi=220-------------------------------------------------

# Create a song type variable

# First, extract recording ID
songtype <- gsub(".wav", "", sapply(as.character(params$sound.files), function(x){
  strsplit(x, split = "-", fixed = TRUE)[[1]][3]
  }, USE.NAMES = FALSE))

# Now change IDs for letters representing song types
songtype <- gsub("154070|154072", "A", songtype)
songtype <- gsub("154129|154161", "B", songtype)
songtype <- gsub("154123", "C", songtype)
songtype <- gsub("154138", "D", songtype)

# Add song type as a variable representing symbol type
plot(pcascor[, 1], pcascor[, 2], col = as.numeric(params$sound.files), 
pch = as.numeric(as.factor(songtype)), 
     cex = 1, xlab = "PC1", ylab = "PC2")

# Add song type labels 
x <- tapply(pcascor[, 1], songtype, mean)
y <- tapply(pcascor[, 2], songtype, mean)

text(x, y, unique(songtype), cex = 1)



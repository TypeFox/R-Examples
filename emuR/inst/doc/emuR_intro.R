## ----results='hide', message=FALSE, warning=FALSE------------------------
# load the package
library(emuR)

# create demo data in folder provided by the tempdir() function
create_emuRdemoData(dir = tempdir())

# get the path to a .tpl file of a legacyEmuDB that is part of the demo data
tplPath = file.path(tempdir(), "emuR_demoData", "legacy_ae", "ae.tpl")

# convert this legacyEmuDB to the emuDB format
convert_legacyEmuDB(emuTplPath = tplPath, targetDir = tempdir())

## ------------------------------------------------------------------------
# remove the newly generated emuDB as we will not be needing it 
# throughout the rest of this vignette
unlink(file.path(tempdir(), "ae_emuDB"), recursive = TRUE)

## ----results='hide', message=FALSE, warning=FALSE------------------------
# get the path to a folder containing .wav & .TextGrid files that is part of the demo data
path2folder = file.path(tempdir(), "emuR_demoData", "TextGrid_collection")

# convert this TextGridCollection to the emuDB format
convert_TextGridCollection(path2folder, dbName = "myTGcolDB", 
                           targetDir = tempdir())

## ------------------------------------------------------------------------
# remove the newly generated emuDB as we will not be needing it 
# throughout the rest of this vignette
unlink(file.path(tempdir(), "myTGcolDB_emuDB"), recursive = TRUE)

## ------------------------------------------------------------------------
# get the path to a folder containing .wav & .par files that is part of the demo data
path2folder = file.path(tempdir(), "emuR_demoData", "BPF_collection")

# convert this BPFCollection to the emuDB format
convert_BPFCollection(path2folder, dbName = 'myBPF-DB', 
                      targetDir = tempdir(), verbose = F)

## ------------------------------------------------------------------------
# remove the newly generated emuDB as we will not be needing it 
# throughout the rest of this vignette
unlink(file.path(tempdir(), "myBPF-DB_emuDB"), recursive = TRUE)

## ------------------------------------------------------------------------
# get the path to emuDB called 'ae' that is part of the demo data
path2folder = file.path(tempdir(), "emuR_demoData", "ae_emuDB")

# load emuDB into current R session
ae = load_emuDB(path2folder, verbose = FALSE)

## ------------------------------------------------------------------------
summary(ae)

## ----eval=FALSE----------------------------------------------------------
#  serve(ae)

## ----message=FALSE, warning=FALSE----------------------------------------
sl = query(ae, query = "Phonetic==n")

head(sl)

## ------------------------------------------------------------------------
# calculate durations
d = dur(sl)

# calculate mean and by doing so answering the question
mean(d)

## ----results='hide', message=FALSE, warning=FALSE------------------------
# query emuDB
sl = query(ae, query = "Phonetic==I|o:|u:|V|@")

## ----results='hide', message=FALSE, warning=FALSE------------------------
# get formant values for those segments
td = get_trackdata(ae, sl,
                   onTheFlyFunctionName = "forest",
                   resultType = "emuRtrackdata")

## ------------------------------------------------------------------------
class(td)

## ---- fig.height = 5, fig.width = 5--------------------------------------
# load package
library(ggplot2)

# scatter plot of F1 and F2 values using ggplot
ggplot(td, aes(x=T2, y=T1, label=td$labels)) + 
  geom_text(aes(colour=factor(labels))) + 
  scale_y_reverse() + scale_x_reverse() + 
  labs(x = "F2(Hz)", y = "F1(Hz)") +
  guides(colour=FALSE)

## ------------------------------------------------------------------------
sibil = query(ae,"Phonetic==s|z|S|Z")

head(sibil)

## ------------------------------------------------------------------------
words = requery_hier(ae, sibil, level = "Word")

head(words)

## ------------------------------------------------------------------------
words = requery_hier(ae, sibil, level = "Text")

head(words)

## ------------------------------------------------------------------------
# get left context by off-setting the annotational units in sibil one unit to the left
leftContext = requery_seq(ae, sibil, offset = -1)

head(leftContext)

## ----eval=FALSE----------------------------------------------------------
#  # get right context by off-setting the annotational units in sibil one unit to the right
#  rightContext = requery_seq(ae, sibil, offset = 1)

## ------------------------------------------------------------------------
rightContext = requery_seq(ae, sibil, 
                           offset = 1, 
                           ignoreOutOfBounds = TRUE)

head(rightContext)

## ------------------------------------------------------------------------
sibil = query(ae,"Phonetic=~'[szSZ]'")

## ----results='hide', message=FALSE, warning=FALSE------------------------
dftTd = get_trackdata(ae, 
                      seglist = sibil,
                      onTheFlyFunctionName = 'dftSpectrum')

## ----results='hide'------------------------------------------------------
# execute this to show 16 spectra calculated from the first segment in sibil (an 's')
# (console output will not be shown here as it is very lengthy)
dftTd[1]

## ------------------------------------------------------------------------
dftTdRelFreq = dftTd[, 1000:10000]

## ------------------------------------------------------------------------
dftTdRelFreqMom = fapply(dftTdRelFreq, moments, minval = T)

## ------------------------------------------------------------------------
dftTdRelFreqMom[1]

## ----fig.height = 5, fig.width = 5---------------------------------------
dplot(dftTdRelFreqMom[, 1], 
      sibil$labels,
      normalise = TRUE, 
      xlab = "Normalized Time [%]", 
      ylab = "1st spectral moment [Hz]")

## ----fig.height = 5, fig.width = 5---------------------------------------
dplot(dftTdRelFreqMom[,1],
      sibil$labels,
      normalise = TRUE,
      average = TRUE,
      xlab = "Normalized Time [%]",
      ylab = "1st spectral moment [Hz]")

## ------------------------------------------------------------------------
# cut out the middle 60% portion
dftTdRelFreqMomMid = dcut(dftTdRelFreqMom, 
                          left.time = 0.2, 
                          right.time = 0.8, 
                          prop = T)

# display original moments of the first segment
dftTdRelFreqMom[1]

# display 60% portion moments of the first segment
dftTdRelFreqMomMid[1]

## ------------------------------------------------------------------------
meanFirstMoments = trapply(dftTdRelFreqMomMid[,1],
                           fun = mean,
                           simplify = T)

# display resulting vector
meanFirstMoments

## ----fig.height = 5, fig.width = 5---------------------------------------
boxplot(meanFirstMoments ~ sibil$labels)

## ----echo=FALSE, results='hide', message=FALSE, warning=FALSE------------
# disconnect to avoid file locking to sqliteDB that causes unlink
# to fail under windows
DBI::dbDisconnect(ae$connection)

## ----results='hide', message=FALSE, warning=FALSE------------------------
# remove emuR_demoData as we will not be needing it 
# throughout the rest of this vignette
unlink(file.path(tempdir(), "emuR_demoData"), recursive = TRUE)


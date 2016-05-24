## ----setup, include = FALSE, cache = FALSE---------------------
options(useFancyQuotes = FALSE)
library(knitr)
opts_chunk$set(comment = "#", out.height = "4.0in", out.width = "4.0in", fig.height = 5.2, fig.width = 5.2, fig.align = "center", cache = TRUE)
options(tidy = TRUE, width = 65)

## ----echo = FALSE, warning = FALSE, error = FALSE, message = FALSE----
library(tuneR)
library(monitoR)
#sf <- list.files('/home/sasha/Dropbox/UVMacoustic/monitoR/R', full.names = TRUE)
#sf <- sf[!grepl("eventEval.R", sf)]
#for(ff in sf) source(ff)
#df <- list.files('/home/sasha/Dropbox/UVMacoustic/monitoR/data', full.names = TRUE)
#for(ff in df) load(ff)

## ----eval = FALSE----------------------------------------------
#  library(monitoR)

## --------------------------------------------------------------
data(survey)
survey

## ----out.width = "5.8in", fig.width = 7------------------------
viewSpec(survey)

## ----eval = FALSE----------------------------------------------
#  setWavPlayer("play")
#  play(survey)

## --------------------------------------------------------------
data(btnw)
data(oven)
btnw
oven

## --------------------------------------------------------------
viewSpec(btnw)

## --------------------------------------------------------------
viewSpec(oven)

## --------------------------------------------------------------
btnw.fp <- file.path(tempdir(), "btnw.wav")
oven.fp <- file.path(tempdir(), "oven.wav")
survey.fp <- file.path(tempdir(), "survey2010-12-31_120000_EST.wav")
writeWave(btnw, btnw.fp)
writeWave(oven, oven.fp)
writeWave(survey, survey.fp)

## --------------------------------------------------------------
wct1 <- makeCorTemplate(btnw.fp)

## --------------------------------------------------------------
wct1

## ----fig.keep = "none"-----------------------------------------
wct1 <- makeCorTemplate(btnw.fp, name = "w1")
wct1

## --------------------------------------------------------------
wct2 <- makeCorTemplate(btnw.fp, t.lim = c(1.5, 2.1), frq.lim = c(4.2, 5.6), name = "w2")
wct2

## ----fig.keep = 'none'-----------------------------------------
oct1 <- makeCorTemplate(oven.fp, t.lim = c(1, 4), frq.lim = c(1, 11), name = "o1")

## ----fig.keep = 'none'-----------------------------------------
oct2 <- makeCorTemplate(oven.fp, t.lim = c(1, 4), frq.lim = c(1, 11), dens = 0.1, name = "o2")

## --------------------------------------------------------------
ctemps <- combineCorTemplates(wct1, wct2, oct1, oct2)
ctemps

## ----eval = FALSE----------------------------------------------
#  plot(ctemps)

## --------------------------------------------------------------
cscores <- corMatch(survey.fp, ctemps)

## --------------------------------------------------------------
cscores

## --------------------------------------------------------------
cdetects <- findPeaks(cscores)

## --------------------------------------------------------------
cdetects

## --------------------------------------------------------------
getDetections(cdetects)

## --------------------------------------------------------------
ctemps

## ----out.width = "5.8in", fig.width = 7------------------------
plot(cdetects)

## --------------------------------------------------------------
templateCutoff(ctemps)

## --------------------------------------------------------------
templateCutoff(ctemps)[2:4] <- c(0.3, 0.2, 0.2)

## --------------------------------------------------------------
templateCutoff(ctemps) <- c(w2 = 0.3, o1 = 0.2, o2 = 0.2)
ctemps

## --------------------------------------------------------------
templateCutoff(ctemps) <- c(w2 = 0.3, default = 0.2)
ctemps

## --------------------------------------------------------------
templateCutoff(cdetects) <- c(w2 = 0.3, default = 0.2)

## --------------------------------------------------------------
cdetects

## ----out.width = "5.8in", fig.width = 7, fig.keep = "none"-----
plot(cdetects)

## --------------------------------------------------------------
cdetects <- cdetects[c("w2", "o2")]
cdetects

## ----out.width = "5.8in", fig.width = 7, fig.keep = 'none'-----
plot(cdetects)

## --------------------------------------------------------------
wbt1 <- makeBinTemplate(btnw.fp, amp.cutoff = -40, name = "w1")

## ----fig.keep = 'none'-----------------------------------------
wbt2 <- makeBinTemplate(btnw.fp, amp.cutoff = -30, t.lim = c(1.5, 2.1), frq.lim = c(4.2, 5.6), buffer = 2, name = "w2")

## ----fig.keep = 'none'-----------------------------------------
obt1 <- makeBinTemplate(oven.fp, amp.cutoff = -20, t.lim = c(1, 4), frq.lim = c(1, 11), name = "o1")
obt2 <- makeBinTemplate(oven.fp, amp.cutoff = -17, t.lim = c(1, 4), frq.lim = c(1, 11), buffer = 2, name = "o2")

## --------------------------------------------------------------
btemps <- combineBinTemplates(wbt1, wbt2, obt1, obt2)
btemps

## --------------------------------------------------------------
bscores <- binMatch(survey.fp, btemps)

## --------------------------------------------------------------
bdetects <- findPeaks(bscores)

## ----out.width = "5.8in", fig.width = 7, fig.keep = 'none'-----
plot(bdetects)

## --------------------------------------------------------------
bdetects <- bdetects[-1]

## --------------------------------------------------------------
templateCutoff(bdetects) <- c(w2 = 7, default = 4)

## ----out.width = "5.8in", fig.width = 7, fig.keep = 'none'-----
plot(bdetects)

## --------------------------------------------------------------
getDetections(bdetects)


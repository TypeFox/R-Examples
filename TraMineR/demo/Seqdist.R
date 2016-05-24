## ============================
## Examples for the chapter
## 'Measuring similarities and
## distances between sequences'
## in TraMineR User's Guide
## ============================

##require(grDevices); require(graphics)
##oask <- devAskNewPage(dev.interactive(orNone = TRUE))

require(TraMineR)

## ----
## MPOS
## ----
data(famform)
famform.seq <- seqdef(famform)
famform.seq

## number of matching positions
seqmpos(famform.seq[1,],famform.seq[2,])
seqmpos(famform.seq[2,],famform.seq[4,])

## ---
## LCP
## ---
data(famform)
famform.seq <- seqdef(famform)
famform.seq

seqLLCP(famform.seq[1,],famform.seq[2,])
seqLLCP(famform.seq[3,],famform.seq[4,])
seqLLCP(famform.seq[3,],famform.seq[5,])

seqdist(famform.seq,method="LCP")

seqdist(famform.seq,method="LCP",norm=TRUE)

## Deriving proximity measures from normalized distances
1-seqdist(famform.seq,method="LCP",norm=TRUE)

## ---
## LCS
## ---
data(biofam)
biofam.seq <- seqdef(biofam,10:25)
biofam.lcs <- seqdist(biofam.seq,method="LCS")

## --
## OM
## --
costs <- seqsubm(biofam.seq,method="TRATE")
costs

biofam.om <- seqdist(biofam.seq, method="OM", indel=3, sm=costs)

object.size(biofam.om)/1024^2

round(biofam.om[1:10,1:10],1)

## ----------
## LCS <=> OM
## ----------
ccosts <- seqsubm(biofam.seq,method="CONSTANT",cval=2)
ccosts

biofam.om2 <- seqdist(biofam.seq,method="OM",indel=1,sm=ccosts)
biofam.om2[1:10,1:10]

all.equal(biofam.om2,biofam.lcs)

## ---
## DHD
## ---

dhdcosts <- seqsubm(biofam.seq, method="TRATE", time.varying=TRUE)

biofam.dhd <- seqdist(biofam.seq, method="DHD", sm=dhdcosts)


##devAskNewPage(oask)

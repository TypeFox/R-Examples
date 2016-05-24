### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### demo/lsa_landauer.r
### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
# generate the files of the famous Landauer example

ldir = tempfile()
dir.create(ldir)
write( c("human", "interface", "computer"), file=paste(ldir, "c1", sep="/"))
write( c("survey", "user", "computer", "system", "response", "time"), file=paste(ldir, "c2", sep="/"))
write( c("EPS", "user", "interface", "system"), file=paste(ldir, "c3", sep="/"))
write( c("system", "human", "system", "EPS"), file=paste(ldir, "c4", sep="/"))
write( c("user", "response", "time"), file=paste(ldir, "c5", sep="/"))
write( c("trees"), file=paste(ldir, "m1", sep="/"))
write( c("graph", "trees"), file=paste(ldir, "m2", sep="/"))
write( c("graph", "minors", "trees"), file=paste(ldir, "m3", sep="/"))
write( c("graph", "minors", "survey"), file=paste(ldir, "m4", sep="/"))

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
# generate doc term matrix from landauer files

dtm = textmatrix(ldir, minWordLength=1)
dtm

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
# make a space, reconstruct original

landauerOriginalSpace = lsa(dtm, dims=dimcalc_raw())
X = as.textmatrix(landauerOriginalSpace)

# X should be equal to dtm (beside rounding errors)
all( (round(X,2) == dtm) == TRUE)

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
# reduce dimensionality (Y shall be 
# the recalculated 'reduced' matrix)

landauerSpace = lsa(dtm, dims=2)
Y = as.textmatrix(landauerSpace)
round(Y,2)

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
# now read in again the landauer sample (but 
# with the vocabulary of the existing matrix)

pdocs = textmatrix(ldir, vocabulary=rownames(dtm))

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
# now calc a pseudo SVD on the basis of dtm's SVD

Y2 = fold_in(pdocs, landauerSpace)
round(Y2,2)

# Y and Y2 should be the same (as well as 
# dtm and pdocs should be equal)

all( (round(Y,2) == round(Y2,2)) == TRUE)

# calc pearson doc2doc correlation

rawCor = cor(dtm)
lsaCor = cor(Y)

# you should clearly see, that the "computer" documents (starting with "C")
# can in lsaCor be much better be differentiated from the "math" documents
# (starting with "m"). Moreover, the computer and math documents respectively
# have become more similar within their group.

round(rawCor,2)
round(lsaCor,2)

# clean up
unlink(ldir, recursive=TRUE)


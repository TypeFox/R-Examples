library(longmemo)

options(digits = 5) # since we want to "R CMD Rdiff" the output
data(NileMin)
fF <- FEXPest(NileMin, order.poly=3, pvalmax= .5, verbose=TRUE)
fF
plot(fF)

data(videoVBR)
for(max.poly in 0:3) {
    cat("-------------------------------------------\n",
        "max.poly= ", max.poly,":\n------------\n",sep='')
    fv <- FEXPest(videoVBR, order.poly=max.poly, pvalmax= .5, verbose=TRUE)
    print(fv["H"])
}

fv
plot(fv, type = "o", cex = 0.5)

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

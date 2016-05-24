### R code from vignette source 'protr.Rnw'

###################################################
### code chunk number 1: prelim
###################################################
protr.version = '1.1-0'
now.date = strftime(Sys.Date(), "%Y-%m-%d")


###################################################
### code chunk number 2: extractAAC
###################################################
require(protr)
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
extractAAC(x)


###################################################
### code chunk number 3: extractDC
###################################################
dc = extractDC(x)
head(dc, n = 30L)


###################################################
### code chunk number 4: extractTC
###################################################
tc = extractTC(x)
head(tc, n = 36L)


###################################################
### code chunk number 5: extractMoreau1
###################################################
moreau = extractMoreauBroto(x)
head(moreau, n = 36L)


###################################################
### code chunk number 6: extractMoreau2
###################################################
# Define 3 custom properties
myprops = data.frame(AccNo = c("MyProp1", "MyProp2", "MyProp3"),
                     A = c(0.62,  -0.5, 15),  R = c(-2.53,   3, 101),
                     N = c(-0.78,  0.2, 58),  D = c(-0.9,    3, 59),
                     C = c(0.29,    -1, 47),  E = c(-0.74,   3, 73),
                     Q = c(-0.85,  0.2, 72),  G = c(0.48,    0, 1),
                     H = c(-0.4,  -0.5, 82),  I = c(1.38, -1.8, 57),
                     L = c(1.06,  -1.8, 57),  K = c(-1.5,    3, 73),
                     M = c(0.64,  -1.3, 75),  F = c(1.19, -2.5, 91),
                     P = c(0.12,     0, 42),  S = c(-0.18, 0.3, 31),
                     T = c(-0.05, -0.4, 45),  W = c(0.81, -3.4, 130),
                     Y = c(0.26,  -2.3, 107), V = c(1.08, -1.5, 43))

# Use 4 properties in the AAindex database, and 3 cutomized properties
moreau2 = extractMoreauBroto(x, customprops = myprops,
                             props = c('CIDH920105', 'BHAR880101',
                                       'CHAM820101', 'CHAM820102',
                                       'MyProp1', 'MyProp2', 'MyProp3'))
head(moreau2, n = 36L)


###################################################
### code chunk number 7: extractMoran
###################################################
# Use the 3 custom properties defined before
# and 4 properties in the AAindex database
moran = extractMoran(x, customprops = myprops,
                     props = c('CIDH920105', 'BHAR880101',
                               'CHAM820101', 'CHAM820102',
                               'MyProp1', 'MyProp2', 'MyProp3'))
head(moran, n = 36L)


###################################################
### code chunk number 8: extractGeary
###################################################
# Use the 3 custom properties defined before
# and 4 properties in the AAindex database
geary = extractGeary(x, customprops = myprops,
                     props = c('CIDH920105', 'BHAR880101',
                               'CHAM820101', 'CHAM820102',
                               'MyProp1', 'MyProp2', 'MyProp3'))
head(geary, n = 36L)


###################################################
### code chunk number 9: extractCTDC
###################################################
extractCTDC(x)


###################################################
### code chunk number 10: extractCTDT
###################################################
extractCTDT(x)


###################################################
### code chunk number 11: extractCTDD
###################################################
extractCTDD(x)


###################################################
### code chunk number 12: extractCTriad
###################################################
ctriad = extractCTriad(x)
head(ctriad, n = 65L)


###################################################
### code chunk number 13: extractSOCN
###################################################
extractSOCN(x)


###################################################
### code chunk number 14: extractQSO
###################################################
extractQSO(x)


###################################################
### code chunk number 15: extractPAAC
###################################################
extractPAAC(x)


###################################################
### code chunk number 16: extractAPAAC
###################################################
extractAPAAC(x)


###################################################
### code chunk number 17: extractDescScales
###################################################
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
descscales = extractDescScales(x, propmat = 'AATopo',
                               index = c(37:41, 43:47),
                               pc = 5, lag = 7, silent = FALSE)


###################################################
### code chunk number 18: extractDescScales2
###################################################
length(descscales)
head(descscales, 15)


###################################################
### code chunk number 19: extractBLOSUM
###################################################
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
blosum = extractBLOSUM(x, submat = 'AABLOSUM62',
                       k = 5, lag = 7, scale = TRUE, silent = FALSE)


###################################################
### code chunk number 20: extractBLOSUM2
###################################################
length(blosum)
head(blosum, 15)


###################################################
### code chunk number 21: protcheck
###################################################
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
# A real sequence
protcheck(x)
# An artificial sequence
protcheck(paste(x, 'Z', sep = ''))


###################################################
### code chunk number 22: protseg
###################################################
protseg(x, aa = 'M', k = 5)



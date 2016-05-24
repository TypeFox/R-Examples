### R code from vignette source 'qualityTools.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: qualityTools.rnw:176-177
###################################################
library(qualityTools)


###################################################
### code chunk number 2: qualityTools.rnw:181-185
###################################################
#create artificial defect data set
defects = c(rep("E", 62), rep("B", 15), rep("F", 3), rep("A", 10),
rep("C",20), rep("D", 10))
paretoChart(defects)


###################################################
### code chunk number 3: qualityTools.rnw:262-266
###################################################
x = c(9.991, 10.013, 10.001, 10.007, 10.010, 10.013, 10.008, 10.017, 10.005,
10.005, 10.002, 10.017, 10.005, 10.002,  9.996, 10.011, 10.009, 10.006,
10.008, 10.003, 10.002, 10.006, 10.010, 9.992, 10.013)
cg(x, target = 10.003, tolerance = c(9.903, 10.103))


###################################################
### code chunk number 4: qualityTools.rnw:288-296
###################################################
#create a gage RnR design
design = gageRRDesign(Operators=3, Parts=10, Measurements=2, randomize=FALSE)
#set the response	
response(design) = c(23,22,22,22,22,25,23,22,23,22,20,22,22,22,24,25,27,28,
23,24,23,24,24,22,22,22,24,23,22,24,20,20,25,24,22,24,21,20,21,22,21,22,21,
21,24,27,25,27,23,22,25,23,23,22,22,23,25,21,24,23)
#perform Gage RnR
gdo = gageRR(design)


###################################################
### code chunk number 5: GageRR
###################################################
#visualization of Gage RnR
plot(gdo)


###################################################
### code chunk number 6: qualityTools.rnw:306-307
###################################################
#visualization of Gage RnR
plot(gdo)


###################################################
### code chunk number 7: NormPcr
###################################################
set.seed(1234)
#generate some data
norm = rnorm(20, mean = 20)
#generate some data
weib = rweibull(20, shape = 2, scale = 8)
#process capability
pcr(norm, "normal", lsl = 17, usl = 23)


###################################################
### code chunk number 8: WeibPcr
###################################################
#process cabapility
pcr(weib, "weibull", usl = 20)


###################################################
### code chunk number 9: qualityTools.rnw:406-407
###################################################
set.seed(1234)
#generate some data
norm = rnorm(20, mean = 20)
#generate some data
weib = rweibull(20, shape = 2, scale = 8)
#process capability
pcr(norm, "normal", lsl = 17, usl = 23)


###################################################
### code chunk number 10: qualityTools.rnw:412-413
###################################################
#process cabapility
pcr(weib, "weibull", usl = 20)


###################################################
### code chunk number 11: qualityTools.rnw:425-426
###################################################
pcr(weib, "weibull", usl = 20)


###################################################
### code chunk number 12: qqPlot
###################################################
par(mfrow = c(1,2))
qqPlot(weib, "weibull"); qqPlot(weib, "normal")


###################################################
### code chunk number 13: qualityTools.rnw:437-438
###################################################
par(mfrow = c(1,2))
qqPlot(weib, "weibull"); qqPlot(weib, "normal")


###################################################
### code chunk number 14: ppPlot
###################################################
par(mfrow = c(1,2))
ppPlot(norm, "weibull"); ppPlot(norm, "normal")


###################################################
### code chunk number 15: qualityTools.rnw:452-453
###################################################
par(mfrow = c(1,2))
ppPlot(norm, "weibull"); ppPlot(norm, "normal")


###################################################
### code chunk number 16: qualityTools.rnw:491-497
###################################################
set.seed(1234)
fdo = facDesign(k = 3, centerCube = 4) #fdo - factorial design object
names(fdo) = c("Factor 1", "Factor 2", "Factor 3") #optional
lows(fdo) = c(80, 120, 1) #optional
highs(fdo) = c(120, 140, 2) #optional
summary(fdo) #information about the factorial design


###################################################
### code chunk number 17: qualityTools.rnw:503-505
###################################################
#set first value
yield = simProc(x1 = 120, x2 = 140, x3 = 2)


###################################################
### code chunk number 18: qualityTools.rnw:517-521
###################################################
yield = c(simProc(120,140, 1),simProc(80,140, 1),simProc(120,140, 2),
simProc(120,120, 1),simProc(90,130, 1.5),simProc(90,130, 1.5),
simProc(80,120, 2),simProc(90,130, 1.5),simProc(90,130, 1.5),
simProc(120,120, 2),simProc(80,140, 2),simProc(80,120, 1))


###################################################
### code chunk number 19: qualityTools.rnw:527-528
###################################################
response(fdo) = yield	#assign yield to the factorial design object


###################################################
### code chunk number 20: effPlot
###################################################
effectPlot(fdo, classic = TRUE)


###################################################
### code chunk number 21: iaPlot
###################################################
interactionPlot(fdo)


###################################################
### code chunk number 22: qualityTools.rnw:545-546
###################################################
effectPlot(fdo, classic = TRUE)


###################################################
### code chunk number 23: qualityTools.rnw:551-552
###################################################
interactionPlot(fdo)


###################################################
### code chunk number 24: qualityTools.rnw:564-566
###################################################
lm.1 = lm(yield ~ A*B*C, data = fdo)
summary(lm.1)


###################################################
### code chunk number 25: ParetPlot
###################################################
par(mfrow = c(1,2))
paretoPlot(fdo)
normalPlot(fdo)


###################################################
### code chunk number 26: qualityTools.rnw:579-580
###################################################
par(mfrow = c(1,2))
paretoPlot(fdo)
normalPlot(fdo)


###################################################
### code chunk number 27: wirePlot
###################################################
par(mfrow = c(1,2))
wirePlot(A, B, yield, data = fdo)
contourPlot(A, B, yield, data = fdo)


###################################################
### code chunk number 28: qualityTools.rnw:595-596
###################################################
par(mfrow = c(1,2))
wirePlot(A, B, yield, data = fdo)
contourPlot(A, B, yield, data = fdo)


###################################################
### code chunk number 29: qualityTools.rnw:612-613
###################################################
fdo.frac = fracDesign(k = 3, gen = "C = AB", centerCube = 4)


###################################################
### code chunk number 30: qualityTools.rnw:618-619
###################################################
summary(fdo.frac)


###################################################
### code chunk number 31: qualityTools.rnw:632-633
###################################################
aliasTable(fdo.frac)


###################################################
### code chunk number 32: qualityTools.rnw:638-639
###################################################
confounds(fdo.frac)


###################################################
### code chunk number 33: fracChoose
###################################################
fracChoose()


###################################################
### code chunk number 34: qualityTools.rnw:649-650
###################################################
fracChoose()


###################################################
### code chunk number 35: qualityTools.rnw:660-661
###################################################
fdo1 = facDesign(k = 3, centerCube = 2, replicates = 2)


###################################################
### code chunk number 36: qualityTools.rnw:668-671
###################################################
set.seed(1234)
y2 = rnorm(12, mean = 20)
response(fdo) = data.frame(yield, y2) 


###################################################
### code chunk number 37: multiresp
###################################################
par(mfrow = c(1,2))
wirePlot(A, B, yield, data = fdo, form = "yield~A+B+C+A*B")
contourPlot(A, B, y2, data = fdo, form = "y2~A+B+C+A*B")


###################################################
### code chunk number 38: qualityTools.rnw:683-684
###################################################
par(mfrow = c(1,2))
wirePlot(A, B, yield, data = fdo, form = "yield~A+B+C+A*B")
contourPlot(A, B, y2, data = fdo, form = "y2~A+B+C+A*B")


###################################################
### code chunk number 39: multiresp2
###################################################
par(mfrow = c(1,2))
wirePlot(A,B,y2, data = fdo, factors = list(C=-1), form = "y2~A*B*C")
wirePlot(A,B,y2, data = fdo, factors = list(C=1), form = "y2~A*B*C")


###################################################
### code chunk number 40: qualityTools.rnw:699-700
###################################################
par(mfrow = c(1,2))
wirePlot(A,B,y2, data = fdo, factors = list(C=-1), form = "y2~A*B*C")
wirePlot(A,B,y2, data = fdo, factors = list(C=1), form = "y2~A*B*C")


###################################################
### code chunk number 41: qualityTools.rnw:708-711
###################################################
fits(fdo) = lm(yield ~ A+B, data = fdo)
fits(fdo) = lm(y2 ~ A*B*C, data = fdo)
fits(fdo)


###################################################
### code chunk number 42: qualityTools.rnw:718-720
###################################################
sao =steepAscent(factors=c("A","B"),response="yield",data=fdo,steps=20)
sao


###################################################
### code chunk number 43: steepAsc
###################################################
predicted = simProc(sao[,5], sao[,6])
response(sao) = predicted 
plot(sao, type = "b", col = 2)


###################################################
### code chunk number 44: qualityTools.rnw:732-733
###################################################
predicted = simProc(sao[,5], sao[,6])
response(sao) = predicted 
plot(sao, type = "b", col = 2)


###################################################
### code chunk number 45: qualityTools.rnw:747-753
###################################################
#set the seed for randomization of the runs
set.seed(1234)
fdo2 = facDesign(k = 2, centerCube = 3)
names(fdo2) = c("Factor 1", "Factor 2") 
lows(fdo2) = c(134, 155)
highs(fdo2) = c(155, 175)


###################################################
### code chunk number 46: qualityTools.rnw:758-762
###################################################
yield = c(simProc(134,175), simProc(144.5,165.5), simProc(155,155),
simProc(144.5,165.5), simProc(155,175), simProc(144.5,165.5),
simProc(134,155))
response(fdo2) = yield


###################################################
### code chunk number 47: qualityTools.rnw:769-771
###################################################
rsdo = starDesign(data = fdo2)
rsdo


###################################################
### code chunk number 48: qualityTools.rnw:776-779
###################################################
yield2 = c(yield, simProc(130,165), simProc(155,165), simProc(144,155),
simProc(144,179),simProc(144,165),simProc(144,165),simProc(144,165))
response(rsdo) = yield2


###################################################
### code chunk number 49: qualityTools.rnw:784-785
###################################################
lm.3 = lm(yield2 ~ A*B + I(A^2) + I(B^2), data = rsdo)


###################################################
### code chunk number 50: respSur
###################################################
par(mfrow=c(1,2))
wirePlot(A,B,yield2,form="yield2~A*B+I(A^2)+I(B^2)",data=rsdo,theta=-70)
contourPlot(A,B,yield2,form="yield2~A*B+I(A^2)+I(B^2)",data=rsdo)


###################################################
### code chunk number 51: qualityTools.rnw:797-798
###################################################
par(mfrow=c(1,2))
wirePlot(A,B,yield2,form="yield2~A*B+I(A^2)+I(B^2)",data=rsdo,theta=-70)
contourPlot(A,B,yield2,form="yield2~A*B+I(A^2)+I(B^2)",data=rsdo)


###################################################
### code chunk number 52: qualityTools.rnw:807-811
###################################################
A = seq(40, 210, length = 100)
B = seq(90, 190, length = 100)
C = seq(90, 190, length = 100)
filled.contour(A, B,outer(A,B, simProc, noise = FALSE), xlab = "Factor 1", ylab = "Factor 2", color.palette = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))) 


###################################################
### code chunk number 53: qualityTools.rnw:821-822
###################################################
fdo = rsmDesign(k = 3, alpha = 1.633, cc = 0, cs = 6)


###################################################
### code chunk number 54: qualityTools.rnw:827-828
###################################################
fdo = randomize(fdo, so = TRUE)


###################################################
### code chunk number 55: desTab
###################################################
rsdo = rsmChoose()


###################################################
### code chunk number 56: qualityTools.rnw:838-839
###################################################
rsdo = rsmChoose()


###################################################
### code chunk number 57: qualityTools.rnw:849-851
###################################################
fdo3 = facDesign(k = 6)
rsdo = starDesign(alpha = "orthogonal", data = fdo3)


###################################################
### code chunk number 58: qualityTools.rnw:860-861
###################################################
randomize(fdo, random.seed = 123)


###################################################
### code chunk number 59: qualityTools.rnw:865-866
###################################################
randomize(fdo, so = TRUE)


###################################################
### code chunk number 60: qualityTools.rnw:903-906
###################################################
d1 = desirability(y1, 120, 170, scale = c(1,1), target = "max")
d3 = desirability(y3, 400, 600, target = 500)
d1


###################################################
### code chunk number 61: desirb
###################################################
par(mfrow = c(1,2))
plot(d1, col = 2); plot(d3, col = 2)


###################################################
### code chunk number 62: qualityTools.rnw:917-918
###################################################
par(mfrow = c(1,2))
plot(d1, col = 2); plot(d3, col = 2)


###################################################
### code chunk number 63: qualityTools.rnw:931-939
###################################################
ddo = rsmDesign(k = 3, alpha = 1.633, cc = 0, cs = 6)
ddo = randomize(ddo, so = TRUE)
#optional
names(ddo) = c("silica", "silan", "sulfur")
#optional
highs(ddo) = c(1.7, 60, 2.8)
#optional
lows(ddo) = c(0.7, 40, 1.8)


###################################################
### code chunk number 64: qualityTools.rnw:944-952
###################################################
y1 = c(102, 120, 117, 198, 103, 132, 132, 139, 102, 154, 96, 163, 116,
 153, 133, 133, 140, 142, 145, 142)
y2 = c(900, 860, 800, 2294, 490, 1289, 1270, 1090, 770, 1690, 700, 1540,
 2184, 1784, 1300, 1300, 1145, 1090, 1260, 1344)
y3 = c(470, 410, 570, 240, 640, 270, 410, 380, 590, 260, 520, 380, 520,
 290, 380, 380, 430, 430, 390, 390)
y4 = c(67.5, 65, 77.5, 74.5, 62.5, 67, 78, 70, 76, 70, 63, 75, 65, 71,
 70, 68.5, 68, 68, 69, 70)


###################################################
### code chunk number 65: qualityTools.rnw:957-958
###################################################
response(ddo) = data.frame(y1, y2, y3, y4)[c(5,2,3,8,1,6,7,4,9:20),]


###################################################
### code chunk number 66: qualityTools.rnw:963-965
###################################################
d2 = desirability(y2, 1000, 1300, target = "max")
d4 = desirability(y4, 60, 75, target = 67.5)


###################################################
### code chunk number 67: qualityTools.rnw:970-971
###################################################
desires(ddo)=d1; desires(ddo)=d2; desires(ddo)=d3; desires(ddo)=d4


###################################################
### code chunk number 68: qualityTools.rnw:976-980
###################################################
fits(ddo) = lm(y1 ~ A+B+C+A:B+A:C+B:C+I(A^2)+I(B^2)+I(C^2), data = ddo)
fits(ddo) = lm(y2 ~ A+B+C+A:B+A:C+B:C+I(A^2)+I(B^2)+I(C^2), data = ddo)
fits(ddo) = lm(y3 ~ A+B+C+A:B+A:C+B:C+I(A^2)+I(B^2)+I(C^2), data = ddo)
fits(ddo) = lm(y4 ~ A+B+C+A:B+A:C+B:C+I(A^2)+I(B^2)+I(C^2), data = ddo)


###################################################
### code chunk number 69: qualityTools.rnw:985-986
###################################################
optimum(ddo, type = "optim")


###################################################
### code chunk number 70: qualityTools.rnw:995-1003
###################################################
mdo = mixDesign(3,2, center = FALSE, axial = FALSE, randomize = FALSE,
replicates  = c(1,1,2,3))
names(mdo) = c("polyethylene", "polystyrene", "polypropylene")

#set response (i.e. yarn elongation)
elongation = c(11.0, 12.4, 15.0, 14.8, 16.1, 17.7, 16.4, 16.6, 8.8, 10.0,
10.0, 9.7, 11.8, 16.8, 16.0)  
response(mdo) = elongation


###################################################
### code chunk number 71: qualityTools.rnw:1008-1009
###################################################
mdo


###################################################
### code chunk number 72: respAcont
###################################################
par(mfrow=c(1,2))
contourPlot3(A, B, C, elongation, data = mdo, form = "quadratic")
wirePlot3(A, B, C, elongation, data=mdo, form="quadratic", theta=-170)


###################################################
### code chunk number 73: qualityTools.rnw:1021-1022
###################################################
par(mfrow=c(1,2))
contourPlot3(A, B, C, elongation, data = mdo, form = "quadratic")
wirePlot3(A, B, C, elongation, data=mdo, form="quadratic", theta=-170)


###################################################
### code chunk number 74: qualityTools.rnw:1040-1046
###################################################
set.seed(1234)
tdo = taguchiDesign("L9_3")
values(tdo) = list(A  = c(20, 40, 60), B = c("mateial 1", "material 2",
"material 3"), C = c(1,2,3))
names(tdo) = c("Factor 1", "Factor 2", "Factor 3", "Factor 4")
summary(tdo)


###################################################
### code chunk number 75: taguchi
###################################################
response(tdo) = rnorm(9)
effectPlot(tdo, col = 2)


###################################################
### code chunk number 76: qualityTools.rnw:1057-1058
###################################################
response(tdo) = rnorm(9)
effectPlot(tdo, col = 2)


###################################################
### code chunk number 77: qualityTools.rnw:1095-1096
###################################################
sessionInfo()



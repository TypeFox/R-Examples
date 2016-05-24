### R code from vignette source 'PRROC.Rnw'

###################################################
### code chunk number 1: PRROC.Rnw:33-35
###################################################
library(PRROC)
set.seed(127)


###################################################
### code chunk number 2: PRROC.Rnw:46-48
###################################################
fg<-rnorm(300);
bg<-rnorm(500,-2);


###################################################
### code chunk number 3: PRROC.Rnw:54-56
###################################################
roc<-roc.curve(scores.class0 = fg, scores.class1 = bg)
pr<-pr.curve(scores.class0 = fg, scores.class1 = bg)


###################################################
### code chunk number 4: PRROC.Rnw:60-61
###################################################
roc


###################################################
### code chunk number 5: PRROC.Rnw:66-67
###################################################
pr


###################################################
### code chunk number 6: PRROC.Rnw:72-74
###################################################
roc<-roc.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
pr<-pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)


###################################################
### code chunk number 7: PRROC.Rnw:78-80
###################################################
roc
pr


###################################################
### code chunk number 8: roc
###################################################
plot(roc)


###################################################
### code chunk number 9: PRROC.Rnw:91-92
###################################################
plot(roc)


###################################################
### code chunk number 10: pr
###################################################
plot(pr)


###################################################
### code chunk number 11: PRROC.Rnw:103-104
###################################################
plot(pr)


###################################################
### code chunk number 12: PRROC.Rnw:114-116
###################################################
x<-c(fg,bg);
lab<-c(rep(1,length(fg)),rep(0,length(bg)))


###################################################
### code chunk number 13: PRROC.Rnw:121-123
###################################################
roc<-roc.curve(scores.class0 = x, weights.class0 = lab);
pr<-pr.curve(scores.class0 = x, weights.class0 = lab);


###################################################
### code chunk number 14: PRROC.Rnw:128-130
###################################################
roc
pr


###################################################
### code chunk number 15: PRROC.Rnw:145-146
###################################################
wfg<- c(runif(300,min=0.5,max=1),runif(500,min=0,max=0.5))


###################################################
### code chunk number 16: PRROC.Rnw:151-153
###################################################
hist(wfg[301:800],col=2,xlim=c(0,1),main="Weights",xlab="foreground weight");
hist(wfg[1:300],col=3,add=T);


###################################################
### code chunk number 17: PRROC.Rnw:160-162
###################################################
wroc<-roc.curve(scores.class0 = x, weights.class0 = wfg, curve = TRUE)
wpr<-pr.curve(scores.class0 = x, weights.class0 = wfg, curve = TRUE)


###################################################
### code chunk number 18: PRROC.Rnw:168-172
###################################################
wroc<-roc.curve(scores.class0 = x, scores.class1 = x, 
  weights.class0 = wfg, weights.class1 = 1-wfg, curve = TRUE)
wpr<-pr.curve(scores.class0 = x, scores.class1 = x, 
  weights.class0 = wfg,weights.class1 = 1-wfg, curve = TRUE)


###################################################
### code chunk number 19: wroc
###################################################
plot(wroc)


###################################################
### code chunk number 20: PRROC.Rnw:181-182
###################################################
plot(wroc)


###################################################
### code chunk number 21: wpr
###################################################
plot(wpr)


###################################################
### code chunk number 22: PRROC.Rnw:192-193
###################################################
plot(wpr)


###################################################
### code chunk number 23: PRROC.Rnw:202-206
###################################################
wpr<-pr.curve(scores.class0 = x, weights.class0 = wfg, curve = TRUE, 
  max.compute = T, min.compute = T, rand.compute = T)
wroc<-roc.curve(scores.class0 = x, weights.class0 = wfg, curve = TRUE, 
  max.compute = T, min.compute = T, rand.compute = T)


###################################################
### code chunk number 24: PRROC.Rnw:211-213
###################################################
wpr
wroc


###################################################
### code chunk number 25: wpr2
###################################################
plot(wpr,max.plot = TRUE, min.plot = TRUE, rand.plot = TRUE, 
  fill.area = TRUE)


###################################################
### code chunk number 26: PRROC.Rnw:224-225
###################################################
plot(wpr,max.plot = TRUE, min.plot = TRUE, rand.plot = TRUE, 
  fill.area = TRUE)


###################################################
### code chunk number 27: PRROC.Rnw:232-238
###################################################
y<-c(rnorm(300,sd=2),rnorm(500,-5,sd=2))

wpr2<-pr.curve(scores.class0 = y, weights.class0 = wfg, curve = TRUE, 
  max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
wroc2<-roc.curve(scores.class0 = y, weights.class0 = wfg, curve = TRUE, 
  max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)


###################################################
### code chunk number 28: PRROC.Rnw:243-245
###################################################
plot(wpr, max.plot = TRUE, min.plot = TRUE, rand.plot = TRUE, 
  fill.area = T, color=2, auc.main = FALSE);


###################################################
### code chunk number 29: PRROC.Rnw:249-250
###################################################
plot(wpr2, add = TRUE, color = 3);


###################################################
### code chunk number 30: PRROC.Rnw:255-257
###################################################
plot(wpr, max.plot = TRUE, min.plot = TRUE, rand.plot = TRUE, fill.area = T, color=2, auc.main = FALSE);
plot(wpr2, add = TRUE, color = 3);


###################################################
### code chunk number 31: plot1
###################################################
plot(wpr,scale.color = heat.colors(100));


###################################################
### code chunk number 32: PRROC.Rnw:273-274
###################################################
plot(wpr,scale.color = heat.colors(100));


###################################################
### code chunk number 33: plot2
###################################################
plot(wpr, auc.main = FALSE, main = "My classifier")


###################################################
### code chunk number 34: PRROC.Rnw:284-285
###################################################
plot(wpr, auc.main = FALSE, main = "My classifier")


###################################################
### code chunk number 35: plot6
###################################################
plot(wpr, legend = FALSE)


###################################################
### code chunk number 36: PRROC.Rnw:295-296
###################################################
plot(wpr, legend = FALSE)


###################################################
### code chunk number 37: plot3
###################################################
plot(wpr, color=3, lty="dotted");


###################################################
### code chunk number 38: PRROC.Rnw:306-307
###################################################
plot(wpr, color=3, lty="dotted");


###################################################
### code chunk number 39: plot4
###################################################
plot(wpr,legend=1);


###################################################
### code chunk number 40: PRROC.Rnw:317-318
###################################################
plot(wpr,legend=1);


###################################################
### code chunk number 41: plot5
###################################################
plot(wpr, rand.plot = TRUE, fill.area = TRUE, 
  fill.color = rgb(0.8,1,0.8), maxminrand.col = "blue" );


###################################################
### code chunk number 42: PRROC.Rnw:329-330
###################################################
plot(wpr, rand.plot = TRUE, fill.area = TRUE, 
  fill.color = rgb(0.8,1,0.8), maxminrand.col = "blue" );


###################################################
### code chunk number 43: PRROC.Rnw:338-339
###################################################
curve.points<-wpr$curve


###################################################
### code chunk number 44: PRROC.Rnw:343-344
###################################################
curve.points[1:5,]


###################################################
### code chunk number 45: plot6
###################################################
plot(curve.points[,1],curve.points[,2],
	 xlab="Recall",ylab="Precision",t="l")


###################################################
### code chunk number 46: PRROC.Rnw:357-358
###################################################
plot(curve.points[,1],curve.points[,2],
	 xlab="Recall",ylab="Precision",t="l")


###################################################
### code chunk number 47: PRROC.Rnw:364-365
###################################################
library(ggplot2)


###################################################
### code chunk number 48: plot7
###################################################
( 
ggplot(data.frame(wpr$curve),aes(x=X1,y=X2,color=X3)) 
	+ geom_line() 
	+ labs(x="Recall",y="Precision",
		   title=format(wpr$auc.integral,digits=3),
		   colour="Threshold") 
	+ scale_colour_gradient2(low="red", mid="orange",high="yellow")
)


###################################################
### code chunk number 49: PRROC.Rnw:380-381
###################################################
( 
ggplot(data.frame(wpr$curve),aes(x=X1,y=X2,color=X3)) 
	+ geom_line() 
	+ labs(x="Recall",y="Precision",
		   title=format(wpr$auc.integral,digits=3),
		   colour="Threshold") 
	+ scale_colour_gradient2(low="red", mid="orange",high="yellow")
)



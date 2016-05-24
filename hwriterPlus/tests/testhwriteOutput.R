require(hwriterPlus)
### Test capture of output only first
### Annette Dobson (1990) "An Introduction to Generalized Linear Models".
### Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2,10,20, labels=c("Ctl","Trt"))
weight <- c(ctl, trt)
plants <- data.frame(weight = weight, group = group)
strOutput <- capture.output(str(plants))
lm.D9 <- lm(weight ~ group, data = plants)
anovaOutput <- capture.output(anova(lm.D9))
hwrite("Analysis of Weight Data", heading = 1,
       center = FALSE, br = TRUE)
hwriteOutput(strOutput, fontSize = "8pt")
hwriteOutput(anovaOutput)

### Test capture of a whole session
tmpFile <- tempfile("Session.txt")
script(tmpFile)
clotting <-
    data.frame(
               u = c(5,10,15,20,30,40,60,80,100),
               lot1 = c(118,58,42,35,27,25,21,19,18),
               lot2 = c(69,35,26,21,18,16,13,12,12)
               )
clotting
coef(glm(lot1 ~ log(u), data=clotting, family=Gamma))
q()
sessionOut <- readLines(tmpFile)
sessionOut
hwriteOutput(sessionOut, center = FALSE, br = TRUE)

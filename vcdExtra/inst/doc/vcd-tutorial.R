### R code from vignette source 'vcd-tutorial.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
set.seed(1071)
#library(vcd)
library(vcdExtra)
library(ggplot2)
#data(Titanic)
data(HairEyeColor)
data(PreSex)
data(Arthritis)
art <- xtabs(~Treatment + Improved, data = Arthritis)
if(!file.exists("fig")) dir.create("fig")


###################################################
### code chunk number 2: case-form
###################################################
names(Arthritis)      # show the variables
str(Arthritis)        # show the structure
head(Arthritis,5)     # first 5 observations, same as Arthritis[1:5,] 


###################################################
### code chunk number 3: frequency-form
###################################################
# Agresti (2002), table 3.11, p. 106
GSS <- data.frame(
  expand.grid(sex=c("female", "male"), 
              party=c("dem", "indep", "rep")),
  count=c(279,165,73,47,225,191))
GSS
names(GSS)
str(GSS)
sum(GSS$count)


###################################################
### code chunk number 4: table-form1
###################################################
str(HairEyeColor)                      # show the structure
sum(HairEyeColor)                      # number of cases
sapply(dimnames(HairEyeColor), length) # table dimension sizes


###################################################
### code chunk number 5: table-form2
###################################################
## A 4 x 4 table  Agresti (2002, Table 2.8, p. 57) Job Satisfaction
JobSat <- matrix(c(1,2,1,0, 3,3,6,1, 10,10,14,9, 6,7,12,11), 4, 4)
dimnames(JobSat) = list(income=c("< 15k", "15-25k", "25-40k", "> 40k"),
                satisfaction=c("VeryD", "LittleD", "ModerateS", "VeryS"))
JobSat


###################################################
### code chunk number 6: table-form3
###################################################
JobSat <- as.table(JobSat)
str(JobSat)


###################################################
### code chunk number 7: relevel (eval = FALSE)
###################################################
## dimnames(JobSat)$income<-c(7.5,20,32.5,60)
## dimnames(JobSat)$satisfaction<-1:4


###################################################
### code chunk number 8: reorder1
###################################################
HairEyeColor <- HairEyeColor[, c(1,3,4,2), ]
str(HairEyeColor)


###################################################
### code chunk number 9: reorder2 (eval = FALSE)
###################################################
## Arthritis <- read.csv("arthritis.txt",header=TRUE)
## Arthritis$Improved <- ordered(Arthritis$Improved, levels=c("None", "Some", "Marked"))


###################################################
### code chunk number 10: Arthritis
###################################################
mosaic(art, gp = shading_max, split_vertical = TRUE, main="Arthritis: [Treatment] [Improved]")


###################################################
### code chunk number 11: reorder3
###################################################
UCB <- aperm(UCBAdmissions, c(2, 1, 3))
dimnames(UCB)[[2]] <- c("Yes", "No")
names(dimnames(UCB)) <- c("Sex", "Admit?", "Department")
ftable(UCB)


###################################################
### code chunk number 12: structable
###################################################
structable(HairEyeColor)                   # show the table: default
structable(Hair+Sex ~ Eye, HairEyeColor)   # specify col ~ row variables


###################################################
### code chunk number 13: structable1 (eval = FALSE)
###################################################
## HSE < - structable(Hair+Sex ~ Eye, HairEyeColor)   # save structable object
## mosaic(HSE)                                        # plot it


###################################################
### code chunk number 14: setup
###################################################
 n=500
 A <- factor(sample(c("a1","a2"), n, rep=TRUE))
 B <- factor(sample(c("b1","b2"), n, rep=TRUE))
 C <- factor(sample(c("c1","c2"), n, rep=TRUE))
 mydata <- data.frame(A,B,C)


###################################################
### code chunk number 15: table-ex1
###################################################
# 2-Way Frequency Table
attach(mydata)
mytable <- table(A,B)   # A will be rows, B will be columns
mytable                 # print table

margin.table(mytable, 1) # A frequencies (summed over B)
margin.table(mytable, 2) # B frequencies (summed over A)

prop.table(mytable)    # cell percentages
prop.table(mytable, 1) # row percentages
prop.table(mytable, 2) # column percentages


###################################################
### code chunk number 16: table-ex2
###################################################
# 3-Way Frequency Table
mytable <- table(A, B, C)
ftable(mytable)


###################################################
### code chunk number 17: xtabs-ex1
###################################################
# 3-Way Frequency Table
mytable <- xtabs(~A+B+C, data=mydata)
ftable(mytable)    # print table
summary(mytable)   # chi-square test of indepedence


###################################################
### code chunk number 18: xtabs-ex2
###################################################
(GSStab <- xtabs(count ~ sex + party, data=GSS))
summary(GSStab)


###################################################
### code chunk number 19: dayton1
###################################################
str(DaytonSurvey)
head(DaytonSurvey)


###################################################
### code chunk number 20: dayton2
###################################################
# data in frequency form
# collapse over sex and race
Dayton.ACM.df <- aggregate(Freq ~ cigarette+alcohol+marijuana, 
                           data=DaytonSurvey, FUN=sum)
Dayton.ACM.df


###################################################
### code chunk number 21: dayton3
###################################################
# in table form
Dayton.tab <- xtabs(Freq~cigarette+alcohol+marijuana+sex+race, data=DaytonSurvey)
structable(cigarette+alcohol+marijuana ~ sex+race, data=Dayton.tab)


###################################################
### code chunk number 22: dayton4
###################################################
# collapse over sex and race
Dayton.ACM.tab <- apply(Dayton.tab, MARGIN=1:3, FUN=sum)
Dayton.ACM.tab <- margin.table(Dayton.tab, 1:3)   # same result
structable(cigarette+alcohol ~ marijuana, data=Dayton.ACM.tab)


###################################################
### code chunk number 23: dayton5 (eval = FALSE)
###################################################
## Dayton.ACM.df <- ddply(DaytonSurvey, .(cigarette, alcohol, marijuana), 
##                        plyr::summarise, Freq=sum(Freq))


###################################################
### code chunk number 24: collapse1
###################################################
# create some sample data in frequency form
sex <- c("Male", "Female")
age <- c("10-19", "20-29",  "30-39", "40-49", "50-59", "60-69")
education <- c("low", 'med', 'high')
data <- expand.grid(sex=sex, age=age, education=education)
counts <- rpois(36, 100)   # random Possion cell frequencies
data <- cbind(data, counts)

# make it into a 3-way table
t1 <- xtabs(counts ~ sex + age + education, data=data)
structable(t1)


###################################################
### code chunk number 25: collapse2
###################################################
# collapse age to 3 levels, education to 2 levels
t2 <- collapse.table(t1, 
         age=c("10-29", "10-29",  "30-49", "30-49", "50-69", "50-69"),
         education=c("<high", "<high", "high"))
structable(t2)


###################################################
### code chunk number 26: convert-ex1
###################################################
as.data.frame(GSStab)


###################################################
### code chunk number 27: convert-ex2
###################################################
Art.tab <-with(Arthritis, table(Treatment, Sex, Improved))
str(Art.tab)
ftable(Art.tab)


###################################################
### code chunk number 28: convert-ex3
###################################################
Art.df <- expand.dft(Art.tab)
str(Art.df)


###################################################
### code chunk number 29: tv1
###################################################
tv.data<-read.table(system.file("doc","extdata","tv.dat",package="vcdExtra"))
head(tv.data,5)


###################################################
### code chunk number 30: tv2 (eval = FALSE)
###################################################
## tv.data<-read.table("C:/R/data/tv.dat")


###################################################
### code chunk number 31: tv2
###################################################
TV <- array(tv.data[,5], dim=c(5,11,5,3))                                        
dimnames(TV) <- list(c("Monday","Tuesday","Wednesday","Thursday","Friday"), 
                c("8:00","8:15","8:30","8:45","9:00","9:15","9:30",         
                  "9:45","10:00","10:15","10:30"),                            
                c("ABC","CBS","NBC","Fox","Other"), c("Off","Switch","Persist"))
names(dimnames(TV))<-c("Day", "Time", "Network", "State")


###################################################
### code chunk number 32: tv2a (eval = FALSE)
###################################################
## TV <- xtabs(V5 ~ ., data=tv.data)
## dimnames(TV) <- list(Day=c("Monday","Tuesday","Wednesday","Thursday","Friday"), 
##                 Time=c("8:00","8:15","8:30","8:45","9:00","9:15","9:30",         
##                        "9:45","10:00","10:15","10:30"),                            
##                 Network=c("ABC","CBS","NBC","Fox","Other"), 
##                 State=c("Off","Switch","Persist"))


###################################################
### code chunk number 33: tv3
###################################################
TV <- TV[,,1:3,]     # keep only ABC, CBS, NBC
TV <- TV[,,,3]       # keep only Persist -- now a 3 way table
structable(TV)


###################################################
### code chunk number 34: tv4
###################################################
TV.df <- as.data.frame.table(TV)
levels(TV.df$Time) <- c(rep("8:00-8:59",4),rep("9:00-9:59",4), rep("10:00-10:44",3))
TV2 <- xtabs(Freq ~ Day + Time + Network, TV.df)
structable(Day ~ Time+Network,TV2)


###################################################
### code chunk number 35: xtabs-ex2
###################################################
# 2-Way Cross Tabulation
library(gmodels)
CrossTable(GSStab,prop.t=FALSE,prop.r=FALSE,prop.c=FALSE)


###################################################
### code chunk number 36: chisq
###################################################
(HairEye <- margin.table(HairEyeColor, c(1, 2)))
chisq.test(HairEye)


###################################################
### code chunk number 37: fisher
###################################################
fisher.test(GSStab)


###################################################
### code chunk number 38: vcd-tutorial.Rnw:840-841 (eval = FALSE)
###################################################
## fisher.test(HairEye)


###################################################
### code chunk number 39: mantel1
###################################################
## UC Berkeley Student Admissions
mantelhaen.test(UCBAdmissions)


###################################################
### code chunk number 40: mantel2
###################################################
oddsratio(UCBAdmissions, log=FALSE)
lor <- oddsratio(UCBAdmissions)  # capture log odds ratios
summary(lor)
woolf_test(UCBAdmissions) 


###################################################
### code chunk number 41: fourfold1
###################################################
col <- c("#99CCFF", "#6699CC", "#F9AFAF", "#6666A0", "#FF0000", "#000080")
fourfold(UCB,mfrow=c(2,3), color=col)


###################################################
### code chunk number 42: fourfold2 (eval = FALSE)
###################################################
## cotabplot(UCB, panel = cotab_fourfold)


###################################################
### code chunk number 43: doubledecker
###################################################
doubledecker(Admit ~ Dept + Gender, data=UCBAdmissions[2:1,,])


###################################################
### code chunk number 44: oddsratio0 (eval = FALSE)
###################################################
## plot(lor, xlab="Department", ylab="Log Odds Ratio (Admit | Gender)")


###################################################
### code chunk number 45: oddsratio
###################################################
plot(lor, xlab="Department", ylab="Log Odds Ratio (Admit | Gender)")


###################################################
### code chunk number 46: jobsat
###################################################
JobSat


###################################################
### code chunk number 47: cmh1
###################################################
CMHtest(JobSat, rscores=c(7.5,20,32.5,60))


###################################################
### code chunk number 48: assoc1
###################################################
assocstats(GSStab)


###################################################
### code chunk number 49: gamma
###################################################
GKgamma(JobSat)


###################################################
### code chunk number 50: kappa
###################################################
(K <- Kappa(SexualFun))
confint(K)


###################################################
### code chunk number 51: agreesex
###################################################
agree <- agreementplot(SexualFun, main="Is sex fun?")
unlist(agree)


###################################################
### code chunk number 52: ca1
###################################################
library(ca)
ca(HairEye)


###################################################
### code chunk number 53: ca-haireye0 (eval = FALSE)
###################################################
## plot(ca(HairEye), main="Hair Color and Eye Color")
## title(xlab="Dim 1 (89.4%)", ylab="Dim 2 (9.5%)")


###################################################
### code chunk number 54: ca-haireye
###################################################
plot(ca(HairEye), main="Hair Color and Eye Color")
title(xlab="Dim 1 (89.4%)", ylab="Dim 2 (9.5%)")


###################################################
### code chunk number 55: loglm-hec1
###################################################
library(MASS)
## Independence model of hair and eye color and sex.  
hec.1 <- loglm(~Hair+Eye+Sex, data=HairEyeColor)
hec.1


###################################################
### code chunk number 56: loglm-hec2
###################################################
## Conditional independence
hec.2 <- loglm(~(Hair + Eye) * Sex, data=HairEyeColor)
hec.2


###################################################
### code chunk number 57: loglm-hec3
###################################################
## Joint independence model.  
hec.3 <- loglm(~Hair*Eye + Sex, data=HairEyeColor)
hec.3


###################################################
### code chunk number 58: loglm-anova
###################################################
anova(hec.1, hec.2, hec.3)


###################################################
### code chunk number 59: mental1
###################################################
str(Mental)
xtabs(Freq ~ mental+ses, data=Mental)   # display the frequency table


###################################################
### code chunk number 60: mental2
###################################################
indep <- glm(Freq ~ mental + ses, family = poisson, data = Mental)  # independence model


###################################################
### code chunk number 61: mental3
###################################################
# Use integer scores for rows/cols 
Cscore <- as.numeric(Mental$ses)
Rscore <- as.numeric(Mental$mental)	


###################################################
### code chunk number 62: mental4
###################################################
# column effects model (ses)
coleff <- glm(Freq ~ mental + ses + Rscore:ses, family = poisson, data = Mental)

# row effects model (mental)
roweff <- glm(Freq ~ mental + ses + mental:Cscore, family = poisson, data = Mental)

# linear x linear association
linlin <- glm(Freq ~ mental + ses + Rscore:Cscore, family = poisson, data = Mental)


###################################################
### code chunk number 63: mental4a
###################################################
# compare models using AIC, BIC, etc
vcdExtra::LRstats(glmlist(indep, roweff, coleff, linlin))


###################################################
### code chunk number 64: mental5
###################################################
anova(indep, linlin, coleff, test="Chisq")	
anova(indep, linlin, roweff, test="Chisq")	


###################################################
### code chunk number 65: mental6
###################################################
CMHtest(xtabs(Freq~ses+mental, data=Mental))


###################################################
### code chunk number 66: mental7
###################################################
RC1 <- gnm(Freq ~ mental + ses + Mult(mental,ses), data=Mental, 
             family=poisson, , verbose=FALSE)
RC2 <- gnm(Freq ~ mental+ses + instances(Mult(mental,ses),2), data=Mental, 
             family=poisson, verbose=FALSE)
anova(indep, RC1, RC2, test="Chisq")


###################################################
### code chunk number 67: Arthritis1 (eval = FALSE)
###################################################
## mosaic(art, gp = shading_max, split_vertical = TRUE, 
##        main="Arthritis: [Treatment] [Improved]")


###################################################
### code chunk number 68: art1
###################################################
summary(art)


###################################################
### code chunk number 69: art21
###################################################
mosaic(art, gp = shading_max, split_vertical = TRUE, 
       main="Arthritis: gp = shading_max")


###################################################
### code chunk number 70: art22
###################################################
mosaic(art, gp = shading_Friendly, split_vertical = TRUE, 
       main="Arthritis: gp = shading_Friendly")


###################################################
### code chunk number 71: hec-mosaic (eval = FALSE)
###################################################
## # mosaic plots, using plot.loglm() method
## plot(hec.1, main="model: [Hair][Eye][Sex]")
## plot(hec.2, main="model: [HairSex][EyeSex]")
## plot(hec.3, main="model: [HairEye][Sex]")


###################################################
### code chunk number 72: hec1
###################################################
plot(hec.1, main="model: [Hair][Eye][Sex]")


###################################################
### code chunk number 73: hec2
###################################################
plot(hec.2, main="model: [HairSex][EyeSex]")


###################################################
### code chunk number 74: hec3
###################################################
plot(hec.3, main="model: [HairSex][EyeSex]")


###################################################
### code chunk number 75: TV-mosaic0 (eval = FALSE)
###################################################
## dimnames(TV2)$Time <- c("8", "9", "10")     # re-level for mosaic display
## mosaic(~ Day + Network + Time, data=TV2, expected=~Day:Time + Network, 
##          legend=FALSE, gp=shading_Friendly)


###################################################
### code chunk number 76: TV-mosaic
###################################################
dimnames(TV2)$Time <- c("8", "9", "10")     # re-level for mosaic display
mosaic(~ Day + Network + Time, data=TV2, expected=~Day:Time + Network, 
         legend=FALSE, gp=shading_Friendly)


###################################################
### code chunk number 77: TV-mosaic1 (eval = FALSE)
###################################################
## mosaic(~ Day + Network + Time, data=TV2, 
##          expected=~Day:Time + Day:Network + Time:Network, 
##          legend=FALSE, gp=shading_Friendly)


###################################################
### code chunk number 78: mental-plots (eval = FALSE)
###################################################
## long.labels <- list(set_varnames = c(mental="Mental Health Status", ses="Parent SES"))
## mosaic(indep, ~ses+mental, residuals_type="rstandard",  
##      labeling_args = long.labels, labeling=labeling_residuals,
##      main="Mental health data: Independence")
## 
## mosaic(linlin, ~ses+mental, residuals_type="rstandard", 
##      labeling_args = long.labels, labeling=labeling_residuals, suppress=1, 
##      gp=shading_Friendly, main="Mental health data: Linear x Linear")


###################################################
### code chunk number 79: mental-plots1
###################################################
long.labels <- list(set_varnames = c(mental="Mental Health Status", ses="Parent SES"))
mosaic(indep, ~ses+mental, residuals_type="rstandard",  
     labeling_args = long.labels, labeling=labeling_residuals,
     main="Mental health data: Independence")


###################################################
### code chunk number 80: mental-plots2
###################################################
long.labels <- list(set_varnames = c(mental="Mental Health Status", ses="Parent SES"))
mosaic(linlin, ~ses+mental, residuals_type="rstandard", 
     labeling_args = long.labels, labeling=labeling_residuals, suppress=1, 
     gp=shading_Friendly, main="Mental health data: Linear x Linear")


###################################################
### code chunk number 81: mental-RC1 (eval = FALSE)
###################################################
## Mental$mental <- C(Mental$mental, treatment)
## Mental$ses <- C(Mental$ses, treatment)
## RC1model <- gnm(Freq ~ mental + ses + Mult(mental, ses),
##                 family = poisson, data = Mental)
## mosaic(RC1model, residuals_type="rstandard", labeling_args = long.labels, 
##        labeling=labeling_residuals, suppress=1, gp=shading_Friendly,
##        main="Mental health data: RC(1) model")


###################################################
### code chunk number 82: reorder3a
###################################################
UCB <- aperm(UCBAdmissions, c(2, 1, 3))
names(dimnames(UCB)) <- c("Sex", "Admit?", "Department")
dimnames(UCB)[[2]] <- c("Yes", "No")


###################################################
### code chunk number 83: reorder3b
###################################################
vnames <- list(set_varnames = c(Admit="Admission", Gender="Sex", Dept="Department"))

lnames <- list(Admit = c("Yes", "No"),
              Gender = c("Males", "Females"),
              Dept = LETTERS[1:6])

mosaic(UCBAdmissions, labeling_args=vnames, set_labels=lnames)


###################################################
### code chunk number 84: TV-mosaic2
###################################################
dimnames(TV2)$Time <- c("8", "9", "10")     # re-level for mosaic display
mosaic(~ Day + Network + Time, data=TV2, expected=~Day:Time + Network, 
         legend=FALSE, gp=shading_Friendly, 
         labeling_args=list(abbreviate=c(Day=3)) )


###################################################
### code chunk number 85: spine1
###################################################
(spine(Improved ~ Age, data = Arthritis, breaks = 3))
(spine(Improved ~ Age, data = Arthritis, breaks = "Scott"))


###################################################
### code chunk number 86: spine2
###################################################
(spine(Improved ~ Age, data = Arthritis, breaks = 3))


###################################################
### code chunk number 87: spine3
###################################################
(spine(Improved ~ Age, data = Arthritis, breaks = "Scott"))


###################################################
### code chunk number 88: cdplot
###################################################
cdplot(Improved ~ Age, data = Arthritis)
with(Arthritis, rug(jitter(Age), col="white", quiet=TRUE))


###################################################
### code chunk number 89: cdplot1
###################################################
cdplot(Improved ~ Age, data = Arthritis)
with(Arthritis, rug(jitter(Age), col="white", quiet=TRUE))


###################################################
### code chunk number 90: donner1
###################################################
data(Donner, package="vcdExtra")
str(Donner)


###################################################
### code chunk number 91: donner2a (eval = FALSE)
###################################################
## # separate linear fits on age for M/F
## ggplot(Donner, aes(age, survived, color = sex)) +
##   geom_point(position = position_jitter(height = 0.02, width = 0)) +
##   stat_smooth(method = "glm", method.args = list(family = binomial), formula = y ~ x,
##            alpha = 0.2, size=2, aes(fill = sex))


###################################################
### code chunk number 92: donner2b (eval = FALSE)
###################################################
## # separate quadratics
## ggplot(Donner, aes(age, survived, color = sex)) +
##   geom_point(position = position_jitter(height = 0.02, width = 0)) +
##   stat_smooth(method = "glm", method.args = list(family = binomial), formula = y ~ poly(x,2),
##            alpha = 0.2, size=2, aes(fill = sex))


###################################################
### code chunk number 93: donner3a
###################################################
ggplot(Donner, aes(age, survived, color = sex)) +
  geom_point(position = position_jitter(height = 0.02, width = 0)) +
  stat_smooth(method = "glm", method.args = list(family = binomial), formula = y ~ x,
           alpha = 0.2, size=2, aes(fill = sex))



###################################################
### code chunk number 94: donner3b
###################################################
# separate quadratics
ggplot(Donner, aes(age, survived, color = sex)) +
  geom_point(position = position_jitter(height = 0.02, width = 0)) +
  stat_smooth(method = "glm", method.args = list(family = binomial), formula = y ~ poly(x,2),
           alpha = 0.2, size=2, aes(fill = sex))



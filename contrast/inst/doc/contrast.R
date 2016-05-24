### R code from vignette source 'contrast.Rnw'

###################################################
### code chunk number 1: loadLib
###################################################
library(contrast)
library(MASS)
library(nlme)
library(lattice)
library(Hmisc)
library(geepack)
options(useFancyQuotes = FALSE, width = 80)


###################################################
### code chunk number 2: example1Data
###################################################
example1 <- expand.grid(
   tissue = c("Whole Blood", "PBMC"),
   trmt = c("Vehicle", "Compound1", "Compound2"),   
   rep = 1:3)
example1 <- data.frame(
   expression = 
   c(8.029, 7.650, 7.304, 7.240, 7.322, 7.727, 
     7.554, 7.198, 7.575, 7.399, 8.211, 7.321,
     7.548, 7.538, 7.565, 7.304, 7.998, 8.067, 
     7.755, 7.415, 7.850, 7.217, 8.030, 7.432),
   diet =
   c("high fat",  "high fat",  "high fat",   "low fat",  
      "low fat",   "low fat",   "low fat",   "low fat", 
      "low fat",   "low fat",  "high fat",   "low fat",  
      "low fat",   "low fat",   "low fat",   "low fat", 
     "high fat",  "high fat",  "high fat",  "high fat",  
     "high fat",  "high fat",  "high fat",  "high fat"),
   group =  
   c("treatment",     "vehicle",     "vehicle",  "vehicle",   
      "vehicle",   "treatment",   "treatment",  "vehicle",
      "treatment",    "vehicle",   "treatment",  "vehicle",    
      "treatment",  "treatment",   "treatment",  "vehicle",    
      "treatment",  "treatment",     "vehicle",  "vehicle",    
      "treatment",    "vehicle",   "treatment",  "vehicle"))
      
   designSummary1 <- as.data.frame(table(example1$diet, example1$group))
   names(designSummary1)[1:2] <- c("Diet", "Group")


###################################################
### code chunk number 3: designSummary1
###################################################
latex(
      designSummary1, 
      ctable = TRUE, 
      rowname=NULL, 
      caption = "A summary of the diet experimental design",
      where = "h",
      file = "")


###################################################
### code chunk number 4: example1Plot
###################################################
plotTheme <- col.whitebg()
plotTheme$superpose.line$lty <- rep(1, 8)

trellis.par.set(plotTheme, warn = FALSE)

out <- stripplot(
                 expression ~ group, 
                 example1,  
                 groups = diet,
                 as.Table = TRUE,
                 type = c("p", "a", "g"),
                 auto.key = list(columns = 2),
                 ylab = "Expression")  
print(out)



###################################################
### code chunk number 5: example1LinearMod
###################################################
lmFit1 <- lm(expression ~ (group + diet)^2, data = example1)
summary(lmFit1)


###################################################
### code chunk number 6: example1Contrast
###################################################
highFatDiff <- contrast(lmFit1, 
                        list(diet = "low fat", group = "vehicle"),
                        list(diet = "low fat", group = "treatment"))
print(highFatDiff, X = TRUE)


###################################################
### code chunk number 7: example1ContrastStat
###################################################
basicTestStat <- highFatDiff$testStat


###################################################
### code chunk number 8: eachTest
###################################################
eachTrmtEffect <- contrast(lmFit1, 
                           list(diet = levels(example1$diet), group = "vehicle"),
                           list(diet = levels(example1$diet), group = "treatment"))  
print(eachTrmtEffect, X = TRUE)


###################################################
### code chunk number 9: meanEffect
###################################################
meanTrmtEffect <- contrast(lmFit1, 
                           list(diet = levels(example1$diet), group = "vehicle"),
                           list(diet = levels(example1$diet), group = "treatment"),
                           type = "average")  
  
print(meanTrmtEffect, X = TRUE)


###################################################
### code chunk number 10: example1Sand
###################################################
highFatDiffSAND <- contrast(lmFit1, 
                            list(diet = "low fat", group = "vehicle"),
                            list(diet = "low fat", group = "treatment"),
                            covType = "HC3")
print(highFatDiffSAND)


###################################################
### code chunk number 11: lmTest
###################################################
lmTest <- contrast(lmFit1, 
                   list(diet = "high fat", group = "vehicle"),
                   list(diet = "high fat", group = "treatment"))
if(abs(lmTest$Pvalue- summary(lmFit1)$coefficients[2,4]) > 0.0001)
stop("lm test case failed")      


###################################################
### code chunk number 12: example1GenLinearMod
###################################################
glmFit1 <- glm(2^expression ~ (group + diet)^2, 
               data = example1, 
               family = gaussian(link = "log"))
summary(glmFit1)
highFatDiff <- contrast(glmFit1, 
                        list(diet = "low fat", group = "vehicle"),
                        list(diet = "low fat", group = "treatment"))
print(highFatDiff, X = TRUE)


###################################################
### code chunk number 13: glmTest
###################################################
glmTest <- contrast(glmFit1, 
                    list(diet = "high fat", group = "vehicle"),
                    list(diet = "high fat", group = "treatment"))
if(abs(glmTest$Pvalue - summary(glmFit1)$coefficients[2,4]) > 0.0001)
  stop("glm test case failed")       


###################################################
### code chunk number 14: example2Data
###################################################
 example2 <- data.frame(
   subject =
      c("donor1", "donor1", "donor1", "donor1", "donor1", "donor1", "donor1", 
        "donor2", "donor2", "donor2", "donor3", "donor2", "donor2", "donor2", 
        "donor2", "donor3", "donor3", "donor3", "donor3", "donor3", "donor3", 
        "donor3", "donor1"),
   config = 
      c("A", "D", "C", "D", "C", "D", "C", 
        "B", "A", "D", "B", "C", "D", "D", 
        "C", "A", "D", "C", "D", "D", "C", 
        "C", "B"),
   day = c(1, 1, 1, 2, 2, 4, 4, 1, 1, 1, 1, 1, 2, 4, 4, 1, 1, 1, 2, 4, 2, 4, 1),
   expression = 
      c(9.364, 9.533, 9.267, 9.407, 9.463, 9.419, 9.507, 
        9.226, 9.153, 9.112, 9.657, 8.973, 9.096, 9.232, 
        9.261, 9.383, 9.365, 9.242, 9.361, 9.345, 9.444, 
        9.551, 9.601))
   example2Table <- as.data.frame(table(example2$day, example2$config))
   names(example2Table) <- c("Day", "Configuration", "Number of Donors")


###################################################
### code chunk number 15: designSummary2
###################################################
latex(
      example2Table, 
      ctable = TRUE, 
      rowname=NULL, 
      caption = "A summary of the stem cell experimental design",
      where = "h",
      file = "")


###################################################
### code chunk number 16: design2factor
###################################################
example2$group <- factor(paste(example2$day, ":", example2$config, sep = ""))
print(table(example2$group))


###################################################
### code chunk number 17: design2gls
###################################################
glsFit <-  gls(expression ~ group, 
               data = example2, 
               corCompSymm(form = ~ 1 | subject))
summary(glsFit)


###################################################
### code chunk number 18: design2glsCont
###################################################
print(
      contrast(glsFit, 
               list(group = "4:C"),
               list(group = "4:D")),
      X = TRUE)     


###################################################
### code chunk number 19: example2Plot
###################################################
out2 <- stripplot(
                  expression ~ factor(day), 
                  example2,  
                  groups = factor(config),
                  as.Table = TRUE,
                  type = c("p", "a", "g"),
                  xlab = "Day",
                  ylab = "Log2 Expression",
                  auto.key = list(columns = 4))   
print(out2)


###################################################
### code chunk number 20: glsTest
###################################################
glsTest <- contrast(glsFit, 
                    list(group = "1:A"),
                    list(group = "4:D"))
if(abs(glsTest$Pvalue - summary(glsFit)$tTable[8,4])  > 0.0001)
  stop("gls test case failed")      


###################################################
### code chunk number 21: design2lme
###################################################
lmeFit <-  lme(expression ~ group, 
               data = example2, 
               random = ~1|subject)
summary(lmeFit)

print(
      contrast(lmeFit, 
               list(group = "4:C"),
               list(group = "4:D")),
      X = TRUE)        


###################################################
### code chunk number 22: lmeTest
###################################################
lmeTest <- contrast(lmeFit, 
                    list(group = "1:A"),
                    list(group = "4:D"))
if(abs(lmeTest$Pvalue - summary(lmeFit)$tTable[8,5])  > 0.0001)
  stop("lme test case failed")      


###################################################
### code chunk number 23: design2lme
###################################################
geeFit <-  geese(2^expression ~ group,
                 data = example2,
                 id = subject,
                 family = gaussian(link = "log"),
                 corstr = "exchangeable")
summary(geeFit)

print(
      contrast(geeFit, 
               list(group = "4:C"),
               list(group = "4:D")),
      X = TRUE)   


###################################################
### code chunk number 24: geeTest
###################################################
geeTest <- contrast(geeFit, 
                    list(group = "1:A"),
                    list(group = "4:D"))
if(abs(geeTest$Pvalue - summary(geeFit)$mean[8,4])  > 0.0001)
  stop("gee test case failed")      


###################################################
### code chunk number 25: ex1FC
###################################################
eachTrmtEffect <- contrast(lmFit1, 
                           list(diet = levels(example1$diet), group = "vehicle"),
                           list(diet = levels(example1$diet), group = "treatment"),
                           fcfunc = function(u)2^u)  
print(eachTrmtEffect, X = TRUE)
eachTrmtEffect$foldChange



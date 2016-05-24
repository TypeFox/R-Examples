### R code from vignette source 'introduction.Rnw'

###################################################
### code chunk number 1: introduction.Rnw:40-43
###################################################
require(hamlet)
data(vcapwide)
vcapwide[1:2,]


###################################################
### code chunk number 2: introduction.Rnw:48-50
###################################################
data(vcaplong)
vcaplong[1:3,]


###################################################
### code chunk number 3: introduction.Rnw:70-72
###################################################
ex <- read.table(file="example.csv", sep=";", dec=",", stringsAsFactors=F, header=T)
ex


###################################################
### code chunk number 4: introduction.Rnw:101-104
###################################################
d <- dist(ex[,2:4]) # By default Euclidean distance
d <- as.matrix(d)
d <- round(d, 2) # distance matrix d


###################################################
### code chunk number 5: introduction.Rnw:107-109
###################################################
require(xtable)
print(xtable(d, caption="Euclidean distance matrix D for 18 animals", label="tab:eucld"), scalebox=0.7)


###################################################
### code chunk number 6: introduction.Rnw:125-129
###################################################
sol <- match.bb(d, g=3)
submatches <- paste("Submatch_", LETTERS[1:6][sol$solution], sep="")
names(submatches) <- names(sol$solution)
submatches


###################################################
### code chunk number 7: introduction.Rnw:146-149
###################################################
ex[,"Submatch"] <- submatches
set.seed(1) # for reproducibility
ex[,"AllocatedGroups"] <- match.allocate(ex[,"Submatch"])


###################################################
### code chunk number 8: introduction.Rnw:152-154
###################################################
require(xtable)
print(xtable(ex, caption="The result table in variable \\texttt{ex} after performing the optimal matching and allocation.", label="tab:extable"), scalebox=0.8)


###################################################
### code chunk number 9: introduction.Rnw:170-172
###################################################
boxplot(PSA.week.10..ug.l. ~ AllocatedGroups, data = ex, range=0, 
xlab="Group", ylab="PSA week 10 ul/g")


###################################################
### code chunk number 10: introduction.Rnw:181-182
###################################################
mixplot(ex[,2:5], pch=16)


###################################################
### code chunk number 11: introduction.Rnw:189-190
###################################################
mixplot(ex[,c(2:4,6)], pch=16)


###################################################
### code chunk number 12: introduction.Rnw:200-201
###################################################
heatmap(d)


###################################################
### code chunk number 13: introduction.Rnw:219-224
###################################################
veh <- vcapwide[vcapwide[,"Group"]=="Vehicle",
	c("Submatch","PSAWeek10","BWWeek10","PSAWeek14")]
mdv <- vcapwide[vcapwide[,"Group"]=="MDV",
	c("Submatch","PSAWeek10","BWWeek10","PSAWeek14")]
t.test(veh[,"PSAWeek14"], mdv[,"PSAWeek14"])


###################################################
### code chunk number 14: introduction.Rnw:229-231
###################################################
veh <- veh[order(veh[,"Submatch"]),]
mdv <- mdv[order(mdv[,"Submatch"]),]


###################################################
### code chunk number 15: introduction.Rnw:234-237
###################################################
mat1 <- cbind(Veh.PSAWeek10 = veh[,"PSAWeek10"], MDV.PSAWeek10 = mdv[,"PSAWeek10"])
rownames(mat1) <- veh[,"Submatch"]
print(xtable(mat1, caption="Submatches in the real VCaP experiment, per PSA at week 10 in tumors allocated to the Vehicle and MDV groups", label="tab:mat1"), scalebox=0.9)


###################################################
### code chunk number 16: introduction.Rnw:240-243
###################################################
mat2 <- cbind(Veh.BWWeek10 = veh[,"BWWeek10"], MDV.BWWeek10 = mdv[,"BWWeek10"])
rownames(mat2) <- veh[,"Submatch"]
print(xtable(mat2, caption="Submatches in the real VCaP experiment, per body weight at week 10 in tumors allocated to the Vehicle and MDV groups", label="tab:mat2"), scalebox=0.9)


###################################################
### code chunk number 17: introduction.Rnw:249-250
###################################################
t.test(veh[,"PSAWeek14"], mdv[,"PSAWeek14"], paired=TRUE)



### R code from vignette source 'arules.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: arules.Rnw:73-76
###################################################
options(width = 75)
### for sampling
set.seed <- 1234


###################################################
### code chunk number 2: arules.Rnw:1141-1142
###################################################
library("arules")


###################################################
### code chunk number 3: epub1
###################################################
data("Epub")
Epub


###################################################
### code chunk number 4: epub2
###################################################
summary(Epub)


###################################################
### code chunk number 5: arules.Rnw:1167-1169
###################################################
year <- strftime(as.POSIXlt(transactionInfo(Epub)[["TimeStamp"]]), "%Y")
table(year)


###################################################
### code chunk number 6: arules.Rnw:1177-1180
###################################################
Epub2003 <- Epub[year == "2003"]
length(Epub2003)
image(Epub2003)


###################################################
### code chunk number 7: epub
###################################################
print(image(Epub2003))


###################################################
### code chunk number 8: arules.Rnw:1206-1207
###################################################
transactionInfo(Epub2003[size(Epub2003) > 20])


###################################################
### code chunk number 9: arules.Rnw:1219-1220
###################################################
inspect(Epub2003[1:5])


###################################################
### code chunk number 10: arules.Rnw:1226-1227
###################################################
as(Epub2003[1:5], "list")


###################################################
### code chunk number 11: arules.Rnw:1233-1235
###################################################
EpubTidLists <- as(Epub, "tidLists")
EpubTidLists


###################################################
### code chunk number 12: arules.Rnw:1242-1243
###################################################
as(EpubTidLists[1:3], "list") 


###################################################
### code chunk number 13: data
###################################################
data("AdultUCI")
dim(AdultUCI)
AdultUCI[1:2,]


###################################################
### code chunk number 14: arules.Rnw:1290-1292
###################################################
AdultUCI[["fnlwgt"]] <- NULL
AdultUCI[["education-num"]] <- NULL


###################################################
### code chunk number 15: arules.Rnw:1308-1323
###################################################
AdultUCI[[ "age"]] <- ordered(cut(AdultUCI[[ "age"]], c(15,25,45,65,100)),
    labels = c("Young", "Middle-aged", "Senior", "Old"))

AdultUCI[[ "hours-per-week"]] <- ordered(cut(AdultUCI[[ "hours-per-week"]],
      c(0,25,40,60,168)),
    labels = c("Part-time", "Full-time", "Over-time", "Workaholic"))
			    
AdultUCI[[ "capital-gain"]] <- ordered(cut(AdultUCI[[ "capital-gain"]],
      c(-Inf,0,median(AdultUCI[[ "capital-gain"]][AdultUCI[[ "capital-gain"]]>0]),Inf)),
    labels = c("None", "Low", "High"))

AdultUCI[[ "capital-loss"]] <- ordered(cut(AdultUCI[[ "capital-loss"]],
      c(-Inf,0,
	median(AdultUCI[[ "capital-loss"]][AdultUCI[[ "capital-loss"]]>0]),Inf)),
    labels = c("none", "low", "high"))


###################################################
### code chunk number 16: coerce
###################################################
Adult <- as(AdultUCI, "transactions")
Adult


###################################################
### code chunk number 17: summary
###################################################
summary(Adult)


###################################################
### code chunk number 18: itemFrequencyPlot (eval = FALSE)
###################################################
## itemFrequencyPlot(Adult, support = 0.1, cex.names=0.8)


###################################################
### code chunk number 19: arules.Rnw:1365-1366
###################################################
itemFrequencyPlot(Adult, support = 0.1, cex.names=0.8)


###################################################
### code chunk number 20: apriori
###################################################
rules <- apriori(Adult, 
                 parameter = list(support = 0.01, confidence = 0.6))
rules


###################################################
### code chunk number 21: summary
###################################################
summary(rules)


###################################################
### code chunk number 22: rules
###################################################
rulesIncomeSmall <- subset(rules, subset = rhs %in% "income=small" & lift > 1.2)
rulesIncomeLarge <- subset(rules, subset = rhs %in% "income=large" & lift > 1.2)


###################################################
### code chunk number 23: subset
###################################################
inspect(head(rulesIncomeSmall, n = 3, by = "confidence"))
inspect(head(rulesIncomeLarge, n = 3, by = "confidence"))


###################################################
### code chunk number 24: write_rules (eval = FALSE)
###################################################
## write(rulesIncomeSmall, file = "data.csv", sep = ",", col.names = NA)


###################################################
### code chunk number 25: pmml (eval = FALSE)
###################################################
## write.PMML(rulesIncomeSmall, file = "data.xml")


###################################################
### code chunk number 26: arules.Rnw:1490-1493
###################################################
data("Adult")
fsets <- eclat(Adult, parameter = list(support = 0.05), 
	control = list(verbose=FALSE))


###################################################
### code chunk number 27: arules.Rnw:1501-1508
###################################################
singleItems <- fsets[size(items(fsets)) == 1]

## Get the col numbers we have support for
singleSupport <- quality(singleItems)$support
names(singleSupport) <- unlist(LIST(items(singleItems),
	    decode = FALSE))
head(singleSupport, n = 5)


###################################################
### code chunk number 28: arules.Rnw:1517-1524
###################################################
itemsetList <- LIST(items(fsets), decode = FALSE)

allConfidence <- quality(fsets)$support / 
    sapply(itemsetList, function(x) 
    max(singleSupport[as.character(x)]))

quality(fsets) <- cbind(quality(fsets), allConfidence)


###################################################
### code chunk number 29: arules.Rnw:1528-1529
###################################################
summary(fsets)


###################################################
### code chunk number 30: arules.Rnw:1537-1540
###################################################
fsetsEducation <- subset(fsets, subset = items %pin% "education")
inspect(sort(fsetsEducation[size(fsetsEducation)>1], 
	by = "allConfidence")[1 : 3])


###################################################
### code chunk number 31: arules.Rnw:1553-1555
###################################################
data("Adult")
Adult


###################################################
### code chunk number 32: arules.Rnw:1564-1570
###################################################
supp <- 0.05
epsilon <- 0.1
c <- 0.1

n <- -2 * log(c)/ (supp * epsilon^2)
n


###################################################
### code chunk number 33: arules.Rnw:1579-1580
###################################################
AdultSample <- sample(Adult, n, replace = TRUE)


###################################################
### code chunk number 34: itemFrequencyPlot2 (eval = FALSE)
###################################################
## itemFrequencyPlot(AdultSample, population = Adult, support = supp,
##     cex.names = 0.7)


###################################################
### code chunk number 35: arules.Rnw:1598-1599
###################################################
itemFrequencyPlot(AdultSample, population = Adult, support = supp,
    cex.names = 0.7)


###################################################
### code chunk number 36: itemFrequencyPlot3 (eval = FALSE)
###################################################
## itemFrequencyPlot(AdultSample, population = Adult, 
##     support = supp, lift = TRUE, 
##     cex.names = 0.9)


###################################################
### code chunk number 37: arules.Rnw:1627-1628
###################################################
itemFrequencyPlot(AdultSample, population = Adult, 
    support = supp, lift = TRUE, 
    cex.names = 0.9)


###################################################
### code chunk number 38: arules.Rnw:1639-1646
###################################################
time <- system.time(itemsets <- eclat(Adult, 
    parameter = list(support = supp), control = list(verbose = FALSE)))
time

timeSample <- system.time(itemsetsSample <- eclat(AdultSample, 
    parameter = list(support = supp), control = list(verbose = FALSE)))
timeSample


###################################################
### code chunk number 39: arules.Rnw:1655-1657
###################################################
# speed up
time[1] / timeSample[1]


###################################################
### code chunk number 40: arules.Rnw:1665-1667
###################################################
itemsets
itemsetsSample


###################################################
### code chunk number 41: arules.Rnw:1675-1678
###################################################
match <- match(itemsets, itemsetsSample, nomatch = 0)
## remove no matches
sum(match > 0) / length(itemsets)


###################################################
### code chunk number 42: arules.Rnw:1687-1689
###################################################
summary(quality(itemsets[match == 0])$support)
summary(quality(itemsetsSample[-match])$support)


###################################################
### code chunk number 43: arules.Rnw:1698-1704
###################################################
supportItemsets <- quality(itemsets[which(match > 0)])$support
supportSample <- quality(itemsetsSample[match])$support

accuracy <- 1 - abs(supportSample - supportItemsets) / supportItemsets

summary(accuracy)



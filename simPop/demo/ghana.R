# -----------------------
# Authors: Andreas Alfons
#          K.U.Leuven
# -----------------------

## initializations
library("simPopulation")
data("ghanaS")
seed <- 1234

## simulate household structure
ghanaP <- simStructure(ghanaS, hid = "hhid", w = "weight", strata = "region", 
    additional = c("age", "sex", "relate"), keep = FALSE)

## categorize age
breaks <- c(min(ghanaS$age), seq(6, 18, 2), seq(20, 80, 5), max(ghanaS$age))
ghanaS$ageCat <- cut(ghanaS$age, breaks = breaks, include.lowest = TRUE)
ghanaP$ageCat <- cut(ghanaP$age, breaks = breaks, include.lowest = TRUE)

## simulate categorical variables conditional on household head
basic <- c("ageCat", "sex", "hsize")
additional <- c("nation", "ethnic", "religion")
ghanaP <- simRelation(ghanaS, ghanaP, hid = "hhid", w = "weight", 
    strata = "region", basic = basic, additional = additional)

## simulate additional categorical variables
basic <- c(basic, additional)
additional <- c("highest_degree", "occupation")
# structural zeros for "highest_degree"
hd14 <- c("bece", "mslc", "ssce", "other", "tech/prof cert")
hd18 <- c("hnd", "gce 'a' level", "gce 'o' level", "tech/prof dip", "voc/comm")
hd21 <- c("bachelor", "teacher trng a", "teacher trng b")
limitHD <- list(ageCat=list("[0,6]" = "none", "(6,8]" = "none", 
        "(8,10]" = "none", "(10,12]" = "none", "(12,14]" = c("none", hd14), 
        "(14,16]" = c("none", hd14), "(16,18]" = c("none", hd14, hd18), 
        "(18,20]" = c("none", hd14, hd18), 
        "(20,25]" = c("none", hd14, hd18, hd21)))
# structural zeros for "occupation"
o8 <- c("craft and related trades workers", "elementary occupations", 
    "service workers and shop and market sales workers", 
    "skilled agricultural and fishery workers")
o18 <- c("armed forces and other security personnel", "clerks", 
    "plant and machine operators and assemblers", 
    "technicians and associate professionals")
limitO <- list(ageCat=list("[0,6]" = "none", "(6,8]" =c ("none", o8), 
        "(8,10]" = c("none", o8), "(10,12]" = c("none", o8), 
        "(12,14]" = c("none", o8), "(14,16]" = c("none", o8), 
        "(16,18]" = c("none", o8, o18), "(18,20]" = c("none", o8, o18), 
        "(20,25]" = c("none", o8, o18)))
# simulate "highest_degree" and "occupation"
ghanaP <- simCategorical(ghanaS, ghanaP, w = "weight", 
    strata = "region", basic = basic, additional = additional, 
    limit = list(highest_degree = limitHD, occupation = limitO))

## combine larger household sizes for better readable mosaic plots
breaks <- c(0:10, 15, 20, max(ghanaS$hsize))
labs <- c(1:10, "(11,15]", "(15,20]", "(20,29]")
ghanaS$hsizeCat <- cut(ghanaS$hsize, breaks = breaks, labels = labs)
ghanaP$hsizeCat <- cut(as.numeric(as.character(ghanaP$hsize)), 
    breaks = breaks, labels = labs)

## create mosaic plots for sex, region and household size
abbRegion <- c("W", "C", "GA", "V", "E", "A", "BA", "N", "UE", "UW")
abbHsize <- c(1:10, "15", "20", ">")
nam <- c(sex = "Sex", region = "Region", hsizeCat = "Household size")
lab <- labeling_border(set_labels = list(region = abbRegion, 
        hsizeCat = abbHsize), set_varnames = nam)
spMosaic(c("sex", "region", "hsizeCat"), "weight", ghanaS, ghanaP, 
    labeling = lab)

## create mosaic plots for sex, ethnicity and occupation
abbEthnic <- c("A", "O", "E", "GD", "Gs", "Gn", "Gm", "M", "MD")
abbOccupation <- c("F", "Cl", "Cr", "E", "O", "N", "M", "P", "S", "A", "T")
nam <- c(sex = "Sex", ethnic = "Ethnicity", occupation = "Occupation")
lab <- labeling_border(set_labels = list(ethnic = abbEthnic, 
        occupation = abbOccupation), set_varnames = nam)
spMosaic(c("sex", "ethnic", "occupation"), "weight", ghanaS, ghanaP, 
    labeling = lab)

## simulate annual income
basic <- c(basic, additional)
probs <- c(0.05, 0.1, seq(from=0.2, to=0.8, by=0.2), 
    0.9, 0.95, 0.975, 0.99, 0.995, 0.999)
# structural zeros for annual income
limitI <- list(occupation = list(none = "0"))
censorI <- list(
    "(1.8e+03,2.6e+03]" = list(ageCat = c("[0,6]","(6,8]","(8,10]")), 
    "(2.6e+03,4.16e+03]" = list(ageCat = c("[0,6]", "(6,8]", "(8,10]", 
            "(10,12]", "(12,14]")), 
    "(4.16e+03,7.68e+03]" = list(ageCat = c("[0,6]", "(6,8]", "(8,10]", 
            "(10,12]", "(12,14]")), 
    "(7.68e+03,1.3e+04]" = list(ageCat = c("[0,6]", "(6,8]", "(8,10]", 
            "(10,12]", "(12,14]", "(14,16]")), 
    "(1.3e+04,3.84e+04]" = list(ageCat = c("[0,6]", "(6,8]", "(8,10]", 
            "(10,12]", "(12,14]", "(14,16]", "(16,18]")), 
    "(3.84e+04,1e+06]" = list(ageCat = c("[0,6]", "(6,8]", "(8,10]", "(10,12]", 
            "(12,14]", "(14,16]", "(16,18]", "(18,20]", "(20,25]"),
        occupation = c("armed forces and other security personnel", 
            "clerks", "elementary occupations")))
# simulate "income"
ghanaP <- simContinuous(ghanaS, ghanaP, w = "weight", strata = "region", 
    basic = basic, additional = "income", upper = 1000000, probs = probs, 
    limit = limitI, censor = censorI, keep = FALSE)
ghanaP$income <- round(ghanaP$income, 2)  # round to two decimal places

## plot cumulative distribution functions of annual income for the 
## main parts of the data
subset <- which(ghanaS[, "income"] > 0)
q <- quantileWt(ghanaS[subset, "income"], ghanaS[subset, "weight"], 
    probs = 0.99)
spCdfplot("income", "weight", dataS = ghanaS, dataP = ghanaP, xlim = c(0, q))

## create box plots of annual income
spBwplot("income", "weight", dataS = ghanaS, dataP = ghanaP, pch = "|")

## create box plots of annual income conditional on region
spBwplot("income", "weight", "region", dataS = ghanaS, dataP = ghanaP, 
    pch = "|", layout = c(1, 10))

## create box plots of annual income conditional on gender
spBwplot("income", "weight", "sex", dataS = ghanaS, dataP = ghanaP, 
    pch = "|", layout = c(1, 2))

## create box plots of annual income conditional on ethnicity
spBwplot("income", "weight", "ethnic", dataS = ghanaS, dataP = ghanaP, 
    pch = "|", layout = c(1, 9))

## create box plots of annual income conditional on occupation
spBwplot("income", "weight", "occupation", dataS = ghanaS, dataP = ghanaP, 
    pch = "|", layout = c(1, 11))

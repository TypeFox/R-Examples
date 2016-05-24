# -----------------------
# Authors: Andreas Alfons
#          K.U.Leuven
# -----------------------

## initializations
library("simPopulation")
data("eusilcS")
seed <- 1234

## wrapper function for EU-SILC
eusilcMP <- simEUSILC(eusilcS, upper = 200000, equidist = FALSE, seed = seed)
eusilcTR <- simEUSILC(eusilcS, method = "twostep", seed = seed)
eusilcTN <- simEUSILC(eusilcS, method = "twostep", residuals = FALSE, 
    seed = seed)

## simulate basic household structure
eusilcP <- simStructure(eusilcS, hid = "db030", w = "db090", 
    strata = "db040", additional = c("age", "rb090"))

## categorize age
breaks <- c(min(eusilcS$age), seq(15, 80, 5), max(eusilcS$age))
eusilcS$ageCat <- as.character(cut(eusilcS$age, 
        breaks = breaks, include.lowest = TRUE))
eusilcP$ageCat <- as.character(cut(eusilcP$age, 
        breaks = breaks, include.lowest = TRUE))

## simulate additional categorical variables
basic <- c("ageCat", "rb090", "hsize")
eusilcP <- simCategorical(eusilcS, eusilcP, w = "rb050", strata = "db040", 
    basic = basic, additional = c("pl030", "pb220a"))

## create mosaic plots of gender, region and household size
abb <- c("B", "LA", "Vi", "C", "St", "UA", "Sa", "T", "Vo")
nam <- c(rb090 = "Gender", db040 = "Region", hsize = "Household size")
lab <- labeling_border(set_labels = list(db040 = abb), set_varnames = nam)
spMosaic(c("rb090", "db040", "hsize"), "rb050", eusilcS, eusilcP, 
    labeling = lab)

## create mosaic plots of gender, economic status and citizenship
nam <- c(rb090 = "Gender", pl030 = "Economic status", pb220a = "Citizenship")
lab <- labeling_border(abbreviate = c(FALSE, FALSE, TRUE), set_varnames = nam)
spMosaic(c("rb090", "pl030", "pb220a"), "rb050", eusilcS, eusilcP, 
    labeling = lab)

## simulate personal net income with different approaches
seedP <- .Random.seed
basic <- c(basic, "pl030", "pb220a")
eusilcMP <- simContinuous(eusilcS, eusilcP, w = "rb050", strata = "db040", 
    basic = basic, additional = "netIncome", upper = 200000, equidist = FALSE, 
    seed = seedP)
seedMP <- .Random.seed
eusilcTR <- simContinuous(eusilcS, eusilcP, w = "rb050", strata = "db040", 
    basic = basic, additional = "netIncome", method = "lm", seed = seedP)
seedTR <- .Random.seed
eusilcTN <- simContinuous(eusilcS, eusilcP, w = "rb050", strata = "db040", 
    basic = basic, additional = "netIncome", method = "lm", residuals = FALSE, 
    seed = seedP)
seedTN <- .Random.seed

## plot cumulative distribution functions of personal net income for the 
## main parts of the data
subset <- which(eusilcS[, "netIncome"] > 0)
q <- quantileWt(eusilcS[subset, "netIncome"], eusilcS[subset, "rb050"], 
    probs = 0.99)
listP <- list(MP = eusilcMP, TR = eusilcTR, TN = eusilcTN)
spCdfplot("netIncome", "rb050", dataS = eusilcS, dataP = listP, xlim = c(0, q))

## create box plots of personal net income
spBwplot("netIncome", "rb050", dataS = eusilcS, dataP = listP, pch = "|")

## create box plots of personal net income conditional on gender
spBwplot("netIncome", "rb050", "rb090", dataS = eusilcS, 
    dataP = listP, pch = "|", layout = c(1, 2))

## create box plots of personal net income conditional on citizenship
spBwplot("netIncome", "rb050", "pb220a", dataS = eusilcS, 
    dataP = listP, pch = "|", layout = c(1, 3))

## create box plots of personal net income conditional on region
spBwplot("netIncome", "rb050", "db040", dataS = eusilcS, 
    dataP = listP, pch = "|", layout = c(1, 9))

## create box plots of personal net income conditional on economic status
spBwplot("netIncome", "rb050", "pl030", dataS = eusilcS, 
    dataP = listP, pch = "|", layout = c(1, 7))

## categorize personal net income
breaks <- getBreaks(eusilcS$netIncome, eusilcS$rb050, 
    upper = Inf, equidist = FALSE)
eusilcS$netIncomeCat <- getCat(eusilcS$netIncome, breaks)
eusilcMP$netIncomeCat <- getCat(eusilcMP$netIncome, breaks)
eusilcTR$netIncomeCat <- getCat(eusilcTR$netIncome, breaks)
eusilcTN$netIncomeCat <- getCat(eusilcTN$netIncome, breaks)

## split personal net income into components
components <- c("py010n", "py050n", "py090n", 
    "py100n", "py110n", "py120n", "py130n", "py140n")
eusilcMP <- simComponents(eusilcS, eusilcMP, w = "rb050", 
    total = "netIncome", components = components, 
    conditional = c("netIncomeCat", "pl030"), seed = seedMP)
eusilcTR <- simComponents(eusilcS, eusilcTR, w = "rb050", 
    total = "netIncome", components = components, 
    conditional = c("netIncomeCat", "pl030"), seed = seedTR)
eusilcTN <- simComponents(eusilcS, eusilcTN, w = "rb050", 
    total = "netIncome", components = components, 
    conditional = c("netIncomeCat", "pl030"), seed = seedTN)

## create box plots of the income components
listP <- list(MP = eusilcMP, TR = eusilcTR, TN = eusilcTN)
spBwplot(components, "rb050", dataS = eusilcS, dataP = listP, 
    pch = "|", minRatio = 0.2, layout = c(2, 4))

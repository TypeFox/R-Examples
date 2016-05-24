#####################
## Mosaic Displays ##
#####################

#########################
## Hair Eye Color Data ##
#########################

data(HairEyeColor)

## Basic Mosaic Display ##

HairEye <- margin.table(HairEyeColor, c(1,2))

mosaic(HairEye, main = "Basic Mosaic Display of Hair Eye Color data")

## Hair Mosaic Display with Pearson residuals ##
Hair <- margin.table(HairEyeColor,1)
Hair
mHair <- as.table(rep(mean(margin.table(HairEyeColor, 1)), 4))
names(mHair) <- names(Hair)
mHair

## Pearson residuals from Equiprobability model ##

resid <- (Hair - mHair) / sqrt(mHair)
resid

## First Step in a Mosaic Display ##

mosaic(Hair, residuals = resid, main = "Hair Color Proportions")

## Hair Eye Mosais Display with Pearson residuals ##

mosaic(HairEye, main = " Hair Eye Color with Pearson residuals")

## Show Pearson Residuals ##

(HairEye - loglin(HairEye, c(1, 2), fit = TRUE)$fit) /
  sqrt(loglin(HairEye, c(1, 2), fit = TRUE)$fit)

###################
## UKSoccer Data ##
###################

data(UKSoccer)

## UKSoccer Mosaic Display ##

mosaic(UKSoccer, main = "UK Soccer Scores")

###############################
## Repeat Victimization Data ##
###############################

data(RepVict)

## mosaic(RepVict[-c(4, 7), -c(4, 7)], main = "Repeat Victimization Data")


##################
## 3-Way Tables ##
##################

## Hair Eye Sex Mosais Display with Pearson residuals ##
mosaic(HairEyeColor, main = "Hair Eye Color Sex" )

mosaic(HairEyeColor, expected = ~ Hair * Eye + Sex,
                main = "Model: (Hair Eye) (Sex)" )

mosaic(HairEyeColor, expected = ~ Hair * Sex + Eye*Sex,
               main = "Model: (Hair Sex) (Eye Sex)")


####################
## Premarital Sex ##
####################

data(PreSex)

## Mosaic display for Gender and Premarital Sexual Expirience ##

## (Gender Pre) ##
mosaic(margin.table(PreSex, c(3, 4)), legend = FALSE,
                main = "Gender and Premarital Sex")

## (Gender Pre)(Extra) ##
mosaic(margin.table(PreSex,c(2,3,4)), legend = FALSE,
                expected = ~ Gender * PremaritalSex + ExtramaritalSex ,
                main = "(PreMaritalSex Gender) (Sex)")

## (Gender Pre Extra)(Marital) ##
mosaic(PreSex,
       expected = ~ Gender * PremaritalSex * ExtramaritalSex + MaritalStatus,
       legend = FALSE,
       main = "(PreMarital ExtraMarital) (MaritalStatus)")

## (GPE)(PEM) ##
mosaic(PreSex,
       expected = ~ Gender * PremaritalSex * ExtramaritalSex
       + MaritalStatus * PremaritalSex * ExtramaritalSex,
       legend = FALSE,
       main = "(G P E) (P E M)")

############################
## Employment Status Data ##
############################

data(Employment)

## Employment Status ##
# mosaic(Employment,
#        expected = ~ LayoffCause * EmploymentLength + EmploymentStatus,
#        main = "(Layoff Employment) + (EmployStatus)")


# mosaic(Employment,
#        expected = ~ LayoffCause * EmploymentLength +
#                 LayoffCause * EmploymentStatus,
#        main = "(Layoff EmpL) (Layoff EmplS)")

# ## Closure ##
# mosaic(Employment[,,1], main = "Layoff : Closure")

# ## Replaced ##
# mosaic(Employment[,,2], main = "Layoff : Replaced")


#####################
## Mosaic Matrices ##
#####################

data(UCBAdmissions)

pairs(PreSex)

pairs(UCBAdmissions)

pairs(UCBAdmissions, type = "conditional")

pairs(UCBAdmissions, type = "pairwise", gp = shading_max)




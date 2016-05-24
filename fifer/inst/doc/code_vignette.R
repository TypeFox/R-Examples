##################################################
##################################################
##################################################
##################################################
		# INSTALLATION
##################################################
##################################################
##################################################
##################################################

### first install devtools so I can use their install_github function
install.packages("devtools")
require(devtools)
### install fifer with the install_github function
devtools::install_github("fifer", username="dustinfife")
require(fifer)

##################################################
##################################################
##################################################
##################################################
		# BROWSE THE PACKAGE
##################################################
##################################################
##################################################
##################################################

### look at functions in the fifer package
contents("fifer")
### pull up vignette
vignette("fifer_package")

##################################################
##################################################
##################################################
##################################################
		# DATA MANIPULATION
##################################################
##################################################
##################################################
##################################################

### load the fake medical data
data(fakeMedicalData)
### look at first several rows of fakeMedicalData
head(fakeMedicalData)


### DEMONSTRATE THE r (range) FUNCTION
#################
#################
#################

## see documentation
?r 
### extract column indices between B_regs_10A and B_regs_9E
bregs = r("B_regs_10A", "B_regs_9E", data.names=names(fakeMedicalData))
bregs

### extract column names instead
bregs = r("B_regs_10A", "B_regs_9E", data.names=names(fakeMedicalData), names=T)
bregs


### DEMONSTRATE THE make.null FUNCTION
#################
#################
#################

### see documentation
?make.null
### extract only demographics columns and bregs columns
newData = make.null("ID", "gender", "ethnicity", "age", "disease",
					bregs, 
					data=fakeMedicalData,
					keep=TRUE)
head(newData)
### we could instead drop everything after bregs
newData2 = make.null(
			r("BCI_10A", "TNF_9E", data.names=names(fakeMedicalData)),
			data=fakeMedicalData, keep=FALSE)					
dim(newData)
dim(newData2)

### DEMONSTRATE THE excelMatch FUNCTION
#################
#################
#################

### see documentation
?excelMatch

#### extact the variable names corresponding to Excel Columns AA, CD, and FF
excel.names = excelMatch("AA", "CD", "FF", names=names(fakeMedicalData))

### or, we can extract the column indices instead
excel.names

### now subset the matrix to just those using make.null


### DEMONSTRATE THE subsetString FUNCTION
#################
#################
#################

### see documentation
?subsetString
#### generate random data (normally this would come from importing a file)
#### print the names (so we can see how messy they are)
head(data)
names(data) = subsetString(names(data), sep=".", position=1)
#### be careful!! The pattern may not be consistent across datasets


##################################################
##################################################
##################################################
##################################################
		# DATA ANALYSIS
##################################################
##################################################
##################################################
##################################################

### DEMONSTRATE THE missing.vals FUNCTION
#################
#################
#################

### see documentation
?missing.vals
?missingVals

### summarize missing values
missing.vals(fakeMedicalData)
missing.vals(na.omit(fakeMedicalData))

### DEMONSTRATE THE demographics FUNCTION
#################
#################
#################

### see documentation
?demographics

### summarize missing values
demographics(disease~age + gender + ethnicity, data=fakeMedicalData, latex=FALSE)

### DEMONSTRATE THE make.formula FUNCTION
#################
#################
#################

### see documentation
?make.formula

#### list all the variables I want to use using the r function
### now write the formula

### DEMONSTRATE THE univariate.tests FUNCTION
#################
#################
#################

### see documentation
?univariate.tests

### compute significance tests for each variable in dataset but the ID column

##################################################
##################################################
##################################################
##################################################
		# PLOTTING
##################################################
##################################################
##################################################
##################################################

#### plot the best five biomarkers
best.five = names(sort(p.adjusted)[1:5])

### prepare the layout
par(mfrow = c(3,2))	## old way
auto.layout(5)		## new way
for (i in 1:length(best.five)){
	### set my favorite defaults
	par1()
	### make the formula
	formula = make.formula(best.five[i], "disease")
	### plot density plots
	densityPlotR(formula, data=fakeMedicalData, main="")
}

### do prism-like plots with significance bars
auto.layout(4)
for (i in 1:4){
	### set different defaults
	par2()
	### make the formula
	formula = make.formula(best.five[i], "disease")
	### plot prism plots
	prism.plots(formula, data=fakeMedicalData)
	### show significance bars
	plotSigBars(formula, data=fakeMedicalData, type="tukey")
}

par1()
#### color code according to disease status
###### make the ... for the densityPlotR
                
#### show off spearman.plot
x = rnorm(100)^2
y = rnorm(100)^2                
### induce a correlation of .6 (approximately) with choleski decomp
cor = matrix(c(1, .6, .6, 1), nrow=2)
skewed.data = cbind(x,y)%*%chol(cor)
names(skewed.data) = c("x", "y")
### show original plot
par2()
plot(skewed.data, xlab="x", ylab="y")
#### now show spearman plot
par2()
spearman.plot(skewed.data, xlab="rank(x)", ylab="rank(y)", pch=16)
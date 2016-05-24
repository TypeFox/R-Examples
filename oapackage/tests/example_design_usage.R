# Analysis of data using a design generated with the Orthogonal Array package

# Alan Vazquez
# Pieter Eendebak
# University of Antwerp
# Department of Engineering Management
# July 2015
#
# For more information see http://pietereendebak.nl/oapackage/index.html
#

# This script uses the following R packages:
# (install with install.packages('packagename'))
library(oapackage)
library(reshape2)


#------------------------------------------------------
# 				Section 1: Introduction
#------------------------------------------------------

# We perform an experiment with N runs and k factors. Main effects as well as two-factor interactions are expected to be active.

N <- 28
k <- 4

# The design used for the analysis is constructed using the Doptimize function with optimization of D+2*Ds

D = Doptimize(N, k, nrestarts=30, alpha1=1, alpha2=2, alpha3=0) 

#------------------------------------------------------
# 				Section 2: Creating data set
#------------------------------------------------------

# For ilustrative purposes, we assume that the true model includes the following effects:
#   B, C, D, BC, BD
#
# We use these effects and the parameter estimates to create a 'true model' and evaluate 
# the performance of a 32-run two-level design constructed using the 'Doptimize' function. 

# The effects of the active factors and two-factor interactions are:
betaparam <- c(10., 1., -2., 3., -.3, -.2)

# The order of the parameter estimates follows the vector:
# c(Intercept, B, C, D, BC, BD)

# The assumed variance for the model residuals is  
sigmasq <- 0.04

design <- as.data.frame(D) # transforming matrix to data frame
print(design) # Print design

# Rename the columns of the design
colnames(design) <- c('A', 'B', 'C', 'D')

# The model matrix is constructed as follows:
formtruem <- ~ B + C + D + B:C + B:D 
modelMat <- model.matrix( formtruem, data = design, contrasts = "contr.helmert" )

y0 <- betaparam%*%t(modelMat) # 
y <- y0 + rnorm(N, mean = 0, sd = sqrt(sigmasq)) # 

# Define data set
data <- data.frame(design, 'Y' = t(y))
print('Simulated Data set')
print(data) # print data set

#------------------------------------------------------
# 				Section 3: Evaluation of the Design
#------------------------------------------------------

MEAlias <- abs(t(D)%*%D)/N # Correlation between main effects matrix
colnames(MEAlias) <- c('A', 'B', 'C', 'D')
rownames(MEAlias) <- c('A', 'B', 'C', 'D')
FullModelMat <- model.matrix( ~(A+B+C+D)^2, data = design)
m=1+k+k*(k-1)/2
TwoFatInt <- FullModelMat[,(2+k):m]
METFIalias <- abs(t(D)%*%TwoFatInt)/N # Correlation between ME and TWOFI
rownames(METFIalias) <- c('A', 'B', 'C', 'D')

# Correlation between main effects (neends reshape2 and ggplot2)
library(ggplot2)
print(MEAlias)
cordata <- melt(MEAlias)
colnames(cordata) <- c('X', 'Y', 'Correlation')
# qplot needs ggplot2
q <- qplot(x=X, y=Y, data=cordata, fill=Correlation, geom="tile") +xlab(' ') + ylab(' ')

# Correlation between main effects and two-factor interactions
print(METFIalias)
cordata <- melt(METFIalias)
colnames(cordata) <- c('X', 'Y', 'Correlation')
q <- qplot(x=X, y=Y, data=cordata, fill=Correlation, geom="tile") +xlab(' ') + ylab(' ')

#------------------------------------------------------
# 				Section 4: Analysis of data
#------------------------------------------------------

# Histogram of Y
hist( data[, 'Y'], main= 'Histogram of Y',xlab = 'Y' )

# Fitting model containing main effects only
lm.ME <- lm( Y ~ A+B+C+D, data = data )
print( summary(lm.ME) )   # print the resulting model

# Fitting model containing main effects and two-factor interactions
lm.Full <- lm( Y ~ (A+B+C+D)^2, data = data )
summary( lm.Full )

# Fitting a reduced model
lm.obj.red <- lm( Y ~ A+B+ C+ D+ A:B + B:C+ B:D, data = data )
summary( lm.obj.red )

# Analysis of residuals

# Correlation
Residuals <- residuals( lm.obj.red )
Predicted <- fitted.values( lm.obj.red )
plot(x= 1:N, Residuals, xlab = 'Run order', ylab = 'Residuals')
abline(h=0)

# Normality
hist( Residuals, main= 'Histogram of Residuals',xlab = 'Residuals' )

# Constant Variance
plot(x= Predicted, y = Residuals, xlab = 'Predicted Values', ylab = 'Residuals')
abline(h=0)


# Coefficients
print(lm.ME)
print(lm.Full)
#lm.Full$coefficients
print(betaparam)



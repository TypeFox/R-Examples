require(devtools)
setwd('~/Development/')

remove.packages('RWebLogo', lib='/Library/Frameworks/R.framework/Versions/Current/Resources/library')

document('RWebLogo/')
system('rm -rf RWebLogo/build/RWebLogo_1.0.3.tar.gz')
system('R CMD BUILD RWebLogo')

system('mv RWebLogo_1.0.3.tar.gz RWebLogo/build/RWebLogo_1.0.3.tar.gz')
system('R CMD INSTALL RWebLogo/build/RWebLogo_1.0.3.tar.gz')

detach("package:RWebLogo", unload=TRUE)
require(RWebLogo)


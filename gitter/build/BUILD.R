setwd('~/Development/')
require(devtools)
remove.packages('gitter')

targz = 'gitter_1.1.1.tar.gz'
system(sprintf('rm -rf gitter/build/%s', targz))

# Regenerate Rwd files
document('gitter/')

R_PATH = '/Library/Frameworks/R.framework/Versions/3.2/Resources/bin/R'
# Build 
system(sprintf('%s CMD build gitter', R_PATH))

# Move into build directory
system(sprintf('mv %s gitter/build/%s', targz, targz))

# Install 
system(sprintf('%s CMD INSTALL gitter/build/%s', R_PATH, targz))

#detach("package:gitter", unload=TRUE)

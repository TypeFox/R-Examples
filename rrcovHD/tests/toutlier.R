## VT::15.09.2013 - this will render the output independent
##  from the version of the package
suppressPackageStartupMessages(library(rrcovHD))

data(hemophilia)
hemophilia$gr <- factor(hemophilia$gr)

obj <- OutlierSign1(gr~., data=hemophilia)
getDistance(obj)            # returns an array of distances
getClassLabels(obj, 1)      # returns an array of indices for a given class
getCutoff(obj)              # returns an array of cutoff values (for each class, usually equal)
getFlag(obj)                # returns an 0/1 array of flags
plot(obj, class=2)          # standard plot function

obj <- OutlierSign2(gr~., data=hemophilia)
getDistance(obj)            # returns an array of distances
getClassLabels(obj, 1)      # returns an array of indices for a given class
getCutoff(obj)              # returns an array of cutoff values (for each class, usually equal)
getFlag(obj)                # returns an 0/1 array of flags
plot(obj, class=2)          # standard plot function

obj <- OutlierPCDist(gr~., data=hemophilia)
getDistance(obj)            # returns an array of distances
getClassLabels(obj, 1)      # returns an array of indices for a given class
getCutoff(obj)              # returns an array of cutoff values (for each class, usually equal)
getFlag(obj)                # returns an 0/1 array of flags
plot(obj, class=2)          # standard plot function

obj <- OutlierPCOut(gr~., data=hemophilia)
getDistance(obj)            # returns an array of distances
getClassLabels(obj, 1)      # returns an array of indices for a given class
getCutoff(obj)              # returns an array of cutoff values (for each class, usually equal)
getFlag(obj)                # returns an 0/1 array of flags
plot(obj, class=2)          # standard plot function

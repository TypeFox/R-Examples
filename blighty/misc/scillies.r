# demonstration r script for generating blighty() compatible
# maps from truncated .fig files - read the comments in the file
# functions.r

# get the functions getdata(), invertfeature() and rescale()
source("functions.r")

# absolutely critical the limits are in the form
# min x, max x, min y, max y
# and correspond to the maxima and minima of the
# featureset you are processing
# for the Isles of Scilly are the extremes measured from an OS map
# and must be done manually
limits <- c(317.62, 467.3, 983.01, 1218.57)


# assign a vector of filenames to be processed
# best got by ls --color=never > files.txt and 
# including that
filenames <- c(
"Scillies-Annet.fig",
"Scillies-Bryher.fig",
"Scillies-Eastern-Islands.fig",
"Scillies-Gugh.fig",
"Scillies-Round-Island.fig",
"Scillies-Samson.fig",
"Scillies-St-Agnes.fig",
"Scillies-St-Martins.fig",
"Scillies-St-Marys.fig",
"Scillies-Tean.fig",
"Scillies-Tresco.fig",
"Scillies-White-Island.fig"
)

# another really important bit is to assign what types the features are going
# to be - so far blighty() recognises:
# 1 as a primary area - ie: Scotland, or the Isle of Man or other landmass
# 2 as a secondary area - ie: a lake
# 3 as a linear feature - ie: a river
# there must be the same number of featuretypes as feature files - the order
# must be the same as in the vector of features
featuretypes <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

# test to see if user has updated the featuretypes vector - exit if not
if(length(featuretypes) != length(filenames)){stop("You've forgotten to update the feature types vector")}

noobjects <- length(filenames)

# automatically generate the objectnames from the filenames
obs <- strsplit(make.names(filenames), split=".fig")
objectnames <- rep(0, noobjects)
for(ctr in 1:noobjects){objectnames[ctr] <- obs[[ctr]]}

# import the data from the filenames and assign to
# the objectnames
for(ctr in 1:noobjects)
	{
	tmp <- getdata(filenames[ctr])
	assign(objectnames[ctr], tmp)
	}

# construct the allpoints vector - allpoints just allows the extreme points
# to be extracted so invertfeature() and rescale() can do their stuff
allpoints <- get(objectnames[1])
for(ctr in 2:noobjects){allpoints <- rbind(allpoints, get(objectnames[ctr]))}

# invert feature converts from the xfig norm of having the zero in the nw corner +ve
# being towards the s - to a conventional Cartesian form
for(ctr in 1:noobjects){assign(objectnames[ctr], invertfeature(get(objectnames[ctr]), allpoints))}

# reconstruct allpoints as all the ys have now changed
allpoints <- get(objectnames[1])
for(ctr in 2:noobjects){allpoints <- rbind(allpoints, get(objectnames[ctr]))}

# rescale every point in the featureset to the minima and maxima set out in limits
for(ctr in 1:noobjects){assign(objectnames[ctr], rescale(get(objectnames[ctr]), allpoints, limits, featuretype=featuretypes[ctr]))}

# reconstruct allpoints as all the ys have now changed
allpoints <- get(objectnames[1])[2:nrow(get(objectnames[1])),]
for(ctr in 2:noobjects){allpoints <- rbind(allpoints, get(objectnames[ctr])[2:nrow(get(objectnames[ctr])),])}

# crudly plot them out
#for(ctr in 1:noobjects){points(get(objectnames[ctr]), type="l")}
#for(ctr in 2:3){points(get(objectnames[ctr]), type="l")}
plot(allpoints, type="l")

for(ctr in 1:noobjects)
	{
	write.table(round(get(objectnames[ctr]), digits=4),
		    file=paste(objectnames[ctr], "txt", sep="."),
	    	    row.names=FALSE,
		    col.names=c("x","y"), 
		    append=FALSE,
		    eol="\n",
		    quote=FALSE)
	}


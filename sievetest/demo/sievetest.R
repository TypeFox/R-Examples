library(sievetest)

# Load demo data provided for testing
data(lignite)

# Load the same data from file
fi <- system.file("lignite.csv",package="sievetest",mustWork=T)
my_std <- read.std(file=fi)

# plot RR diagram, must be the same in this case
plot(lignite)
plot(my_std)

# plot RR diagram for particular sample
plot(lignite[1])

# sumarize all samples
summary(lignite)

# create the object by hand
# in microns:
sieve_aperture_size <- c(500, 200, 90, 0)
# in percents:
mass_ppc_retained <- c(1.01, 24, 42.8, 32.190)
# some metadata:
md <- desc.std(Title="Coal powder, learning std")
# sieve test data (std) object
my_std2 <- std(a=sieve_aperture_size, r=mass_ppc_retained, desc=md)

# create and plot new group of samples
new_gr <- c(my_std2,lignite[2])
plot(new_gr)

# accessing regression summary of particular sample
summary(lignite[[2]]$lmfit)

# new std with adjusted weights for linear regression
newdesc <- desc.std(Title="Lignite (weghted lm)",x=lignite[2])
lignite_twk <- tweak.std(x=lignite[2],lmargs=list(weights=c(0,1,1,1)),desc=newdesc)
new_gr2 <- c(lignite[2],lignite_twk)
plot(new_gr2)



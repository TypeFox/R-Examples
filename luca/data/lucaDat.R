require(genetics)
lucaDat<-read.csv("lucaDat.txt")
lucaDat[,"g"]<-genotype(lucaDat[,"g"])

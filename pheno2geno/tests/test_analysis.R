require(pheno2geno)
#setwd("C:/Users/Konrad/Documents/Github/phenotypes2genotypes/tests")
children  <- read.csv(file="offspring_phenotypes.csv",header=TRUE,row.names=1)
parents   <- read.csv(file="parental_phenotypes.csv",header=TRUE,row.names=1)
genotypes <- read.csv(file="genotypes.csv",header=TRUE,row.names=1)
map       <- read.csv(file="map.csv",header=TRUE,row.names=1)

#with parental data
population <- create.population(children,parents,c(0,0,0,0,0,0,1,1,1,1,1,1),genotypes,mapsPhysical=map,verbose=TRUE)
population <- find.diff.expressed(population)
population <- generate.biomarkers(population, threshold=0.001, margin=5, pProb=0.8, verbose=T, debug=2)
population <- scan.qtls(population,verbose=T,step=4, map="physical", epistasis = "ignore")

####THREE WAYS TO ASSIGN CHROMOSOMES
set.seed(101010)
cross_newmap    <- cross.denovo(population,n.chr=16,map="physical",comparisonMethod=sumMajorityCorrelation,reOrder=TRUE,use.orderMarkers=FALSE,verbose=TRUE,debugMode=2)
cross_saturated <- cross.saturate(population,map="physical",verbose=TRUE,debugMode=2)

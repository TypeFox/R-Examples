## cross-checking sources
# the reference sources
ref.source <- thermo$refs$key
# sources of elemental data
element.source <- thermo$element$source
# sources in the primary thermodynamic database
os1 <- thermo$obigt$ref1
os2 <- thermo$obigt$ref2
# sources also in the supplemental database (OBIGT-2.csv)
add.obigt()
os3 <- thermo$obigt$ref1
os4 <- thermo$obigt$ref2
data(thermo)
# all of the thermodynamic data sources - some of them might be NA
obigt.source <- unique(c(os1,os2,os3,os4))
obigt.source <- obigt.source[!is.na(obigt.source)]
# sources of protein compositions
protein.source <- thermo$protein$ref
# sources of stress response proteins
stressfile <- system.file("extdata/abundance/stress.csv", package="CHNOSZ")
stressdat <- read.csv(stressfile, check.names=FALSE, as.is=TRUE)
stress.source <- as.character(stressdat[2,])
# if the sources are all accounted for 
# these all produce character(0)
print("missing these sources for elemental properties:")
print(unique(element.source[!(element.source %in% ref.source)]))
print("missing these sources (1) for thermodynamic properties:")
print(unique(obigt.source[!(obigt.source %in% ref.source)]))
print("missing these sources for protein compositions:")
print(unique(protein.source[!(protein.source %in% ref.source)]))
print("missing these sources for stress response experiments:")
print(unique(stress.source[!(stress.source %in% ref.source)]))
# determine if all the reference sources are cited
my.source <- c(element.source,obigt.source,protein.source,stress.source)
# this should produce character(0)
print("these sources are present but not cited:")
print(ref.source[!(ref.source %in% my.source)])

require(SensMixed)
load(system.file("testdata","senscherry.RData",package="SensMixed"))

testRep <- FALSE

if(testRep){
#convert some variables to factors in TVbo
TVbo <- convertToFactors(TVbo, c("Assessor", "Repeat", "Picture"))
                         
res_paral <- sensmixed(c("Coloursaturation", "Colourbalance"),
                       Prod_effects = c("TVset", "Picture"), replication="Repeat", 
                       individual="Assessor", data=TVbo, error_structure="3-WAY")

stopifnot(nrow(res_paral$random$Chi)==11)

senscherry$Replicate <- as.factor(senscherry$Replicate)
res_cherry<- sensmixed(c("Apple.aroma", "Cherry.aroma"),
                       Prod_effects = c("Fla", "Fib"), replication="Replicate", 
                       individual="Assessor", data = senscherry, 
                       error_structure="3-WAY")

stopifnot(nrow(res_cherry$random$Chi)==11)

res_cherryMAMmult<- sensmixed(c("Apple.aroma", "Cherry.aroma"),
                       Prod_effects = c("Fla", "Fib"), replication="Replicate", 
                       individual="Assessor", data = senscherry, 
                       error_structure="3-WAY", MAM=TRUE)

res_cherryMAMmult$scaling

## with one sample
res_MAM <- sensmixed(c("Apple.aroma", "Cherry.aroma"),
                     Prod_effects=c("Sample"), replication="Replicate", 
                     individual="Assessor", data=senscherry, MAM=TRUE, 
                     reduce.random = FALSE)

res_MAM_PER <- sensmixed(c("Apple.aroma", "Cherry.aroma"),
                     Prod_effects=c("Sample"), replication="Replicate", 
                     individual="Assessor", data=senscherry, MAM_PER=TRUE)


stopifnot(all.equal(res_MAM$scaling$FScaling[,"Apple.aroma"],
                    res_MAM_PER[[3]][,,1]["Scaling","F"], tol=1e-2))

stopifnot(all.equal(res_MAM$scaling$FScaling[,"Cherry.aroma"],
                    res_MAM_PER[[3]][,,2]["Scaling","F"], tol=1e-2))

}


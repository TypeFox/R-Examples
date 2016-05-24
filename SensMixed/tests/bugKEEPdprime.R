require(SensMixed)

## here is an example when all the random effects are eliminated
## TODO: make an appropriate output for the print method
tools::assertError(print(res2 <- sensmixed(c("Noise", "Elasticeffect"),
                  Prod_effects = c("TVset"), replication="Repeat", 
                  individual="Assessor", data=TVbo, parallel = FALSE)))
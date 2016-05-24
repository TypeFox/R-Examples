require(SensMixed)

#convert some variables to factors in TVbo
TVbo <- convertToFactors(TVbo, c("Assessor", "Repeat", "Picture"))

res <- sensmixed(c("Coloursaturation", "Colourbalance"),
                 Prod_effects = c("TVset", "Picture"), replication="Repeat", 
                 individual="Assessor", data=TVbo, keep.effs = "Assessor")

stopifnot(res$step_res[[1]]$rand.table["Assessor","elim.num"] == "kept")
stopifnot(res$step_res[[2]]$rand.table["Assessor","elim.num"] == "kept")

res <- sensmixed(c("Coloursaturation", "Colourbalance"),
                 Prod_effects = c("TVset", "Picture"), replication="Repeat", 
                 individual="Assessor", data=TVbo, 
                 keep.effs = c("Assessor", "TVset:Picture:Assessor"))

stopifnot(res$step_res[[1]]$rand.table["Assessor","elim.num"] == "kept")
stopifnot(res$step_res[[2]]$rand.table["Assessor","elim.num"] == "kept")
stopifnot(res$step_res[[1]]$rand.table["TVset:Picture:Assessor","elim.num"] == "kept")
stopifnot(res$step_res[[2]]$rand.table["TVset:Picture:Assessor","elim.num"] == "kept")

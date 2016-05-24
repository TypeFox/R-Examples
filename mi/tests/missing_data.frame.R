stopifnot(require(mi))

rdf <- rdata.frame(N = 100, n_partial = 2, n_full = 2)
mdf <- missing_data.frame(rdf$obs)

rdf <- rdata.frame(N = 100, n_partial = 6, n_full = 1, 
                   types = c("ordinal", "cont", "count", "binary",
                             "proportion", "positive", "nominal"))
mdf <- missing_data.frame(rdf$obs)
mdf <- missing_data.frame(rdf$obs, favor_positive = TRUE)

rdf <- rdata.frame(N = 100, n_partial = 5, n_full = 1, experiment = TRUE,
                   types = c("treatment", "cont", "count", "binary",
                             "proportion", "positive"))
mdf <- missing_data.frame(rdf$obs, subclass = "experiment", concept = 
                          as.factor(c("treatment", rep("covariate", 4), "outcome")))

rdf <- rdata.frame(N = 100, n_partial = 5, n_full = 0, types = "ordinal")
mdf <- missing_data.frame(rdf$obs, subclass = "allcategorical")

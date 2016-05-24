#tests for cjoint main function

data("immigrationconjoint")
data("immigrationdesign")

results <- amce(Chosen_Immigrant ~  Gender + Education + `Language Skills` + `Country of Origin` + Job + `Job Experience` + `Job Plans` + `Reason for Application` + `Prior Entry`, data=immigrationconjoint, cluster=TRUE, respondent.id="CaseID", design=immigrationdesign)

test_that("Output class",{
    expect_is(results,"amce")
})


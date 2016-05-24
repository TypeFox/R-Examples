context("Score matrix creation functions")

testneurons <- readRDS('testdata/testneurons.rds')

test_that("calculation of probability matrices is correct", {
  prob_mat.expected <- structure(c(0.004617139042639, 0.012043501237727, 0.00937334853836954, 0.00859455400105693, 0.0804939782493811, 0.102967763468973, 0.0157984034711985, 0.00408867132089116, 0.0137123466748255, 0.00967930353517092, 0.00812171446054571, 0.0818568686896782, 0.103969070731232, 0.0157149611993436, 0.00511779267376853, 0.0154368202931605, 0.0103468417100103, 0.00906739354156816, 0.0884209940755987, 0.109949100214168, 0.0137401607654438, 0.00603565766417267, 0.01593747392429, 0.011320334881651, 0.0115706616972158, 0.102077712569188, 0.120379384196034, 0.00956804717269769), class = "table", distbreaks = c(1, 2, 5, 10, 20, 50, 100, 500), dotprodbreaks = c(0, 0.25, 0.5, 0.75, 1), .Dim = c(7L, 4L), .Dimnames = structure(list(c("(1,2]", "(2,5]", "(5,10]", "(10,20]", "(20,50]", "(50,100]", "(100,500]" ), c("(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")), .Names = c("", "")))
  prob_mat.calculated <- calc_prob_mat(calc_dists_dotprods(testneurons[1:3], testneurons[2:5]), distbreaks=c(1, 2, 5, 10, 20, 50, 100, 500), dotprodbreaks=seq(0, 1, 0.25))
  expect_equal(prob_mat.expected, prob_mat.calculated)
})

test_that("calculation of score matrices is correct", {
  prob_mat_1 <- calc_prob_mat(calc_dists_dotprods(testneurons[1:3], testneurons[2:5]), distbreaks=c(1, 2, 5, 10, 20, 50, 100, 500), dotprodbreaks=seq(0, 1, 0.25))
  prob_mat_2 <- calc_prob_mat(calc_dists_dotprods(testneurons[2:4], testneurons[1:4]), distbreaks=c(1, 2, 5, 10, 20, 50, 100, 500), dotprodbreaks=seq(0, 1, 0.25))
  score_mat.expected <- structure(c(5.58670412933517, 3.54392279704006, 1.44948841421627, -0.428948458404255, -0.418965476490948, 0.274653663203646, 6.47383688918051, 2.76488034723549, 3.14571143338443, 1.29716493780848, -0.41784692720237, -0.422732257032202, 0.280304729738837, 6.47001734289776, 5.66093483703527, 2.73528135640762, 1.14285516828316, -0.386306504542588, -0.426345281997543, 0.257917833825592, 6.3731606227297, 5.77988662143569, 2.30679302110432, 0.920372423679295, -0.243946421089446, -0.426020303089466, 0.331047071615182, 6.11215516036325), class = c("scoringmatrix","table"), distbreaks = c(1, 2, 5, 10, 20, 50, 100, 500), dotprodbreaks = c(0, 0.25, 0.5, 0.75, 1), .Dim = c(7L, 4L), .Dimnames = structure(list(c("(1,2]", "(2,5]", "(5,10]", "(10,20]", "(20,50]", "(50,100]", "(100,500]"), c("(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")), .Names = c("", "")))
  score_mat.actual <- calc_score_matrix(prob_mat_1, prob_mat_2, logbase=4, epsilon=2e-6)
  expect_equal(score_mat.expected, score_mat.actual)
})

test_that("neuron_pairs",{
  expect_equal(allpairs<-neuron_pairs(kcs20, n=NA, ignoreSelf=FALSE),
               expand.grid(query=names(kcs20), target=names(kcs20), stringsAsFactors = FALSE, KEEP.OUT.ATTRS=FALSE))
  expect_equal(neuron_pairs(kcs20, n=NA), allpairs[allpairs$target!=allpairs$query,])
  expect_is(df<-neuron_pairs(kcs20, n=10), 'data.frame')
  expect_equal(nrow(df), 10L)
  expect_is(df$query, 'character')
})

test_that("calc_dists_dotprods",{
  # check that we can show progress bar
  expect_output(dd<-calc_dists_dotprods(kcs20[1:5],.progress='text'), '==')
  # list of length 20, with each element 2 long
  expect_is(dd, 'list')
  baseline_dims <- structure(c(297L, 2L, 374L, 2L, 355L, 2L, 369L, 2L, 284L, 2L,
    374L, 2L, 355L, 2L, 369L, 2L, 284L, 2L, 297L, 2L, 355L, 2L, 369L,
    2L, 284L, 2L, 297L, 2L, 374L, 2L, 369L, 2L, 284L, 2L, 297L, 2L,
    374L, 2L, 355L, 2L), .Dim = c(2L, 20L), .Dimnames = list(NULL, 1:20))
  # calc_dists_dotprods produces the right format output
  expect_equal(sapply(dd, dim), baseline_dims)
  expect_is(dd[[1]]$nndists, 'numeric')
})

test_that("create_scoringmatrix wrapper function passes on arguments", {
  smat.expected <- structure(c(0, -8.3182108909519, -11.1402936474687, -0.420254621476833, -0.181960193657874, 0.282642512418894, 0.675549776384255, -6.45574517467861, -8.66516528624987, -11.5271663713186, -0.580326857349317, -0.173610788907039, 0.312655458083253, 0.675566810158546, 0, -9.44129133098147, -11.8319278406857, -0.819298784485771, -0.164647753087521, 0.305617207199979, 0.675569104424631, 0, -9.98495948997948, -12.3859031633819, -1.26000963427592, -0.243916155861412, 0.344935126998013, 0.675538251178387), class = c("scoringmatrix", "table"), distbreaks = c(1, 2, 5, 10, 20, 50, 100, 500), dotprodbreaks = c(0, 0.25, 0.5, 0.75, 1), .Dim = c(7L, 4L), .Dimnames = structure(list(c("(1,2]", "(2,5]", "(5,10]", "(10,20]", "(20,50]", "(50,100]", "(100,500]" ), c("(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")), .Names = c("", "")))
  smat.calculated <- create_scoringmatrix(testneurons[1:3], testneurons[1:4],
                                 distbreaks=c(1, 2, 5, 10, 20, 50, 100, 500),
                                 dotprodbreaks=seq(0, 1, 0.25),
                                 non_matching_subset=neuron_pairs(testneurons[1:4]))
  expect_equal(smat.expected, smat.calculated)
})

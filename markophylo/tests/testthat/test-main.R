context("Model fit output testing")

#right value, right class
test_that("Data 1", {
  data(simdata1)
  changetodf <- as.data.frame(simdata1$data)
  model1 <- markophylo::estimaterates(usertree = simdata1$tree, userphyl = changetodf, alphabet = c(1, 2), rootprob = "equal", modelmat = "ARD")  
  expect_equal(model1$results$wop$par, c(0.568075375715259, 1.0111462398653), tolerance = 0.1, scale = 1 )
  expect_equal(model1$results$wop$se, c(0.0231519067611956, 0.0311917021383594), tolerance = 0.01, scale = 1 )
  expect_equal(model1$bg, list(1:7), tolerance = 0.01, scale = 1 )
  expect_equal(model1$results$wop$objective, 10055.01597621022, tolerance = 0.1, scale = 1 )
  expect_equal(model1$results$wop$df, 2, tolerance = 0.01, scale = 1 )
  expect_equal(model1$results$wop$AIC, -20114.032, tolerance = 0.1, scale = 1 )
  expect_equal(model1$results$wop$BIC, -20127.0663863828, tolerance = 0.1, scale = 1 )
  expect_equal(model1$phyl, 5000, tolerance = 0.1, scale = 1 )
})

test_that("Data 1 - GTR", {
  data(simdata1)
  model1_rev <- estimaterates(usertree = simdata1$tree, userphyl = simdata1$data, 
                              alphabet = c(1, 2), rootprob = "maxlik", 
                              modelmat = "GTR") 
  expect_equal(model1_rev$results$wop$par, c(1.54841557290695, -0.247752943473407), tolerance = 0.1, scale = 1 )
  expect_equal(model1_rev$results$wop$se, c(0.0375002451746298, 0.00553408019028168), tolerance = 0.01, scale = 1 )
  expect_equal(model1_rev$results$wop$revcombined, structure(list(rates = structure(c(0.678789236864455, 0.869626336042494
  ), .Dim = c(2L, 1L, 1L))), .Names = "rates"), tolerance = 0.01, scale = 1 )
  expect_equal(model1_rev$bg, list(1:7), tolerance = 0.01, scale = 1 )
  expect_equal(model1_rev$results$wop$objective, 10057.6252820299, tolerance = 0.1, scale = 1 )
  expect_equal(model1_rev$results$wop$df, 2, tolerance = 0.01, scale = 1 )
  expect_equal(model1_rev$results$wop$AIC, -20119.25, tolerance = 0.1, scale = 1 )
  expect_equal(model1_rev$results$wop$BIC, -20132.2843863828, tolerance = 0.1, scale = 1 )
})

test_that("Data 1 - unobserved data correction", {
  data(simdata1)
  filterall1 <- which(apply(simdata1$data, MARGIN = 1, FUN = 
                              function(x) isTRUE(all.equal(as.vector(x), c(1, 1, 1, 1)))))
  filterall0 <- which(apply(simdata1$data, MARGIN = 1, FUN = 
                              function(x) isTRUE(all.equal(as.vector(x), c(0, 0, 0, 0)))))
  filteredsimdata1 <- simdata1$data[-c(filterall1, filterall0), ]
  model1_f_corrected <- estimaterates(usertree = simdata1$tree, userphyl = filteredsimdata1, unobserved = matrix(c(0, 0, 0, 0, 1, 1, 1, 1), nrow = 2,   byrow = TRUE), alphabet = c(1, 2), rootprob = "equal", modelmat = "ARD")
  expect_equal(model1_f_corrected$results$wop$par, c(0.575499276635208, 0.960368220064815), tolerance = 0.1, scale = 1 )
  expect_equal(model1_f_corrected$results$wop$se, c(0.0273406890152568, 0.0992307442203887), tolerance = 0.01, scale = 1 )
  expect_equal(model1_f_corrected$bg, list(1:7), tolerance = 0.01, scale = 1 )
  expect_equal(model1_f_corrected$results$wop$objective, 7363.33815280957, tolerance = 0.1, scale = 1 )
  expect_equal(model1_f_corrected$results$wop$df, 2, tolerance = 0.01, scale = 1 )
  expect_equal(model1_f_corrected$results$wop$AIC, -14730.676, tolerance = 0.1, scale = 1 )
  expect_equal(model1_f_corrected$results$wop$BIC, -14743.1897336979, tolerance = 0.1, scale = 1 )
})

test_that("Data 2", {
  data(simdata2)
  model2 <- estimaterates(usertree = simdata2$tree, userphyl = simdata2$data, alphabet = c(1, 2),
                          bgtype = "ancestornodes", bg = c(7),
                          rootprob = "equal", modelmat = matrix(c(NA, 1, 2, NA), 2, 2)) 
  expect_equal(model2$results$wop$par, c(1.12572219285505, 2.00102312502631, 0.558698074013352, 1.00619184806976
  ), tolerance = 0.1, scale = 1 )
  expect_equal(model2$results$wop$se, c(0.0720311958948751, 0.105554886394708, 0.0267706299283269, 
                                        0.0352250699131224), tolerance = 0.01, scale = 1 )
  expect_equal(model2$bg, list(c(7, 3, 4), c(1, 2, 5, 6, 7)), tolerance = 0.01, scale = 1 )
  expect_equal(model2$results$wop$objective, 10893.4834080664, tolerance = 0.1, scale = 1 )
  expect_equal(model2$results$wop$df, 4, tolerance = 0.01, scale = 1 )
  expect_equal(model2$results$wop$AIC, -21794.966, tolerance = 0.1, scale = 1 )
  expect_equal(model2$results$wop$BIC, -21821.0347727657, tolerance = 0.1, scale = 1 )
})

test_that("Data 3", {
  data(simdata3)
  model3 <- estimaterates(usertree = simdata3$tree, userphyl = simdata3$data, 
                          alphabet = c("a", "c", "g", "t"), rootprob = "equal", 
                          partition = list(c(1:2500), c(2501:5000)), 
                          modelmat = "ER")
  expect_equal(model3$results$wop$par, c(0.34515494556094, 0.986983659658561), tolerance = 0.1, scale = 1 )
  expect_equal(model3$results$wop$se, c(0.0091547389802298, 0.0240169425220838), tolerance = 0.01, scale = 1 )
  expect_equal(model3$bg, list(1:7), tolerance = 0.01, scale = 1 )
  expect_equal(model3$results$wop$objective, 21659.3772016081, tolerance = 0.1, scale = 1 )
  expect_equal(model3$results$wop$df, 2, tolerance = 0.01, scale = 1 )
  expect_equal(model3$results$wop$AIC, -43322.754, tolerance = 0.1, scale = 1 )
  expect_equal(model3$results$wop$BIC, -43335.7883863828, tolerance = 0.1, scale = 1 )
})

test_that("Data 3", {
  data(simdata3)
  model3 <- estimaterates(usertree = simdata3$tree, userphyl = simdata3$data, 
                          alphabet = c("a", "c", "g", "t"), rootprob = "equal", 
                          partition = list(c(1:2500), c(2501:5000)), 
                          modelmat = "ER")
  expect_equal(model3$results$wop$par, c(0.34515494556094, 0.986983659658561), tolerance = 0.1, scale = 1 )
  expect_equal(model3$results$wop$se, c(0.0091547389802298, 0.0240169425220838), tolerance = 0.01, scale = 1 )
  expect_equal(model3$bg, list(1:7), tolerance = 0.01, scale = 1 )
  expect_equal(model3$results$wop$objective, 21659.3772016081, tolerance = 0.1, scale = 1 )
  expect_equal(model3$results$wop$df, 2, tolerance = 0.01, scale = 1 )
  expect_equal(model3$results$wop$AIC, -43322.754, tolerance = 0.1, scale = 1 )
  expect_equal(model3$results$wop$BIC, -43335.7883863828, tolerance = 0.1, scale = 1 )
})

test_that("Data 3; stationary; ancestornodes", {
  skip_on_cran()
  data(simdata3)
  model3_st <- estimaterates(usertree = simdata3$tree, userphyl = simdata3$data, 
                             bgtype="anc",bg=7,                      alphabet = c("a", "c", "g", "t"), rootprob = "stationary", 
                             partition = list(c(1:2500), c(2501:5000)), 
                             modelmat = "ER")
  expect_equal(model3_st$results$wop$par, c(0.343850888099329, 0.345516097260013, 1.03215162798544, 0.957475743441639
  ), tolerance = 0.1, scale = 1 )
  expect_equal(model3_st$results$wop$se, c(0.0199493856186924, 0.0103990348418819, 0.0404358593196127, 
                                           0.0305548719185289), tolerance = 0.01, scale = 1 )
  expect_equal(model3_st$bg, list(c(7, 3, 4), c(1, 2, 5, 6, 7)), tolerance = 0.01, scale = 1 )
  expect_equal(model3_st$results$wop$objective, 21658.3371702953, tolerance = 0.1, scale = 1 )
  expect_equal(model3_st$results$wop$df, 4, tolerance = 0.01, scale = 1 )
  expect_equal(model3_st$results$wop$AIC, -43324.674, tolerance = 0.1, scale = 1 )
  expect_equal(model3_st$results$wop$BIC, -43350.7427727657, tolerance = 0.1, scale = 1 )
})

test_that("Data 3 - unobserved data correction", {
  data(simdata3)
  a3 <- which(apply(simdata3$data, MARGIN = 1, FUN = function(x) all(x==rep("a", 4))) ) 
  newpart <- list(c(1:(2500 - 268)), c((2500 - 268):(5000 - 367)))
  model3_f <- estimaterates(usertree = simdata3$tree, userphyl = simdata3$data[-a3, ], unobserved = matrix(rep("a", 4), nrow = 1), alphabet = c("a", "c", "g", "t"), rootprob = "maxlik", partition = newpart, modelmat = matrix(c(NA, 1, 1, 1, 1, NA, 1, 1, 1, 1, NA, 1, 1, 1, 1, NA), 4, 4))
  expect_equal(model3_f$results$wop$par, c(0.336100048839514, 0.992773997420928, 0.179678889898291, 0.0746902904518238, 0.0563024363215518), tolerance = 0.1, scale = 1 )
  expect_equal(model3_f$results$wop$se, c(0.0101884583407262, 0.0260178790640918, 0.0135012818511467, 0.00971320301974167, 0.00962553467325101), tolerance = 0.01, scale = 1 )
  expect_equal(model3_f$bg, list(1:7), tolerance = 0.01, scale = 1 )
  expect_equal(model3_f$results$wop$objective, 20390.6316016373, tolerance = 0.1, scale = 1 )
  expect_equal(model3_f$results$wop$df, 5, tolerance = 0.01, scale = 1 )
  expect_equal(model3_f$results$wop$AIC, -40791.264, tolerance = 0.1, scale = 1 )
  expect_equal(model3_f$results$wop$BIC, -40823.4687994271, tolerance = 0.1, scale = 1 )
})

test_that("Data 3 - unobserved data correction - discgamma", {
  skip_on_cran()
  data(simdata3)
  a3 <- which(apply(simdata3$data, MARGIN = 1, FUN = function(x) all(x==rep("a", 4))) ) 
  newpart <- list(c(1:(2500 - 268)), c((2500 - 268):(5000 - 367)))
  model3_f_dg <- estimaterates(usertree = simdata3$tree, userphyl = simdata3$data[-a3, ], unobserved = matrix(rep("a", 4), nrow = 1), alphabet = c("a", "c", "g", "t"), rootprob = "maxlik", ratevar = "discgamma", nocat = 4, partition = newpart, modelmat = matrix(c(NA, 1, 1, 1, 1, NA, 1, 1, 1, 1, NA, 1, 1, 1, 1, NA), 4, 4))
  expect_equal(model3_f_dg$results$wop$par, c(0.338255127691846, 1.01390829271935, 17.7979992412363, 0.197005948985821, 
                                              0.0751489903938522, 0.0574232192702603), tolerance = 0.1, scale = 1 )
  expect_equal(model3_f_dg$results$wop$se, c(0.0106683363012078, 0.0354250729563796, 18.8024988413938, 0.0141464905128635, 
                                             0.00981422798235644, 0.00970538694203939), tolerance = 0.01, scale = 1 )
  expect_equal(model3_f_dg$bg, list(1:7), tolerance = 0.01, scale = 1 )
  expect_equal(model3_f_dg$results$wop$objective, 20390.1790259543, tolerance = 0.1, scale = 1 )
  expect_equal(model3_f_dg$results$wop$df, 6, tolerance = 0.01, scale = 1 )
  expect_equal(model3_f_dg$results$wop$AIC, -40792.358, tolerance = 0.1, scale = 1 )
  expect_equal(model3_f_dg$results$wop$BIC, -40831.0037593125, tolerance = 0.1, scale = 1 )
})

test_that("Data 3 - unobserved data correction - partionspecificgamma", {
  skip_on_cran()
  data(simdata3)
  a3 <- which(apply(simdata3$data, MARGIN = 1, FUN = function(x) all(x==rep("a", 4))) ) 
  newpart <- list(c(1:(2500 - 268)), c((2500 - 268):(5000 - 367)))
  model3_f_pg <- estimaterates(usertree = simdata3$tree, userphyl = simdata3$data[-a3, ], unobserved = matrix(rep("a", 4), nrow = 1), alphabet = c("a", "c", "g", "t"), rootprob = "maxlik", ratevar = "partitionspecificgamma", nocat = 4, partition = newpart, modelmat = matrix(c(NA, 1, 1, 1, 1, NA, 1, 1, 1, 1, NA, 1, 1, 1, 1, NA), 4, 4))
  expect_equal(model3_f_pg$results$wop$par, c(0.33604061006274, 1.02265413630596, 79.7497382822248, 12.8415242617211, 
                                              0.194645772368142, 0.0753286688219997, 0.0576588947412515), tolerance = 0.1, scale = 1 )
  expect_equal(model3_f_pg$results$wop$se[-3], c(0.0112442463243987, 0.0395626597843926, 633.502298790985, 11.8347656260673, 
                                             0.0141500287278521, 0.00981353464923172, 0.0097061504921285)[-3], tolerance = 0.01, scale = 1 )
  expect_equal(model3_f_pg$bg, list(1:7), tolerance = 0.01, scale = 1 )
  expect_equal(model3_f_pg$results$wop$objective, 20390.0331191522, tolerance = 0.1, scale = 1 )
  expect_equal(model3_f_pg$results$wop$df, 7, tolerance = 0.01, scale = 1 )
  expect_equal(model3_f_pg$results$wop$AIC, -40794.066, tolerance = 0.1, scale = 1 )
  expect_equal(model3_f_pg$results$wop$BIC, -40839.1527191979, tolerance = 0.1, scale = 1 )
})

test_that("Data 4", {
  skip_on_cran()
  data(simdata4)
  model4 <- estimaterates(usertree = simdata4$tree, userphyl = simdata4$data, 
                          alphabet = c("a", "c", "g", "t"), rootprob = "maxlik", 
                          ratevar = "discgamma", nocat = 4, 
                          modelmat = matrix(c(NA, 1, 1, 1, 1, NA, 1, 1, 
                                              1, 1, NA, 1, 1, 1, 1, NA), 4, 4))
  expect_equal(model4$results$wop$par, c(0.332227137226387, 1.99516075159322, -0.766684671257096, 0.452226007919348, 
                                         -0.00262980261450702), tolerance = 0.1, scale = 1 )
  expect_equal(model4$results$wop$se, c(0.00852916557446111, 0.313907966035631, 0.00550422638993651, 
                                        0.00827188807406831, 0.00739162398022605), tolerance = 0.01, scale = 1 )
  expect_equal(model4$bg, list(1:7), tolerance = 0.01, scale = 1 )
  expect_equal(model4$results$wop$objective, 17763.3255939325, tolerance = 0.1, scale = 1 )
  expect_equal(model4$results$wop$df, 5, tolerance = 0.01, scale = 1 )
  expect_equal(model4$results$wop$AIC, -35536.652, tolerance = 0.1, scale = 1 )
  expect_equal(model4$results$wop$BIC, -35569.2379659571, tolerance = 0.1, scale = 1 )
})

test_that("Data 4 - unobservable data correction", {
  skip_on_cran()
  data(simdata4)
  filteralla <- which(apply(simdata4$data, MARGIN = 1, FUN = 
                              function(x) isTRUE(all.equal(as.vector(x), c("a", "a", "a", "a")))))
  filterg3 <- which(apply(simdata4$data, MARGIN = 1, FUN = 
                            function(x) table(x)["g"] >= 3) )
  filteredsimdata4 <- simdata4$data[-c(filteralla, filterg3), ]
  
  alphabet <- c("a", "c", "g", "t")
  allpatt <- expand.grid(alphabet, alphabet, alphabet, alphabet) #all possible combinations
  unob_patt_index_1 <- which(apply(allpatt, MARGIN = 1, FUN = function(x) table(x)["g"] >= 3) )
  unob_patt_index_2 <- which(apply(allpatt, MARGIN = 1, FUN = function(x) isTRUE(all.equal(as.vector(x), c("a", "a", "a", "a"))) ) )
  unob_patt_indices <- sort(union(unob_patt_index_1, unob_patt_index_2)) #Vector of ordered indices of unique patterns
  unob_patt <- allpatt[unob_patt_indices, ] #matrix of unique patterns
  
  model4_f <- estimaterates(usertree = simdata4$tree, userphyl = filteredsimdata4, unobserved = unob_patt, alphabet = c("a", "c", "g", "t"), rootprob = "maxlik", 
                             ratevar = "discgamma", nocat = 4, 
                             modelmat = matrix(c(NA, 1, 1, 1, 1, NA, 1, 1, 
                                                 1, 1, NA, 1, 1, 1, 1, NA), 4, 4))  
  expect_equal(model4_f$results$wop$par, c(0.343385130827976, 1.98007213707476, -1.0284043828689, 0.451318540486224, 
                                           -0.0826744351623355), tolerance = 0.1, scale = 1 )
  expect_equal(model4_f$results$wop$se, c(0.0115998199951891, 0.340854164988812, 0.0102994299943886, 
                                          0.0167143642272215, 0.0261293184969688), tolerance = 0.01, scale = 1 )
  expect_equal(model4_f$bg, list(1:7), tolerance = 0.01, scale = 1 )
  expect_equal(model4_f$results$wop$objective, 12267.8003986961, tolerance = 0.1, scale = 1 )
  expect_equal(model4_f$results$wop$df, 5, tolerance = 0.01, scale = 1 )
  expect_equal(model4_f$results$wop$AIC, -24545.6, tolerance = 0.1, scale = 1 )
  expect_equal(model4_f$results$wop$BIC, -24576.5267511159, tolerance = 0.1, scale = 1 )
})

test_that("Data 5", {
  skip_on_cran()
  data(simdata5)
  model5 <- estimaterates(usertree = simdata5$tree, userphyl = simdata5$data, 
                          alphabet = c("a", "c", "g", "t"), rootprob = "maxlik", 
                          partition = list((1:6000)[-seq(3, 6000, by = 3)], 
                                           seq(3, 6000, by = 3) ),
                          ratevar = "partitionspecificgamma", nocat = 4, 
                          modelmat = matrix(c(NA, 1, 1, 1, 1, NA, 1, 1, 
                                              1, 1, NA, 1, 1, 1, 1, NA), 4, 4))
  expect_equal(model5$results$wop$par, c(0.322249198289113, 0.328960343834729, 1.62695446447176, 0.576449199632739, 
                                         -0.672877490443373, 0.471148787388902, 0.0684712857046065), tolerance = 0.1, scale = 1 )
  expect_equal(model5$results$wop$se, c(0.00945113136016822, 0.0166998049159744, 0.247562986205858, 
                                        0.0685778501029639, 0.00509606529227118, 0.00747241692587547, 
                                        0.00674804566758415), tolerance = 0.01, scale = 1 )
  expect_equal(model5$bg, list(1:7), tolerance = 0.01, scale = 1 )
  expect_equal(model5$results$wop$objective, 20718.1881784304, tolerance = 0.1, scale = 1 )
  expect_equal(model5$results$wop$df, 7, tolerance = 0.01, scale = 1 )
  expect_equal(model5$results$wop$AIC, -41450.376, tolerance = 0.1, scale = 1 )
  expect_equal(model5$results$wop$BIC, -41497.2726032375, tolerance = 0.1, scale = 1 )
})

test_that("Data 5 - reversible; partition; gamma rate variation; ancestornodes", {
  skip_on_cran()
  data(simdata5)
  model5_rev <- estimaterates(usertree = simdata5$tree, userphyl = simdata5$data, bgtype = "ancestor", bg = 7,
                          alphabet = c("a", "c", "g", "t"), rootprob = "maxlik", 
                          partition = list((1:6000)[-seq(3, 6000, by = 3)], 
                                           seq(3, 6000, by = 3) ),
                          ratevar = "partitionspecificgamma", nocat = 4, 
                          modelmat = "ER", reversible = TRUE)
  expect_equal(model5_rev$results$wop$par, c(1.30488722743129, 1.36612129968556, 1.37153293367857, 1.48780409387041, 
                                              1.51094469557796, 0.498087383269959, -0.467819376271539, 0.343802869950167, 
                                              0.0501228571351223), tolerance = 0.1, scale = 1 )
  expect_equal(model5_rev$results$wop$se, 
               c(0.0657580021701737, 0.0490980293048962, 0.104781719049689, 
                 0.103531658607197, 0.228501707495941, 0.0603850891554088, 0.00360005765087641, 
                 0.00490800779952262, 0.00445813758717441), tolerance = 0.01, scale = 1 )
  expect_equal(model5_rev$results$wop$revcombined, structure(list(rates = structure(c(0.199932569660071, 0.199932569660071, 
                                                                                      0.199932569660071, 0.450159697986254, 0.450159697986254, 0.450159697986254, 
                                                                                      0.33560081108362, 0.33560081108362, 0.33560081108362, 0.31919414870135, 
                                                                                      0.31919414870135, 0.31919414870135, 0.209314748563489, 0.209314748563489, 
                                                                                      0.209314748563489, 0.471284214260898, 0.471284214260898, 0.471284214260898, 
                                                                                      0.351349455014282, 0.351349455014282, 0.351349455014282, 0.33417288184689, 
                                                                                      0.33417288184689, 0.33417288184689, 0.210143909787185, 0.210143909787185, 
                                                                                      0.210143909787185, 0.473151118521048, 0.473151118521048, 0.473151118521048, 
                                                                                      0.352741260159711, 0.352741260159711, 0.352741260159711, 0.335496645210628, 
                                                                                      0.335496645210628, 0.335496645210628, 0.227958776348698, 0.227958776348698, 
                                                                                      0.227958776348698, 0.513262316834716, 0.513262316834716, 0.513262316834716, 
                                                                                      0.382644760512632, 0.382644760512632, 0.382644760512632, 0.363938240174363, 
                                                                         0.363938240174363, 0.363938240174363), .Dim = c(12L, 2L, 2L))), .Names = "rates"), tolerance = 0.01, scale = 1 )
  expect_equal(model5_rev$bg, list(c(7, 3, 4), c(1, 2, 5, 6, 7)), tolerance = 0.01, scale = 1 )
  expect_equal(model5_rev$results$wop$objective, 20650.9142963678, tolerance = 0.1, scale = 1 )
  expect_equal(model5_rev$results$wop$df, 9, tolerance = 0.01, scale = 1 )
  expect_equal(model5_rev$results$wop$AIC, -41319.828, tolerance = 0.1, scale = 1 )
  expect_equal(model5_rev$results$wop$BIC, -41380.1236327339, tolerance = 0.1, scale = 1 )
})

data("hospitals", package="rDEA")

## choosing inputs and outputs for analysis
X = hospitals[,c('labor', 'capital')]
Y = hospitals[,c('inpatients', 'outpatients')]
W = hospitals[,c('labor_price', 'capital_price')]
firms = 1:10

context("Input DEA")

test_that("input DEA with variable RTS works", {
  dea = dea( XREF=X, YREF=Y, X=X[firms,], Y=Y[firms,], model="input", RTS="variable" )
  expect_equal( length(dea$thetaOpt), length(firms) )
  correct_thetaOpt =  c(0.8832473, 0.7915816, 0.9444157, 0.8648906, 0.7307074, 0.7657450, 0.6904628, 0.6453444, 0.8474180, 0.7395350)
  expect_equal( dea$thetaOpt, correct_thetaOpt, tolerance=1e-5 )
  expect_equal( dea$lambda_sum, rep.int(1, length(firms)) )
})

test_that("input DEA with constant RTS works", {
  dea = dea( XREF=X, YREF=Y, X=X[firms,], Y=Y[firms,], model="input", RTS="constant" )
  expect_equal( length(dea$thetaOpt), length(firms) )
  correct_thetaOpt   = c(0.8762136, 0.7860619, 0.9400625, 0.7466703, 0.2705233, 0.7560322, 0.6860306, 0.6304803, 0.8434420, 0.5824253)
  correct_lambda_sum = c(2.3656910, 2.6184443, 1.7116862, 0.1823766, 0.2130936, 4.5457910, 1.8709218, 0.2056953, 2.6507813, 0.4621705)
  expect_equal( dea$thetaOpt, correct_thetaOpt, tolerance=1e-5 )
  expect_equal( dea$lambda_sum, correct_lambda_sum, tolerance=1e-5 )
})

test_that("input DEA with non-increasing RTS works", {
  dea = dea( XREF=X, YREF=Y, X=X[firms,], Y=Y[firms,], model="input", RTS="non-increasing" )
  expect_equal( length(dea$thetaOpt), length(firms) )
  correct_thetaOpt   = c(0.8832473, 0.7915816, 0.9444157, 0.7466703, 0.2705233, 0.7657450, 0.6904628, 0.6304803, 0.8474180, 0.5824253)
  correct_lambda_sum = c(1.0000000, 1.0000000, 1.0000000, 0.1823766, 0.2130936, 1.0000000, 1.0000000, 0.2056953, 1.0000000, 0.4621705)
  expect_equal( dea$thetaOpt, correct_thetaOpt, tolerance=1e-5 )
  expect_equal( dea$lambda_sum, correct_lambda_sum, tolerance=1e-5 )
})

##################################################################

context("Output DEA")

test_that("output DEA with variable RTS works", {
  dea = dea( XREF=X, YREF=Y, X=X[firms,], Y=Y[firms,], model="output", RTS="variable" )
  expect_equal( length(dea$thetaOpt), length(firms) )
  correct_thetaOpt = c(0.8838455, 0.8005349, 0.9447536, 0.8191386, 0.2893051, 0.7663874, 0.6920300, 0.6376946, 0.8477847, 0.6122914)
  expect_equal( dea$thetaOpt, correct_thetaOpt, tolerance=1e-5 )
  expect_equal( dea$lambda_sum, rep.int(1, length(firms)) )
})

test_that("output DEA with constant RTS works", {
  dea = dea( XREF=X, YREF=Y, X=X[firms,], Y=Y[firms,], model="output", RTS="constant" )
  expect_equal( length(dea$thetaOpt), length(firms) )
  correct_thetaOpt   = c(0.8762136, 0.7860619, 0.9400625, 0.7466703, 0.2705233, 0.7560322, 0.6860306, 0.6304803, 0.8434420, 0.5824253)
  correct_lambda_sum = c(2.6999022, 3.3310918, 1.8208216, 0.2442531, 0.7877087, 6.0126946, 2.7271697, 0.3262517, 3.1428140, 0.7935276)
  expect_equal( dea$thetaOpt, correct_thetaOpt, tolerance=1e-5 )
  expect_equal( dea$lambda_sum, correct_lambda_sum, tolerance=1e-5 )
})

test_that("output DEA with non-increasing RTS works", {
  dea = dea( XREF=X, YREF=Y, X=X[firms,], Y=Y[firms,], model="output", RTS="non-increasing" )
  expect_equal( length(dea$thetaOpt), length(firms) )
  correct_thetaOpt   = c(0.883845535216761, 0.800534940677135, 0.944753605800893, 0.746670299529581, 0.270523322940735, 0.766387433133089, 0.692030034384703, 0.630480263517832, 0.847784674423146, 0.582425261489661)
  correct_lambda_sum = c(1, 1, 1, 0.244253138500423, 0.787708703374778, 1, 1, 0.326251709829963, 1, 0.793527597023825)
  expect_equal( dea$thetaOpt, correct_thetaOpt, tolerance=1e-5 )
  expect_equal( dea$lambda_sum, correct_lambda_sum, tolerance=1e-5 )
})


##################################################################

context("Costmin DEA")

test_that("costmin DEA with variable RTS works", {
  dea = dea( XREF=X, YREF=Y, X=X[firms,], Y=Y[firms,], W=W[firms,], model="costmin", RTS="variable" )
  expect_equal( length(dea$gammaOpt), length(firms) )
  correct_gammaOpt = c(0.075789159498322, 0.0763752531357781, 0.0456014736519581, 0.34132473305433, 0.0960008610681884, 0.311276863261801, 0.0379973995136666, 0.0821201808600898, 0.584907590215816, 0.212884606360443)
  correct_XOpt1 = c(1.63819095477387, 7.57829839704069, 1.52763819095477, 3, 1, 100.225134716623, 1.85951928144261, 2.58982035928144, 307.825474036726, 3)
  correct_XOpt2 = c(227.638190954774, 307.901356350185, 205.527638190955, 51, 100, 380.143005344121, 183.625094516876, 91.6077844311377, 464.860759020028, 51)
  expect_equal( dea$gammaOpt, correct_gammaOpt, tolerance=1e-5 )
  expect_equal( dea$lambda_sum, rep.int(1, length(firms)) )
  expect_equal( dea$XOpt[,1], correct_XOpt1, tolerance=1e-5 )
  expect_equal( dea$XOpt[,2], correct_XOpt2, tolerance=1e-5 )
})

test_that("costmin DEA with constant RTS works", {
  dea = dea( XREF=X, YREF=Y, X=X[firms,], Y=Y[firms,], W=W[firms,], model="costmin", RTS="constant" )
  expect_equal( length(dea$gammaOpt), length(firms) )
  correct_gammaOpt = c(0.074594384391326, 0.0405556342136567, 0.0434047401964453, 0.148724677698832, 0.0157404430859198, 0.0735575478385203, 0.0342464786836621, 0.0603398580980251, 0.11643628270966, 0.0994976197003038)
  correct_XOpt1 = c(1.51020408163265, 2.04081632653061, 1.36054421768707, 0.417204301075269, 0.131764705882353, 2.44897959183673, 1.20705882352941, 0.630588235294118, 6.45591397849462, 0.296470588235294)
  correct_XOpt2 = c(226.530612244898, 306.122448979592, 204.081632653061, 31.2903225806452, 19.7647058823529, 367.34693877551, 181.058823529412, 94.5882352941177, 484.193548387097, 44.4705882352941)
  correct_lambda_sum = c(0.755102040816327, 1.02040816326531, 0.680272108843537, 0.208602150537634, 0.0658823529411765, 1.22448979591837, 0.603529411764706, 0.315294117647059, 3.22795698924731, 0.148235294117647)
  expect_equal( dea$gammaOpt, correct_gammaOpt, tolerance=1e-5 )
  expect_equal( dea$lambda_sum, correct_lambda_sum, tolerance=1e-5 )
  expect_equal( dea$XOpt[,1], correct_XOpt1, tolerance=1e-5 )
  expect_equal( dea$XOpt[,2], correct_XOpt2, tolerance=1e-5 )
})

test_that("costmin DEA with non-increasing RTS works", {
  dea = dea( XREF=X, YREF=Y, X=X[firms,], Y=Y[firms,], W=W[firms,], model="costmin", RTS="non-increasing" )
  expect_equal( length(dea$gammaOpt), length(firms) )
  correct_gammaOpt = c(0.074594384391326, 0.0763752531357779, 0.0434047401964453, 0.148724677698832, 0.0157404430859198, 0.311276863261801, 0.0342464786836621, 0.0603398580980251, 0.584907590215816, 0.0994976197003038)
  correct_XOpt1 = c(1.51020408163265, 7.57829839704067, 1.36054421768707, 0.417204301075269, 0.131764705882353, 100.225134716623, 1.20705882352941, 0.630588235294118, 307.825474036726, 0.296470588235294)
  correct_XOpt2 = c(226.530612244898, 307.901356350185, 204.081632653061, 31.2903225806452, 19.7647058823529, 380.143005344121, 181.058823529412, 94.5882352941177, 464.860759020028, 44.4705882352941)
  correct_lambda_sum = c(0.755102040816327, 1, 0.680272108843537, 0.208602150537634, 0.0658823529411765, 1, 0.603529411764706, 0.315294117647059, 1, 0.148235294117647)
  expect_equal( dea$gammaOpt, correct_gammaOpt, tolerance=1e-5 )
  expect_equal( dea$lambda_sum, correct_lambda_sum, tolerance=1e-5 )
  expect_equal( dea$XOpt[,1], correct_XOpt1, tolerance=1e-5 )
  expect_equal( dea$XOpt[,2], correct_XOpt2, tolerance=1e-5 )
})

test_that("costmin DEA with 1 dimensional inputs", {
  Y1 = hospitals[c('inpatients', 'outpatients')]
  X1 = hospitals[c('labor')]
  W1 = hospitals[c('labor_price')]
  firms = 1:10
  dea = dea( XREF=X1, YREF=Y1, X=X1[firms,], Y=Y1[firms,], W=W1[firms,], model="costmin", RTS="non-increasing" )
})

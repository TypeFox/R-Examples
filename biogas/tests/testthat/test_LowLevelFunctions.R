# test values form the lower level functions
# Charlotte
# 22 July 2015
# Modified 10 April 2016 SDH

# NTS No test yet for rbindf.R and partCO2.R

context( "Tests low level and hidden functions")

# volume molar values:
vmch4 <- 22360.588
vmco2 <- 22263.009

# Functions test with no dependencies 
#------------------------------------

# interp.R
# NTS maybe need more test here on error messages
test_that("interp interpolate and extrapolate corrrectly",{
  expect_equal(interp(c(1, 3 , 4), y = c(1, 3, 4), time.out = 2), c(t2 = 2))
  expect_equal(interp(c(1, 3 , 4), y = c(1, 3, 4), time.out = 0, extrap = TRUE), c(t0 = 1))
  expect_equal(interp(c(1, 3 , 4), y = c(1, 3, 4), time.out = 0, extrap = TRUE, method = 'hyman'), c(t0 = 0))
})
  
# unitConvert.R
# NTS can be improved with more tests for each units
test_that("unitConvet conversion is stable",{
  expect_equal(unitConvert(1, 'atm', to = 'Pa'), 101325)
  expect_equal(unitConvert(0, 'C', to = 'K'), 273.15)
})

# readFormula.R
test_that("readFormula reads correctly the formulas independement from case sensitivity and double letter elements",{
  expect_equal(readFormula('c6h12o6'), c(c = 6, h = 12, o = 6))
  expect_equal(readFormula('C6H12O6'), c(C = 6, H = 12, O = 6))
  expect_equal(readFormula('NaCl'), c(Na = 1, Cl = 1))
  expect_equal(readFormula('CHONNaClSPK'), c(C = 1, H = 1, O = 1, N = 1, Na = 1, Cl = 1, S = 1, P = 1, K = 1))
})

# watVap.R
test_that("results from watVAp is stable", {
  # water vapor at 20 C, 1 atm 
  wV <- 10^(11.20963 - 2354.731/(293.15 + 7.559))
  expect_equal(watVap(293.15),wV)
})

# Functions test with dependencies 
#------------------------------------

# molMass.R -- depends on readFormula
test_that("molMass results with one letter element are correct", {
  mmcellu <- 12.01*6 + 1.008*12 + 16.00*6
  expect_equal(molMass('C6H12O6'), mmcellu, tolerance = 0.01)
})
test_that("molMass results are not influencend by case sensitivity", {
  expect_equal(molMass('c6h12o6'), molMass('C6H12O6'))
})
test_that("molMass results is correct for non organic formula - two letter elements", {
  mmNaCl <- 22.990 + 35.45
  expect_equal(molMass('NaCl'), mmNaCl, tolerance = 0.01)
})
test_that("molMass results is correct for formula mixing one and two letter elements", {
  mm <-  12.01 + 1.008 + 16.00 + 14.007 + 32.06 + 30.974 + 22.990 + 39.098 + 35.45
  expect_equal(molMass('CHONNaClSPK'), mm, tolerance = 0.01)
})
test_that("molMass results is correct for formula with parentheses", {
  mm <-  
  expect_equal(molMass('FeSO4(H2O)7'), 278.0146, tolerance = 0.01)
})

# calcCOD.R -- depends on molMass , readFormula 
# Not sure if the best approach here. molmass cellu as hard value or calculated with molmass
test_that("calcCOD results with one letter elelement are correct ", {
  CODcellu <- (2*6 + 0.5*12 - 1.5*0 - 6)*16.00 /molMass('C6H12O6')
  expect_equal(calcCOD('C6H12O6'), CODcellu,5)
})
test_that("calcCOD results are not influenced by case sensitivity ", {
  expect_equal(calcCOD('c6h12o6'), calcCOD('C6H12O6'))
})
test_that("calcCOD result is null for non organic formula - two letter elements", {
  expect_equal(calcCOD('NaCl'), 0)
})

# stdVol.R --- depends on watVap  
# NTS : could add tests with conversion of units
test_that("standard volume calculation is stable", {
  # convert 100ml at 35 C, rh = 1, 1 atm to 0 C, rh = 0, 1 atm
  rh <- 1
  pH2O <- rh*watVap(273.15 + 35) 
  dv <- 100*(101325 - pH2O)/101325
  stdv <- dv*273.15 /  (35+273.15)
  expect_equal(stdVol(100, temp = 35, pres = 1, unit.temp = 'C', unit.pres = 'atm'), stdv)
})

# vol2mass.R --- depends on stdVol, watVap, molMass
test_that("vol2mass results are stables", {
  volBg <- stdVol(100, temp = 273.15 + 20, pres = 101325, rh = 1, temp.std = 273.15, pres.std = 101325, unit.temp = 'K', unit.pres = 'Pa', std.message = FALSE)
  mmb <- 0.65*molMass('CH4') + (1 - 0.65)*molMass('CO2')
  mvb <- 0.65*vmch4 + (1 - 0.65)*vmco2
  db <- mmb/mvb
  pH2O <- 1*watVap(temp.k = 273.15 + 35) 
  mH2O <- molMass('H2O')*pH2O/(((1.5 *101325) - pH2O)*mvb)
  mass <- volBg*(db + mH2O)

  expect_equal(vol2mass(100, xCH4 = 0.65, temp.hs = 35, temp.vol = 20, pres.hs = 1.5, pres.vol = 1), mass)
})

# mass2vol.R --- depends on molMass, watVap 
# NTS : PB with the test vlues, I don't get why -- I must make a mistake somewhere
test_that("mass2vol non vectorised result is correct", {
  # source('R/watVap.R')
 # Volume of methane if measured mass loss was 3.1 g at 35 C, 1atm, 100% humidity and 65% of methane
  mmb <- 0.65*molMass('CH4') + (1 - 0.65)*molMass('CO2')
  mvBg <- 0.65*vmch4 + (1 - 0.65)*vmco2
  db <- mmb/mvBg
  pH2O <- 1*watVap(temp.k = (273.15 + 35))
  mH2O <- molMass('H2O')*pH2O/((101325 - pH2O)*mvBg)
  # Biogas volume
  VBg <- 3.1/(db + mH2O)
  stdVCH4 <- VBg*0.65*vmch4/mvBg
  
  expect_equal(mass2vol(3.1, xCH4 = 0.65, temp = 35, pres = 1, unit.temp = 'C', unit.pres = 'atm'), stdVCH4 )
})



test_that("test Yahoo", {
  x <- pdfetch_YAHOO(c("^gspc","^ixic"), "adjclose")
  x <- pdfetch_YAHOO(c("^gspc","^ixic"))
  x <- pdfetch_YAHOO(c("vwo"))
})

test_that("test FRED", {
  x <- pdfetch_FRED(c("MDOTHFRAFDICTP1T4FR","UNRATE","GDPCA"))
})

test_that("test ECB", {
  x <- pdfetch_ECB("FM.B.U2.EUR.4F.KR.DFR.CHG")
  x <- pdfetch_ECB("IEAQ.Q.I6.N.V.B10.Z.S1M.A1.S.1.X.E.Z")
  x <- pdfetch_ECB("EXR.H.AUD.EUR.SP00.A")
  x <- pdfetch_ECB(c("EXR.D.E1.EUR.EN00.A", "IEAQ.Q.I6.N.V.B10.Z.S1M.A1.S.1.X.E.Z"))
})

test_that("test Eurostat", {
  pdfetch_EUROSTAT_DSD("namq_gdp_c")
  pdfetch_EUROSTAT_DSD("cdh_e_fos")
  pdfetch_EUROSTAT_DSD("irt_euryld_d")
  x <- pdfetch_EUROSTAT("cdh_e_fos", FREQ="A", Y_GRAD="TOTAL", FOS07=c("FOS1","FOS2"))
  x <- pdfetch_EUROSTAT("namq_gdp_c", FREQ="Q", S_ADJ="SWDA", UNIT="MIO_EUR", INDIC_NA="B1GM", GEO=c("DE","UK"))
  x <- pdfetch_EUROSTAT("irt_euryld_d", startPeriod=as.Date("2014-01-15"), MATURITY="Y1", FREQ="D", CURV_TYP="YCSR_RT")
  
  
})

test_that("test World Bank", {
  x <- pdfetch_WB("NY.GDP.MKTP.CD", c("BR","MX"))
})

test_that("test Bank of England", {
  x <- pdfetch_BOE(c("LPMVWYR", "LPMVWYR"), "2012-01-01")
})

test_that("test US Bureau of Labor Statistics", {
  x <- pdfetch_BLS(c("EIUIR","EIUIR100"), 1990, 2014)
  x <- pdfetch_BLS(c("LUU0203161800"), 2000, 2010)
  x <- pdfetch_BLS(c("ENU0100110010"), 2000, 2014)
  x <- pdfetch_BLS(c("BDS0000000000000000110001LQ5"), 2000, 2010)
})

test_that("test INSEE", {
  x <- pdfetch_INSEE(c("000810635"))
  x <- pdfetch_INSEE(c("001625866","001625866x","001616357"))
  x <- pdfetch_INSEE(c("001625866x"))
})

test_that("test ONS", {
  x <- pdfetch_ONS(c("K5CB"), "emp")
  x <- pdfetch_ONS(c("K5BZ","K54L"), "emp")
  x <- pdfetch_ONS(c("LF24","LF2G"), "lms")
})

test_that("test EIA", {
  x <- pdfetch_EIA(c("ELEC.GEN.ALL-AK-99.A","ELEC.GEN.ALL-AK-99.Q"), EIA_KEY)
  x <- pdfetch_EIA(c("PET.EMM_EPM0_PTE_NUS_DPG.M"), EIA_KEY)
  x <- pdfetch_EIA(c("PET.EMM_EPM0_PTE_NUS_DPG.W"), EIA_KEY)
  x <- pdfetch_EIA(c("PET.RWTC.D"), EIA_KEY)
})

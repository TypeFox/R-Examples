context("data_selection performance")

data(ecdata)
data(mfidata)
mfi <- mfidata[mfidata$plate=="plate_1",]
mfi <- mfi[mfi$sample!="Background0",]
names(mfi) <- c("platos","celda","analito","muestra","fluor") 

test_that("no input error", {
  expect_that(data_selection(NA), throws_error())
})

test_that("bad input bad output", {
  expect_that(data_selection(mfi, ecdata), throws_error())
})


test_that("no background", {
  aux <-  subset(mfidata, sample!="Background0")
  aux <- aux[aux$plate=="plate_1" & aux$analyte=="FGF",]
  expect_that(data_selection(aux, ecdata)[[1]], is_a("list"))
})

test_that("no controls", {
  aux <-  subset(mfidata, sample%nin%c("Control1","Control2","Control3"))
  aux <- aux[aux$plate=="plate_1" & aux$analyte=="FGF",]
  expect_that(data_selection(aux, ecdata)[[1]], is_a("list"))
})


test_that("no standard", {
  aux <-  subset(mfidata, sample%nin%c(paste0("Standard",1:17)))
  aux <- aux[aux$plate=="plate_1" & aux$analyte=="FGF",]
  expect_that(data_selection(aux, ecdata)[[1]], is_a("list"))
})


test_that("no unknowns", {
  aux <-  subset(mfidata, sample%nin%c(paste0("sid_",1:1000)))
  aux <- aux[aux$plate=="plate_1" & aux$analyte=="FGF",]
  
  expect_that(data_selection(aux, ecdata)[[1]], is_a("list"))
})


test_that("nothing in dataframe", {
  aux <-  subset(mfidata, sample!="Background0")
  aux <-  subset(aux, sample%nin%c("Control1","Control2","Control3"))
  aux <-  subset(aux, sample%nin%c(paste0("Standard",1:17)))
  aux <-  subset(aux, sample%nin%agrep("sid",aux$sample, value=TRUE))
  aux <- aux[aux$plate=="plate_1" & aux$analyte=="FGF",]
  
  expect_that(data_selection(aux, ecdata)[[1]], throws_error())
})


# test_that("flagdata only one batch", {
#  flagfile <- mfidata[mfidata$well=="P1_F10" & mfidata$plate=="plate_1" & mfidata$analyte=="IL12",]
#  expect_that(data_selection(mfidata, ecdata, flagfile)[[1]], throws_error())
# })

test_that("flagdata bad variables link", {
  flagfile <- mfidata[mfidata$well=="P1_F10" & mfidata$plate=="plate_1" & mfidata$analyte=="IL12",]
  aux <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
  expect_that(data_selection(aux, ecdata, flagfile, byvar.flagsfile=c("plate","well", "analyte"))[[1]], 
  gives_warning())
})

test_that("flagdata no variables link", {
  flagfile <- mfidata[mfidata$well=="P1_F10" & mfidata$plate=="plate_1" & mfidata$analyte=="IL12",]
  aux <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
  expect_that(data_selection(aux, ecdata, flagfile, byvar.flagsfile=NULL)[[1]], throws_error())
})

test_that("flagdata variables link", {
  flagfile <- mfidata[mfidata$well=="P1_F10" & mfidata$plate=="plate_1" & mfidata$analyte=="IL12",]
  flagfile$flag <- "bad"
  aux <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="IL12",]  
  expect_that( data_selection(aux, ecdata, flagfile, byvar.flagsfile=c("plate","well","analyte"))[[1]], is_a("list"))
})

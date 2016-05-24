## library(TR8)
## context("tr8 results")

## test_that("ecoflora works", {

##     #skip_on_cran()
##     ecoflora<-tr8(c("Salix alba"),download_list=c("li_form"))
##     eco<-ecoflora@results$li_form
##     expect_identical(as.character(eco[1]),"phanerophyte" )
## })


## test_that("BiolFlor works", {
##     #skip_on_cran()
##     biol<-tr8(c("Salix alba"),download_list=c("li_form_B"))
##     bio<-biol@results$li_form_B

##     expect_identical(as.character(bio[1]),"M (Macrophanerophyte)" )
## })


## test_that("LEDA works", {
##     #skip_on_cran()
##     LEDA<-tr8(c("Salix alba"),download_list=c("growth_form"))
##     led<-LEDA@results$growth_form
##     expect_identical(as.character(led[1]),"Phanerophyte" )
## })


## test_that("Pignatti works", {
##     #skip_on_cran()
##     Pign<-tr8(c("Salix alba"),download_list=c("ell_L_it"))
##     Pign<-Pign@results$ell_L_it
##     expect_identical(as.character(Pign[1]),"5" )
## })



## test_that("Akhmetzhanova works", {
##     #skip_on_cran()
##     Akh<-tr8(c("Abies alba"),download_list=c("Myco_infection"))
##     Akh<-Akh@results$Myco_infection
##     expect_identical(as.character(Akh[1]),"7-NA" )
## })






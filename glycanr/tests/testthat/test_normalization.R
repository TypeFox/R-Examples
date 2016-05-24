library(stringr)
library(dplyr)
context("Normalization functions")

#library(dplyr)
#library(tidyr)
#library(glycanr)
#data(mpiu)
#
#source("../../..//normalizationFunctions_23062015.r")
#
#mw <- mpiu %>% 
#    spread(glycan, value)
#
##
## tanorm
##
#tmp1 <- totalAreaNorm(mw) %>% 
#    select(-normalization)
#
#tmp2 <- tanorm(mpiu) %>% 
#    spread(glycan, value)
#tmp2 <- data.frame(tmp2)
#
#print(all.equal(tmp1, tmp2))
#if(all.equal(tmp1, tmp2)){
#    mpiunorm <- tmp2 %>% 
#        mutate(norm="tanorm")
#}
#    
#
#    
##
## refPeakNorm
##
#tmp1 <- refPeakNorm(mw) %>% 
#    select(-normalization)
#
#tmp2 <- refpeaknorm(mpiu) %>% 
#    spread(glycan, value)
#tmp2 <- data.frame(tmp2)
#
#print(all.equal(tmp1, tmp2))
#if(all.equal(tmp1, tmp2)){
#    mpiunorm <- bind_rows(mpiunorm,
#                          tmp2 %>% mutate(norm="refpeaknorm"))
#}
#    
#    
##
## mediannorm
##
#tmp1 <- medianNorm(mw) %>% 
#    select(-normalization)
#
#tmp2 <- mediannorm(mpiu) %>% 
#    spread(glycan, value)
#tmp2 <- data.frame(tmp2)
#
#print(all.equal(tmp1, tmp2))
#if(all.equal(tmp1, tmp2)){
#    mpiunorm <- bind_rows(mpiunorm,
#                          tmp2 %>% mutate(norm="mediannorm"))
#}
#
##
## medianQuotientNorm
##
#tmp1 <- medianQuotientNorm(mw) %>% 
#    select(-normalization)
#
#tmp2 <- medianquotientnorm(mpiu) %>% 
#    spread(glycan, value)
#tmp2 <- data.frame(tmp2)
#
#print(all.equal(tmp1, tmp2))
#if(all.equal(tmp1, tmp2)){
#    mpiunorm <- bind_rows(mpiunorm,
#                          tmp2 %>% mutate(norm="medianquotientnorm"))
#}
#
##
## quantileNorm
##
#library(preprocessCore)
#tmp1 <- quantileNorm(mw) %>% 
#    select(-normalization)
#
#tmp2 <- quantilenorm(mpiu, transpose=TRUE) %>% 
#    spread(glycan, value)
#tmp2 <- data.frame(tmp2)
#tmp2 <- tmp2[ , c(2,1,3:ncol(tmp2))]
#
#print(all.equal(tmp1, tmp2))
#if(all.equal(tmp1, tmp2)){
#    mpiunorm <- bind_rows(mpiunorm,
#                          tmp2 %>% mutate(norm="quantilenormt"))
#}
#
#mpiunorm <- mpiunorm %>% 
#    gather(glycan, value, contains("GP")) %>% 
#    arrange(norm, Plate, gid, glycan)
#mpiunorm <- data.frame(mpiunorm)
#save(mpiunorm, file="./mpiunorm.RData")

if(.Machine$sizeof.pointer==8){

    data(mpiu)
    data(mpiunorm)


    tanormdf <- mpiunorm %>% 
        filter(norm=="tanorm") %>% 
        select(-norm)

    refpeaknormdf <- mpiunorm %>% 
        filter(norm=="refpeaknorm") %>% 
        select(-norm)

    mediannormdf <- mpiunorm %>% 
        filter(norm=="mediannorm") %>% 
        select(-norm)

    medianquotientnormdf <- mpiunorm %>% 
        filter(norm=="medianquotientnorm") %>% 
        select(-norm)

    quantilenormtdf <- mpiunorm %>% 
        filter(norm=="quantilenormt") %>% 
        select(-norm)

    test_that("tanorm", {
      expect_equal(tanorm(mpiu), tanormdf)
    })


    test_that("refpeaknorm", {
      expect_equal(refpeaknorm(mpiu), refpeaknormdf)
    })


    test_that("mediannorm", {
      expect_equal(mediannorm(mpiu), mediannormdf)
    })


    test_that("medianquotientnorm", {
      expect_equal(medianquotientnorm(mpiu), medianquotientnormdf)
    })


    tmp <- quantilenorm(mpiu, transpose=TRUE)
    # reorder columns
    tmp <- tmp[, names(quantilenormtdf)] %>% 
        arrange(Plate, gid, glycan)
        
    test_that("quantilenormt", {
      expect_equal(tmp, quantilenormtdf)
    })

}

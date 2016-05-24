library(testthat)
library(urlshorteneR)

if (interactive()) {
  test_check("urlshorteneR")
  test_package("urlshorteneR")
}
# bitly_token <- bitly_auth(key = "be03aead58f23bc1aee6e1d7b7a1d99d62f0ede8", secret = "b7e4abaf8b26ec4daa92b1e64502736f5cd78899")
# googl_token <- googl_auth(key = "806673580943-78jdskus76fu7r0m21erihqtltcka29i.apps.googleusercontent.com", secret = "qItL-PZnm8GFxUOYM0zPVr_t")
# 
# saveRDS(bitly_token, file = "tests/testthat/bitly_token.rds")
# saveRDS(googl_token, file = "tests/testthat/googl_token.rds")


## https://gorails.com/setup/ubuntu/15.04
## https://github.com/jennybc/googlesheets/blob/master/internal-projects/10_token-encryption.Rmd

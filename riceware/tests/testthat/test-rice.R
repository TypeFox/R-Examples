context("check that data set is loaded/available")
test_that("check correct dimensions", {
          expect_equal(nrow(wordlist_en), 7776)
          expect_equal(nrow(wordlist_de), 7776)
          expect_equal(nrow(wordlist_es), 7776)
          expect_equal(nrow(wordlist_fr), 7776)
          expect_equal(nrow(wordlist_it), 7776)
          expect_equal(nrow(wordlist_jp), 7776)
          ##expect_equal(nrow(wordlist_nl), 7776)
          expect_equal(nrow(wordlist_sv), 7776)
      })

context("test token")
test_that("check length", {
          expect_false(check_token(character(0)))
          expect_false(check_token(c("11111", "22222")))
          expect_false(check_token("111"))
          expect_false(check_token("111111"))
          expect_true(check_token("11111"))
      })

test_that("check digits", {
          expect_false(check_token("A1111"))
          expect_false(check_token(NA))
          expect_false(check_token(NULL))
          expect_false(check_token("01234"))
          expect_true(check_token("12345"))
          expect_true(check_token("23456"))
      })

not_working <- function() {
    if (capabilities()["http/ftp"]) { ## Are we on a platform that allows internet access?
        if (.Platform$OS.type == "unix") { ## If unix, then we can use nsl to check connection
            res <- is.null(utils::nsl("random.org"))
        } else { ## otherwise just fail
            res <- TRUE
        }
    } else {
        res <- TRUE
    }
    res
}
check_random_org <- function() {
    if(not_working()) {
        skip("can't connect to random.org")
    }
}

context("generate token")
test_that("correct length and pseudorandom", {
          test_len <- lapply(1:10, generate_token, method = "pseudo")
          expect_equal(1:10, sapply(test_len, length))
          expect_equal(1, length(generate_token(1)))
      })

test_that("random numbers", {
              skip_on_cran()
              check_random_org()
              expect_equal(1, length(generate_token(n_words = 1, method = "random")))
          })


context("match token")
test_that("correct matching", {
              expect_equal(match_token("11111", wordlist = wordlist_en, title_case = TRUE), "A")
              expect_equal(match_token("11111", wordlist = wordlist_en, title_case = FALSE), "a")
              expect_equal(match_token("16234", wordlist = wordlist_en, title_case = TRUE), "Cat")
              expect_equal(match_token("16234", wordlist = wordlist_en, title_case = FALSE), "cat")
          })

test_that("fails if incorrect token", {
              expect_error(match_token("cat"), "invalid token")
          })

context("generate passphrase")
test_that("test passphrase", {
           expect_message(generate_passphrase(tokens = c("36156", "35646"),
                                              wordlist = wordlist_en,
                                              title_case = FALSE,
                                              verbose = TRUE), "lava lamp")
            expect_equal(generate_passphrase(tokens = c("36156", "35646"),
                                             wordlist = wordlist_en,
                                             title_case = FALSE,
                                             verbose = FALSE), "lavalamp")
            expect_equal(generate_passphrase(tokens = c("36156", "35646"),
                                             wordlist = wordlist_en,
                                             title_case = TRUE,
                                             verbose = FALSE), "LavaLamp")
           ## German
            expect_equal(generate_passphrase(tokens = c("34454"), wordlist = wordlist_de,
                                             title_case = TRUE, verbose = FALSE), "Katze")
            expect_equal(generate_passphrase(tokens = c("34454"), wordlist = wordlist_de,
                                             title_case = FALSE, verbose = FALSE), "katze")
           ## Spanish
            expect_equal(generate_passphrase(tokens = c("35622"), wordlist = wordlist_es,
                                             title_case = TRUE, verbose = FALSE), "Gato")
            expect_equal(generate_passphrase(tokens = c("35622"), wordlist = wordlist_es,
                                             title_case = FALSE, verbose = FALSE), "gato")
           ## French
            expect_equal(generate_passphrase(tokens = c("21631"), wordlist = wordlist_fr,
                                             title_case = TRUE, verbose = FALSE), "Chaton")
            expect_equal(generate_passphrase(tokens = c("21631"), wordlist = wordlist_fr,
                                             title_case = FALSE, verbose = FALSE), "chaton")
           ## Italian (no cat :( )
           expect_equal(generate_passphrase(tokens = c("32141"), wordlist = wordlist_it,
                                             title_case = TRUE, verbose = FALSE), "Gelato")
           expect_equal(generate_passphrase(tokens = c("32141"), wordlist = wordlist_it,
                                             title_case = FALSE, verbose = FALSE), "gelato")
           ## Japanese
           expect_equal(generate_passphrase(tokens = c("44565"), wordlist = wordlist_jp,
                                             title_case = TRUE, verbose = FALSE), "Neko")
           expect_equal(generate_passphrase(tokens = c("44565"), wordlist = wordlist_jp,
                                             title_case = FALSE, verbose = FALSE), "neko")
           ## Dutch
           #expect_equal(generate_passphrase(tokens = c("53431"), wordlist = wordlist_nl,
           #                                  title_case = TRUE, verbose = FALSE), "Kat")
           #expect_equal(generate_passphrase(tokens = c("53431"), wordlist = wordlist_nl,
           #                                  title_case = FALSE, verbose = FALSE), "kat")

           ## Swedish
           expect_equal(generate_passphrase(tokens = c("33343"), wordlist = wordlist_sv,
                                             title_case = TRUE, verbose = FALSE), "Katt")
           expect_equal(generate_passphrase(tokens = c("33343"), wordlist = wordlist_sv,
                                             title_case = FALSE, verbose = FALSE), "katt")
        })

## ----setup, echo = FALSE, message = FALSE--------------------------------
knitr::opts_chunk$set(tidy = FALSE, comment = "#>")

## ----setseed1, echo = FALSE----------------------------------------------
set.seed(10)

## ----quick-random--------------------------------------------------------
library(lightsout)
random1 <- random_board(3)
random1

random2 <- random_board(5)
random2

## ----quick-newboard------------------------------------------------------
lights_vector <- c(0, 0, 0,
                   0, 1, 0,
                   1, 1, 1)
board1 <- new_board(lights_vector)
board1

## ----solvable------------------------------------------------------------
is_solvable(board1)

## ----play1---------------------------------------------------------------
# Press light at (2,3)
played <- play(board1, 2, 3)
played

## ----variant-------------------------------------------------------------
# Press light at (2,3) in non-classic mode
new_board(lights_vector, classic = FALSE) %>% play(2, 3)

## ----nosolution----------------------------------------------------------
new_board(lights_vector, classic = FALSE) %>% is_solvable()

## ----pressmultiple-------------------------------------------------------
# Press light at (2,3) and then (1,2)
board1 %>% play(2, 3) %>% play(1, 2)

## ----pressmultiple2------------------------------------------------------
# Press light at (2,3) and (1,2)
board1 %>% play(c(2, 1), c(3, 2))

## ----pressmatrix---------------------------------------------------------
# Press light at (2,3) and (1,2)
board1 %>% play(matrix = matrix(nrow = 3, byrow = TRUE,
                                c(0, 1, 0,
                                  0, 0, 1,
                                  0, 0, 0)))

## ----solve1--------------------------------------------------------------
play(board1, 3, 2)

## ----setseed2, echo = FALSE----------------------------------------------
set.seed(15)

## ----board5--------------------------------------------------------------
bigboard <- random_board(size = 5)
bigboard

## ----solvebig------------------------------------------------------------
solution <- solve_board(bigboard)
solution

## ----verify--------------------------------------------------------------
bigboard %>% play(matrix = solution)


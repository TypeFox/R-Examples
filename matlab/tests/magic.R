###
### $Id: magic.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.magic <- function(input, expected) {
    output <- do.call(getFromNamespace("magic", "matlab"), input)
    identical(output, expected)
}

magic.expected.3 <- matrix(c(8, 1, 6,
                             3, 5, 7,
                             4, 9, 2),
                           nrow = 3, ncol = 3, byrow = TRUE)

magic.expected.4 <- matrix(c(16,  2,  3, 13,
                              5, 11, 10,  8,
                              9,  7,  6, 12,
                              4, 14, 15,  1),
                           nrow = 4, ncol = 4, byrow = TRUE)

magic.expected.6 <- matrix(c(35,  1,  6, 26, 19, 24,
                              3, 32,  7, 21, 23, 25,
                             31,  9,  2, 22, 27, 20,
                              8, 28, 33, 17, 10, 15,
                             30,  5, 34, 12, 14, 16,
                              4, 36, 29, 13, 18, 11),
                           nrow = 6, ncol = 6, byrow = TRUE)

magic.expected.10 <- matrix(c(92, 99,  1,  8, 15, 67, 74, 51, 58, 40,
                              98, 80,  7, 14, 16, 73, 55, 57, 64, 41,
                               4, 81, 88, 20, 22, 54, 56, 63, 70, 47,
                              85, 87, 19, 21,  3, 60, 62, 69, 71, 28,
                              86, 93, 25,  2,  9, 61, 68, 75, 52, 34,
                              17, 24, 76, 83, 90, 42, 49, 26, 33, 65,
                              23,  5, 82, 89, 91, 48, 30, 32, 39, 66,
                              79,  6, 13, 95, 97, 29, 31, 38, 45, 72,
                              10, 12, 94, 96, 78, 35, 37, 44, 46, 53,
                              11, 18,100, 77, 84, 36, 43, 50, 27, 59),
                           nrow = 10, ncol = 10, byrow = TRUE)


test.magic(list(n = 3), magic.expected.3) # 'n' odd
test.magic(list(n = 4), magic.expected.4) # 'n' divisible by four
test.magic(list(n = 6), magic.expected.6) # 'n' even but not divisible by four
test.magic(list(n = 10), magic.expected.10)


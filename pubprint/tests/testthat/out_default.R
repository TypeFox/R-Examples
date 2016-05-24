# out.default.value
#############################################################################

# Single value
pubprint:::out.default.value(0.123456789,
                                nsmall = 2,
                                replace0 = FALSE,
                                leading0 = TRUE)
pubprint:::out.default.value(0.123456789,
                                name = "t",
                                nsmall = 2,
                                replace0 = FALSE,
                                leading0 = TRUE)
pubprint:::out.default.value(0.123456789,
                                name = "t",
                                inbracket = "N=1",
                                nsmall = 2,
                                replace0 = FALSE,
                                leading0 = FALSE)

# Number vector
pubprint:::out.default.value(c(0.1,
                                  0.234567891,
                                  0.345678912,
                                  0.456789123,
                                  0.001,
                                  -0.001),
                                nsmall = 2,
                                replace0 = TRUE,
                                leading0 = TRUE)
pubprint:::out.default.value(c(0.123456789,
                                  0.234567891,
                                  0.345678912,
                                  0.456789123),
                                name = "t",
                                nsmall = 2,
                                replace0 = FALSE,
                                leading0 = TRUE)
pubprint:::out.default.value(c(0.123456789,
                                  0.234567891,
                                  0.345678912,
                                  0.456789123,
                                  0.001,
                                  -0.001),
                                name = "t",
                                inbracket = "N=1",
                                nsmall = 2,
                                replace0 = TRUE,
                                leading0 = TRUE)
pubprint:::out.default.value(c(0.123456789,
                                  0.234567891,
                                  0.345678912,
                                  0.456789123),
                                name = c("t", 
                                         "x", 
                                         "s", 
                                         "v"),
                                nsmall = 2,
                                replace0 = FALSE,
                                leading0 = TRUE)
pubprint:::out.default.value(c(0.123456789,
                                  0.234567891,
                                  0.345678912,
                                  0.456789123),
                                name = c("t", 
                                         "x", 
                                         "s", 
                                         "v"),
                                inbracket = c("N=1",
                                              "N=2",
                                              "N=3",
                                              "N=4"),
                                nsmall = 2,
                                replace0 = FALSE,
                                leading0 = TRUE)

# out.default.concat
#############################################################################

# Single string
pubprint:::out.default.concat("Little Test", sep = ";")

# String vector
pubprint:::out.default.concat(c("Little Test", "Second Test"), sep = ";")

# Multiple arguments
pubprint:::out.default.concat("Little Test", "Second Test", sep = ";")

# Multiple vector arguments
pubprint:::out.default.concat(c("Little Test", "Second Test"), c("Test again", "and again"), sep = ";")

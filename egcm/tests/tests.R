# egcm.R  
# Copyright (C) 2014 by Matthew Clegg

#  Tests for the Engle-Granger cointegration models package

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

test <- function(expr, out="") {
    # expr is a string representing an R expression, and
    # out is the output that is expected.  Prints and evaluates
    # expr.  If out is given and it matches the output of
    # evaluating expr, returns TRUE.  Otherwise, returns FALSE.
    
    cat("\n", expr, "\n", sep="")
    if (!missing(out)) cat("Expecting:",out, "", sep="\n")
    val <- eval(parse(text=expr))
    cat("Result:\n")
    print(val)
    if (!missing(out) && 
       (as.character(val) != as.character(out))) {
        return(FALSE)
    }
    return(TRUE)
}

assert <- function (expr, out) {
    # expr is astring representing an R expression,
    # and out is the output that is expected.  Prints
    # and evaluates expr.  If out matches the output of
    # evaluating expr, returns TRUE.  Otherwise, stops
    # the execution with an error message.
    if (!test(expr, out)) {
        stop("Expression does not evaluate to its expected value\n")
    }
}

egcm_tests <- function (include.inexact = FALSE) {
    # If include.inexact is TRUE, then includes those tests for
    # which an exact answer cannot be given.  
    
    if (include.inexact) {
        test("ur_power(pp.test, nrep=100)", 0.4)
        test("ur_power_table(pp.test, nrep=10)")
        test("egcm:::ur_specificity(pp.test, function() cumsum(rnorm(100)))",
           paste("0.01  0.02  0.05   0.1   0.2   0.5",
                 " 0.01  0.02  0.05   0.1   0.2   0.5", sep="\n"))
                 
        test("adf_power(nrep=100)", 0.26)
        test("adf_power_table(nrep=10)")                 
    }
    
    assert("bvr_rho(1:100)", 0.10001)
    assert("bvr.test(1:100+sin(1:100))$p.value", 0.999)
    assert("bvr.test(1:100+sin(1:100),TRUE)$p.value", 0.001)

    if (include.inexact) {
        test("bvr_power(nrep=100)", 0.36)
        test("bvr_power_table(nrep=10)")
    }

    if (include.inexact) {
        test("pgff_rho_ws(1:100)", 1.019777)
        test("pgff_rho_ws(1:100 + sin(1:100))", 1.018697)
        test("pgff_rho_ws(1:100 + sin(1:100), TRUE)", 0.5464542)
    }        
        
    assert("pgff.test(1:100+sin(1:100))$p.value > 0.99", TRUE)
    assert("pgff.test(1:100+sin(1:100),TRUE)$p.value", 0.001)
    
    if (include.inexact) {
        test("pgff_power(nrep=100)", 0.75)
        test("pgff_power_table(nrep=10)")
    }
    
    assert("length(egcm.urtests())", 11)
    assert("egcm.default.urtest()", "pp")
    assert("egcm.set.default.urtest(\"adf\"); egcm.default.urtest()", 
        "adf")
    assert("egcm.set.default.urtest(\"pgff\"); egcm.default.urtest()", 
        "pgff")

    assert("length(egcm.i1tests())", 4)
    assert("egcm.default.i1test()", "pp")
    assert("egcm.set.default.i1test(\"adf\"); egcm.default.i1test()", 
        "adf")
    assert("egcm.set.default.i1test(\"pgff\"); egcm.default.i1test()", 
        "pgff")

    assert("egcm.default.pvalue()", 0.05)
    assert("egcm.set.default.pvalue(0.01); egcm.default.pvalue()", 0.01)
    assert("egcm.set.default.pvalue(0.05); egcm.default.pvalue()", 0.05)
    
    assert("detrend(1:5)", rep(0,5))

test("egcm(1:100, 1:100)",
"Y[i] =   1.0000 X[i] -   0.0000 + R[i], R[i] =  -0.1241 R[i-1] + eps[i], eps ~ N(0,  0.0000^2)
        (0.0000)        (0.0000)                (0.0992)

R[100] = -0.0000 (t = -0.089)")

test("egcm(1:100, (1:100)*2)",
"Y[i] =   2.0000 X[i] -   0.0000 + R[i], R[i] =  -0.1241 R[i-1] + eps[i], eps ~ N(0,  0.0000^2)
        (0.0000)        (0.0000)                (0.0992)

R[100] = -0.0000 (t = -0.089)")

test("egcm(1:100, (1:100)*2, normalize=TRUE)",
"Y[i] =   1.0000 X[i] -   0.0000 + R[i], R[i] =  -0.1241 R[i-1] + eps[i], eps ~ N(0,  0.0000^2)
        (0.0000)        (0.0000)                (0.0992)

R[100] = -0.0000 (t = -0.089)")
          
test("egcm (1:100, (1:100)*2+3)",
"Y[i] =   2.0000 X[i] +   3.0000 + R[i], R[i] =  -0.1241 R[i-1] + eps[i], eps ~ N(0,  0.0000^2)
        (0.0000)        (0.0000)                (0.0992)

R[100] = -0.0000 (t = -0.089)")
         
test("egcm (1:100, (1:100)*2+3, include.const=FALSE)",
"Y[i] =   2.0448 X[i] +   0.0000 + R[i], R[i] =   1.0000 R[i-1] + eps[i], eps ~ N(0,  0.0000^2)
        (0.0026)        (0.0000)                (0.0026)

R[100] = -1.4776 (t = -1.137)")

test("egcm(1:100, 1:100 + sin(1:100))",
"Y[i] =   0.9988 X[i] +   0.0583 + R[i], R[i] =   0.5901 R[i-1] + eps[i], eps ~ N(0,  0.5957^2)
        (0.0025)        (0.3309)                (0.0842)

R[100] = -0.4467 (t = -0.628)

WARNING: The series seem cointegrated but the residuals are not AR(1).")
        
test("egcm(1:100, 1:100 + sin(1:100), debias=FALSE)",
"Y[i] =   0.9988 X[i] +   0.0583 + R[i], R[i] =   0.5452 R[i-1] + eps[i], eps ~ N(0,  0.5949^2)
        (0.0025)        (0.3305)                (0.0842)

R[100] = -0.4467 (t = -0.628)

WARNING: The series seem cointegrated but the residuals are not AR(1).")

test("egcm(1:101, c(1:100 + sin(1:100), 110), robust=FALSE)", 
"Y[i] =   1.0041 X[i] -   0.1211 + R[i], R[i] =   0.4712 R[i-1] + eps[i], eps ~ N(0,  1.0925^2)
        (0.0039)        (0.5923)                (0.1498)

R[101] = 8.7073 (t = 7.665)

WARNING: The series seem cointegrated but the residuals are not AR(1).")

test("egcm(1:101, c(1:100 + sin(1:100), 110), robust=TRUE)",
"Y[i] =   0.9997 X[i] +   0.0284 + R[i], R[i] =   0.5770 R[i-1] + eps[i], eps ~ N(0,  1.1034^2)
        (0.0025)        (0.5710)                (0.0861)

R[101] = 9.0019 (t = 7.874)

WARNING: The series seem cointegrated but the residuals are not AR(1).")

test("egcm(1:25, 1:25+sin(1:25))",
"Y[i] =   0.9818 X[i] +   0.2349 + R[i], R[i] =   0.7471 R[i-1] + eps[i], eps ~ N(0,  0.6132^2)
        (0.0201)        (0.4284)                (0.1711)

R[25] = 0.0890 (t = 0.125)

WARNING: X and Y do not appear to be cointegrated.")

test("egcm(exp(1:25), exp(1:25+sin(1:25)), log=TRUE)",
"Y[i] =   0.9818 X[i] +   0.2349 + R[i], R[i] =   0.7471 R[i-1] + eps[i], eps ~ N(0,  0.6132^2)
        (0.0201)        (0.4284)                (0.1711)

R[25] = 0.0890 (t = 0.125)

WARNING: X and Y do not appear to be cointegrated.")

test("egcm(1:25, 1:25+sin(1:25), urtest=\"adf\")",
"Y[i] =   0.9818 X[i] +   0.2349 + R[i], R[i] =   0.7471 R[i-1] + eps[i], eps ~ N(0,  0.6132^2)
        (0.0201)        (0.4284)                (0.1711)

R[25] = 0.0890 (t = 0.125)

WARNING: The series seem cointegrated but the residuals are not AR(1).")

test("egcm(1:25, 1:25+sin(1:25), p.value=0.2)",
"Y[i] =   0.9818 X[i] +   0.2349 + R[i], R[i] =   0.7471 R[i-1] + eps[i], eps ~ N(0,  0.6132^2)
        (0.0201)        (0.4284)                (0.1711)

R[25] = 0.0890 (t = 0.125)

WARNING: X and Y do not appear to be cointegrated.")

test("try(egcm(1:100, c(1:49,NA,51:100)))",
"Y[i] =   1.0000 X[i] -   0.0000 + R[i], R[i] =  -0.0757 R[i-1] + eps[i], eps ~ N(0,  0.0000^2)
        (0.0000)        (0.0000)                (0.1004)

R[99] = 0.0000 (t = 0.001)")

test("try(egcm(1:100, c(1:49,NA,51:100), na.fail), silent=TRUE)",
"[1] \"Error in na.fail.default(S12) : missing values in object\\n\"
attr(,\"class\")
[1] \"try-error\"
attr(,\"condition\")
<simpleError in na.fail.default(S12): missing values in object>")


    assert("is.cointegrated(egcm(1:100, 1:100))", TRUE)
    assert("is.cointegrated(egcm(1:100, (1:100)^2))", FALSE)
    
    assert("is.ar1(egcm(1:100, 1:100))", TRUE)
    assert("is.ar1(egcm(1:100, 1:100 + sin(1:100)))", FALSE)
    
    assert("acor(1:10)", 1)
    assert("acor(c(1,0,1,0,1,0))", -1)
    assert("acor(c(1,0,-1,0,1,0))", 0)

    if (include.inexact) {
        test("sum(replicate(100, pgff.test(rar1(100, a1=0.5))$p.value < 0.05))", 100)
        test("sum(replicate(100, pgff.test(rar1(100, a1=1))$p.value < 0.05))", 5)
        
    }
    
    # Due to differences in machine precision, we cannot
    # test these exactly, but at least we can try different
    # combinations of parameters and see if the function
    # at least executes.
    
}

library(egcm)
egcm_tests()

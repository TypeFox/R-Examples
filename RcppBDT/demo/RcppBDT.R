
demo.RcppBDT  <- function() {

    require(RcppBDT, quiet=TRUE, warn=FALSE)

    ## this uses the pretty-printing the Rcpp module logic to show all
    ## available functions and their docstring (symbol is not exported)
    ##
    ## print(bdtMod$date)

    ## we use a base object 'bdt' for the functions below
    ## by using the instance stored in the environment
    ##
    ## alternatively could construct a new instance from bdtMod$date, see R/zzz.R

    cat("Demo of setters\n");
    ## conversions from string commented out, see inst/include/RcppBDT.h for details
    ##bdt$fromString("2010-10-02"); 	cat("From 2010-10-02   : ", format(bdt), "\n")
    ##bdt$fromUndelString("20101003");  cat("From 20101003     : ", format(bdt), "\n")
    bdt$setFromUTC(); 			cat("From curr. UTC    : ", format(bdt), "\n")
    bdt$setFromLocalClock();		cat("From curr. local  : ", format(bdt), "\n")
    bdt$setEndOfMonth(); 		cat("end of month      : ", format(bdt), "\n")
    bdt$setFirstOfNextMonth(); 		cat("1st of next Month : ", format(bdt), "\n")
    bdt$addDays(4);                     cat("plus four days    : ", format(bdt), "\n")
    bdt$subtractDays(3);                cat("minus three days  : ", format(bdt), "\n")

    bdt$setIMMDate(12, 2010); 		cat("IMM Date Dec 2010 : ", format(bdt), "\n")
    bdt$setEndOfBizWeek();  		cat("end of biz week   : ", format(bdt), "\n")

    cat("\nDemo of getters\n")
    ## now just functions that return values to R
    cat("From curr. local  : ", format(bdt$getLocalClock()), "\n")
    bdt$setFromLocalClock();
    cat("end of biz week   : ", format(bdt$getEndOfBizWeek()), "\n")
    cat("end of of month   : ", format(bdt$getEndOfMonth()), "\n")
    cat("1st of next month : ", format(bdt$getFirstOfNextMonth()), "\n")

    cat("\nDemo of functions\n")
    cat("IMM Date Dec 2010 : ", format(getIMMDate(Dec, 2010)), "\n")
    cat("3rd Wed Dec 2010  : ", format(getNthDayOfWeek(third, Wed, Dec, 2010)), "\n")
    cat("Last Sat Dec 2010 : ", format(getLastDayOfWeekInMonth(Sat, Dec, 2010)), "\n")
    cat("First Sat Dec 2010: ", format(getFirstDayOfWeekInMonth(Sat, Dec, 2010)), "\n")
    cat("First Wed in 2011 : ", format(getFirstDayOfWeekAfter(Wed, as.Date("2010-12-31"))), "\n")
    cat("Last Wed in 2010  : ", format(getLastDayOfWeekBefore(Wed, as.Date("2010-12-31"))), "\n")
}

demo.RcppBDT()

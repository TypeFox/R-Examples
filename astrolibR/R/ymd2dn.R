ymd2dn = function(yr,m,d) {

                                        #----  Days before start of each month (non-leap year)  -----
    idays = c(0,31,59,90,120,151,181,212,243,273,304,334,366)

                                        #----  Correct for leap year if month ge 3  -------------
    lpyr = (((yr %% 4) == 0) & ((yr %% 100) != 0))  | ((yr %% 400) == 0) & (m >= 3)

    dy = d + idays[m] + lpyr
    return(dy)

}

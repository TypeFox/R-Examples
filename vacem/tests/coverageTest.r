#*  ----------------------------------------------------------------------------
#*  Copyright (C) 2011-2012 - Justin Lessler, Jessica Metcalf
#*  
#*  This program is free software; you can redistribute it and/or modify
#*  it under the terms of the GNU General Public License as published by
#*  the Free Software Foundation; version 2 of the License.
#*  
#*  This program is distributed in the hope that it will be useful,
#*  but WITHOUT ANY WARRANTY; without even the implied warranty of
#*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#*  GNU General Public License for more details.
#*  
#*  You should have received a copy of the GNU General Public License
#*  along with this program; if not, write to the Free Software
#*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#*  
#*  $Id: coverageTest.r 378 2012-01-27 02:35:09Z ken $
#*  ----------------------------------------------------------------------------


# TBD: Add file header
# TBD: Do I need some 'library( pkg )' or 'source( file )' calls here???


# source( "../tests/testUtils.r" )
source( getTestPaths("vacem","tests/testUtils.r") )


# TBD: Do I need/want this 'cat' (or 'print') here?
if ( verbose ) {
  #cat( "\n" )
  cat( "Vacem test cases for coverage functions:\n" )
}


## -----------------------------------------------------------------------------
## Test Data
##
## Define some data for unit test.  This should probably be put in the
## setup function.  Or maybe put this in a separate file that we can
## just read in?

TEST.DATA.CAMPAIGNS <-
"C, date,       v,      N,      age.low, age.high, is.SIA
 1, 2003-01-01, 646166, 775191,       0,     1200,      0
 2, 2004-01-01, 660776, 793461,       0,     1200,      0
 3, 2005-01-01, 718589, 812221,       0,     1200,      0
 4, 2006-01-01, 759222, 891586,       0,     1200,      0
 5, 2007-01-01, 812083, 857899,       0,     1200,      0
"
# camps <- read.table( textConnection(TEST.DATA.CAMPAIGNS), header = TRUE )

# camps <- data.frame(
#       date = c( "2008-01-21", "2009-02-22", "2010-05-31", "2011-07-04" ),
#          N = c(        10000,        15000,        20000,        25000 ),
#          v = c(         1000,         1200,         1400,         1600 ),
#    age.low = c(            0,            0,            0,            0 ),
#   age.high = c(           36,          216,           36,          120 ),
#     is.SIA = c(            0,            1,            0,            1 )
# )
               

TEST.DATA.OBSERVATIONS <-
"C,  date,       age,   y, card.date,  age.at.vac, samp.wt
 1,  2008-10-15, 30.05, 1, 2007-04-05,  11.671232, 0.652967
 2,  2008-10-15, 12.03, 1, 2008-08-07,   9.764382, 0.652967
 3,  2008-11-15, 27.06, 1, 2007-06-15,   9.994519, 0.911245
 4,  2008-11-15, 43.07, 1,         NA,         NA, 0.911245
 5,  2008-09-15, 52.08, 1,         NA,         NA, 1.211556
"
# survey <- read.table( textConnection(TEST.DATA.OBSERVATIONS), header = TRUE )

## -----------------------------------------------------------------------------


# TBD: Add doc
# TBD: Replace verbose switch with a global logging method that
#      the caller puts into the environment and that we can retrieve
#      and write to.  The logging method will decide whether to
#      write out or ignore the message.
#      Note: test functions should not take args, so we need to
#            use the environment to pass parameters
#
test.g <- function () {
  # DEACTIVATED( "Disabling g testing while working on other tests" )
  N <- c( 1000, 1000 )  # target population of 1000 for each campaign
  v <- c(  500,  500 )  # two campaigns of 500 doses
  z <- c( 1, 1 )        # eligibility vector
  w <- c( 0, 0 )        # weight vector
  rho <- 1.0            # rho = 0.9 ==> 90% of population accessible
  toler <- 1.0e-7
  fmt <- paste( "\n",
                " Test: %s\n",
                " Params: N=%s, v=%s, z=%s, w=%s, rho=%s, alpha=%s\n",
                " Expected: %e\n",
                " Actual:   %e (Diff: %g)\n" )
  verbose <- ( getOption( "RUnit" )$verbose && FALSE )

  if ( verbose ) { cat( "\n" ) }

  name <- "perfect efficiency test (alpha = -Inf)"
  if ( verbose ) { cat( sprintf("\t... Running %s ...\n", name) ) }
  alpha <- -Inf
  expect <- 0.75
  actual <- g( z, w, v, N, alpha, rho )
  diff <- ( actual - expect )
  msg <- sprintf( fmt, name, toStr(N), toStr(v), toStr(z), toStr(w),
                       toStr(rho), toStr(alpha), expect, actual, diff )
  checkEquals( actual, expect, msg, toler )

  name <- "random efficiency test (alpha = 0)"
  if ( verbose ) { cat( sprintf("\t... Running %s ...\n", name) ) }
  alpha <- 0
  expect <- 0.6321206
  actual <- g( z, w, v, N, alpha, rho )
  diff <- ( actual - expect )
  msg <- sprintf( fmt, name, toStr(N), toStr(v), toStr(z), toStr(w),
                       toStr(rho), toStr(alpha), expect, actual, diff )
  checkEquals( actual, expect, msg, toler )
  
  if ( verbose ) { cat( " ... " ) }


  # TBD: Add check to confirm 0 <= g <= 1

}  # end of function test.g()



# TBD: Add doc
#
test.z.matrix <- function () {
  fmt <- paste( "\n",
                " Test: %s\n",
                " Params: \n",
                "    campaigns: \n",
                "       %s",
                "    survey: \n",
                "       %s",
                " Expected: %s\n",
                " Actual:   %s\n" )
  verbose <- ( getOption( "RUnit" )$verbose && TRUE )

  if ( verbose ) { cat( "\n" ) }

  # === [ Bad Parameter Tests ] ===========================================

  # Test: obs and/or campaign = NA, NULL, not data frames,
  #       empty data frames
  # TBD

  # Test: invalid campaign$is.SIA values such as:
  #       NA, NULL, character, < 0, > 1, other???
  # TBD

  # Test: invalid campaign$date values such as:
  #       NA, NULL, character, non-date numeric,
  #       date with invalid year, month or day
  #       date with missing year, month or day ???
  #       date in the future, date far in past ???
  # TBD
  # Note: character dates do no work, ie as.numeric returns "1"
  

  # Test: invalid obs$date values such as:
  #       same as invalid campaign$date (see above)
  # TBD

  # Test: invalid obs$age values such as:
  #       NA, NULL, character, -Inf, +Inf, < 0,
  #      > (120 * 12), other???
  # TBD

  # Test: invalid campaign age.low and age.high values
  #       (either and both) such as:
  #       NA, NULL, character, -Inf, +Inf (acceptable?),
  #       < 0, low = high (acceptable?), low > high,
  #       other???
  # TBD
  if ( TRUE ) {
    name <- "invalid campaign age range, low > high"
    date <- as.Date( "2008-06-01" ); min <- 15;  max <- 20
    camps  <- data.frame( date=date, is.SIA=1, age.low=max, age.high=min )
    survey <- data.frame( date=( date + MONTH ),
                          age=c( min, min+1, min+2, max, max+1, max+2 ) )
    # note: ages at campaign = min/max -1/0/+1, ie subtract 1
    #       for month added to survey date
    if ( verbose ) { cat( sprintf("\t... Running %s ...\n", name) ) }
    expect <- matrix( data = 0, nrow = nrow(survey), ncol = nrow(camps) )
    actual <- z.matrix( survey, camps )
    msg <- sprintf( fmt, name, toStr(camps,"\t"), toStr(survey,"\t"),
                               toStr(expect), toStr(actual) )
    checkEquals( actual, expect, msg )
    # FIX: Change 'checkEquals' to 'checkException'
    #      The z.matrix function should raise exception if low > high
  }

  # Test: [Special Case] invalid age at campaign < 0 (calculated)
  #       with age.low < 0 and is.SIA = 1
  # FIX: Test disabled until z.matrix have validation logic
  if ( FALSE ) {
    name <- "invalid age at campaign < 0 (w/ age.low < 0 && is.SIA = TRUE)"
    date <- as.Date( "2008-06-01" ); min <- -10;  max <- 120
    camps  <- data.frame( date=date, is.SIA=1, age.low=min, age.high=max )
    survey <- data.frame( date=( date + 3 * MONTH ),
                          age=c( 1, 2, 3, 4 ) )
    if ( verbose ) { cat( sprintf("\t... Running %s ...\n", name) ) }
    expect <- matrix( data = c( 0, 0, 0, 1 ) )
    actual <- z.matrix( survey, camps )
    msg <- sprintf( fmt, name, toStr(camps,"\t"), toStr(survey,"\t"),
                               toStr(expect), toStr(actual) )
    checkEquals( actual, expect, msg )
    # FIX: Change 'checkEquals' to 'checkException'
    #      The z.matrix function should raise exception if low < 0
  }

#HERE
#  camps <- data.frame(
#        date = c( "2008-01-21", "2009-02-22", "2010-05-31", "2011-07-04" ),
#           N = c(        10000,        15000,        20000,        25000 ),
#           v = c(         1000,         1200,         1400,         1600 ),
#     age.low = c(            0,            0,            0,            0 ),
#    age.high = c(           36,          216,           36,          120 ),
#      is.SIA = c(            1,            1,            1,            1 )
#  )
               

  # === [ Good Parameter Tests ] ==========================================

  # Test: all non-SIA campaigns
  # TBD
  
  # Test: SIA campaigns w/ survey date <= campaign date
  # TBD
  
  # Test: SIA campaigns w/ survey date > campaign date
  #       && age at campaign < 0
  # TBD

  # Test: SIA campaigns w/ survey date > campaign date
  #       && age at campaign < campaign age low
  # TBD

  # Test: SIA campaigns w/ survey date > campaign date
  #       && age at campaign > campaign age high
  # TBD

  # Test: SIA campaigns w/ survey date > campaign date
  #       && age at campaign = campaign age low
  # TBD

  # Test: SIA campaigns w/ survey date > campaign date
  #       && age at campaign = campaign age high
  # TBD

  # Test: SIA campaigns w/ survey date > campaign date
  #       && age at campaign > campaign age low
  #       && age at campaign < campaign age high
  # TBD

  if ( verbose ) { cat( " ... " ) }

}  # end of function test.z.matrix()



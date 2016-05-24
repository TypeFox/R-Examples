cat("************** ECB sdmx   ******************************\n")

require("TSsdmx")

# RJSDMX::sdmxHelp()  # can be useful for finding series identifiers, etc

ecb <- TSconnect("sdmx", dbname="ECB")

#### annual####
z <- TSget('EXR.A.USD.EUR.SP00.A', ecb)
z <- TSget('EXR.A.USD.EUR.SP00.A',start = 2000, end = 2012, ecb)

if(1 != frequency(z)) stop("ECB monthly frequency error.")
if(2000 != start(z))  stop("ECB monthly start error.")


#### monthly ####
z <- TSget('EXR.M.USD.EUR.SP00.A', ecb)
if(! all(c(1999,1) == start(z))) stop("ECB monthly test 1 start error.")

z <- TSget('EXR.M.USD.EUR.SP00.A', start="2008-05", end="2014-07", ecb)
if(! all(c(2008,5) == start(z))) stop("ECB monthly test 2 start error.")
if(! all(c(2014,7) ==   end(z))) stop("ECB monthly test 2 end error.")

z <- TSget('EXR.M.USD.EUR.SP00.A', start=c(2008,5), end=c(2014,7), ecb)
if(! all(c(2008,5) == start(z))) stop("ECB monthly test 3 start error.")
if(! all(c(2014,7) ==   end(z))) stop("ECB monthly test 3 end error.")

z <- TSget('EXR.M.USD.EUR.SP00.A', ecb) 
z <- TSget('EXR.M.USD.EUR.SP00.A', start = "", end = "", ecb) 

#'should give error message  does not exist  on ECB
z <-   try(  TSget('122.ICP.M.U2.N.000000.4.ANR', ecb), silent=TRUE )
if (! grepl('122.ICP.M.U2.N.000000.4.ANR', attr(z,"condition"))) stop(
        'does not exist error not caught properly.')

require("tfplot")

options(TSconnection=ecb)

z <- TSget("ICP.M.U2.N.000000.4.ANR")# annual rates 
if(! all(c(1991,1) == start(z))) stop("ECB monthly test 4 start error.")
##   BUG  ?? tfplot(z)

z <- TSget("ICP.M.U2.N.000000.4.INX")# index 
z <- TSget("BSI.M.U2.Y.U.A21.A.4.U2.2250.Z01.E")

tfplot(z)
plot(z)


#### quarterly data ####

z <- TSget('EXR.Q.USD.EUR.SP00.A')

if(! all(c(1999,1) == start(z))) stop("ECB quarterly  test 1 start error.")

z <- TSget('EXR.Q.USD.EUR.SP00.A', start="2008-Q2", end="2014-Q3")
if(! all(c(2008,2) == start(z))) stop("ECB quarterly  test 2 start error.")
if(! all(c(2014,3) ==   end(z))) stop("ECB quarterly  test 2 end error.")

z <- TSget('EXR.Q.USD.EUR.SP00.A', start=c(2008,2), end=c(2014,3))
if(! all(c(2008,2) == start(z))) stop("ECB quarterly  test 3 start error.")
if(! all(c(2014,3) ==   end(z))) stop("ECB quarterly  test 3 end error.")

# BSI balance sheet indicators
#   FREQ Q
#   U2 Euro area, changing composition.
#   not adjusted 
#   sector T 
#   item  A20 loans   (no A21 cedit for cunsumption)
#   A all maturities
#   data type 1 outstanding
#   count area  U2 Euro area, changing composition.
#   count sector2250 household and non-profit...
#   currency  Z01 all currencies
#   suffix  E euro   (no B average )


z <- TSget("BSI.Q.U2.N.V.*.*.*.*.*.*.*")  #has A,F,R,V Feb 2015

# next is ok Feb 2015
if(! any("BSI.Q.U2.N.V.M30.X.1.U2.2250.Z01.E" %in% seriesNames(z) ))
    stop("available series has changed")  



# next is ok Feb 2015
skey <-c("BSI.Q.U2.N.V.M30.X.I.U2.2240.Z01.E",
         "BSI.Q.U2.N.V.M30.X.I.U2.2250.Z01.E",
         "BSI.Q.U2.N.V.M30.X.1.U2.2250.Z01.E")


z <- TSget(skey[1], ecb)   
  
tfplot(z)

z <- TSget(skey, ecb) 
tfplot(z)
  plot(z)


#### weeky data  ####

# Frequency W. Does not allow the creation of a strictly 
# fetching but then failing translating date with error
#   character string is not in a standard unambiguous format
#
# z <- TSget("ILM.W.U2.C.A010.Z5.Z0Z", ecb)
#
# Frequency W. Does not allow the creation of a strictly 
# fetching but then failing
# z <- TSget("ILM.W.U2.C.A010.Z5.Z0Z", ecb) # zoo BUG warning re underlying regularity


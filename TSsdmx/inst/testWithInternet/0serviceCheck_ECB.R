require("RJSDMX")



############################ "ECB" ############################
  getFlows('ECB')

  getFlows('ECB','*EXR*')

  getCodes('ECB', 'EXR', 'FREQ')
  names(getDimensions('ECB','EXR')) # I think this is also in correct order
  getDimensions('ECB','EXR')

  getDSDIdentifier('ECB','EXR')

#### annual ####
  
  z <- getSDMX("ECB", 'EXR.A.USD.EUR.SP00.A')
  if(1999 != start(z[[1]])) stop("ECB annual retrieval error.")

  z <- getSDMX("ECB", 'EXR.A.USD.EUR.SP00.A', start = "2001", end = "2012")
  if(2001 != start(z[[1]])) stop("start test for ECB annual data failed.")
  if(2012 != end(z[[1]]))   stop(  "end test for ECB annual data failed.")

  
#### monthly ####

  z <- getSDMX("ECB", 'EXR.M.USD.EUR.SP00.A')

  if("Jan 1999" != start(z[[1]])) stop("ECB monthly retrieval error (start check).")
  if(12 != frequency(z[[1]])) stop("ECB monthly retrieval error (frequency check).")

  z <- getSDMX("ECB", 'EXR.M.USD.EUR.SP00.A', start="2008", end="2013")
  if("Jan 2008" != start(z[[1]])) stop("ECB monthly start specification 1 failure.")

  z <- getSDMX("ECB", 'EXR.M.USD.EUR.SP00.A', start="2008-05", end="2014-07")[[1]]
  if("May 2008" != start(z)) stop("ECB monthly start specification 2 failure.")
  if("Jul 2014" != end(z))   stop("ECB monthly  end  specification 2 failure.")

  z <- getSDMX("ECB", 'EXR.M.USD.EUR.SP00.A', start="2008-Q1", end="2014-Q2")[[1]]
  if("Jan 2008" != start(z)) stop("ECB monthly start specification 3 failure.")
  if("Jun 2014" != end(z))   stop("ECB monthly  end  specification 3 failure.")

  z <- getSDMX("ECB", 'EXR.M.USD.EUR.SP00.A',
                       start="2008-01-01", end="2014-01-31")[[1]]
  if("Jan 2008" != start(z)) stop("ECB monthly start specification 4 failure.")
  if("Jan 2014" != end(z))   stop("ECB monthly  end  specification 4 failure.")


#### quarterly ####

  z <- getSDMX("ECB", 'EXR.Q.USD.EUR.SP00.A')

  if("1999 Q1" != start(z[[1]]))
                   stop("ECB quarterly retrieval error (start check).")
  if(4 != frequency(z[[1]]))
                   stop("ECB quarterly retrieval error (frequency check).")

  z <- getSDMX("ECB", 'EXR.Q.USD.EUR.SP00.A', start="2008-Q2", end="2014-Q3")[[1]]
  if("2008 Q2" != start(z)) stop("ECB quarterly start specification 1 failure.")
  if("2014 Q3" != end(z))   stop("ECB quarterly  end  specification 1 failure.")

#### weeky data  ####

# "Frequency W. 

  z <- getSDMX("ECB", "ILM.W.U2.C.A010.Z5.Z0Z")
  
  if(start(z[[1]]) != "1998-W53") stop("ECB weeky retrieval changed start date.")


  # this would make sense but does not work
  #z <- getSDMX("ECB", "ILM.W.U2.C.A010.Z5.Z0Z",start="2008-W1",end="2013-W1")[[1]] 

  # I'm not sure if these are correct
  #z <- getSDMX("ECB", "ILM.W.U2.C.A010.Z5.Z0Z",start="2008-Q1",end="2013-Q4")[[1]] 

start(z)  # "2008-W02"
end(z)    # "2014-W01"

  
  #if(start(z) != "1998-W53") stop("ECB weeky retrieval changed start date.")

  #z <- getSDMX("ECB", "ILM.W.U2.C.A010.Z5.Z0Z", start="2008", end="2013")[[1]] 
  #if(start(z) != "1998-W53") stop("ECB weeky retrieval changed start date.")

# Dates in above are in form 1998-W53. These might be converted to dates with
# require(ISOweek)
#    # assume Wednesday, weekday=3, but there may be more information
#    # available in the SDMX
#    dt <- as.Date(ISOweek::ISOweek2date(paste(times,"-3", sep="")))
#    tmp_ts <- zoo(values, order.by = dt)
#    }

#### daily data  ####

# select years
  z <- getSDMX('ECB', 'EXR.D.USD.EUR.SP00.A', '2000', '2001')[[1]]
  if("2000-01-03" != start(z)) stop("ECB daily start specification 1 failure.")
  if("2001-12-31" != end(z))   stop("ECB daily  end  specification 1 failure.")

 frequency(z) # check this

# select months
  z <- getSDMX('ECB', 'EXR.D.USD.EUR.SP00.A', '2000-01', '2000-12')[[1]]
  if("2000-01-03" != start(z)) stop("ECB daily start specification 2 failure.")
  if("2000-12-29" != end(z))   stop("ECB daily  end  specification 2 failure.")

# select quarters
  z <- getSDMX('ECB', 'EXR.D.USD.EUR.SP00.A', '2000-Q1', '2000-Q2')[[1]]
  if("2000-01-03" != start(z)) stop("ECB daily start specification 3 failure.")
  if("2000-06-30" != end(z))   stop("ECB daily  end  specification 3 failure.")

# select days
  z <- getSDMX('ECB', 'EXR.D.USD.EUR.SP00.A', '2000-01-01', '2000-01-31')[[1]]
  if("2000-01-03" != start(z)) stop("ECB daily start specification 4 failure.")
  if("2000-01-31" != end(z))   stop("ECB daily  end  specification 4 failure.")



## These get mixed monthly and annual frequency, but standard R time series
##  representations do not handle that very well.
##z1 <- getSDMX('ECB', 'EXR.A|M.USD.EUR.SP00.A')
##z2 <- getSDMX('ECB', 'EXR.A+M.USD.EUR.SP00.A')
## get mixed all available frequencies
##z <- getSDMX('ECB', 'EXR.*.USD.EUR.SP00.A')


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


#z <- getSDMX('ECB', "BSI.Q.U2.N.T.A21.A.1.U2.2250.Z01.B") #not found
#z <- getSDMX('ECB', "BSI.Q.U2.N.T.A20.A.1.U2.2250.Z01.E") 
#z <- getSDMX('ECB', "BSI.Q.U2.N.T.*.*.*.*.*.*.*")  #not found

z <- getSDMX('ECB', "BSI.Q.U2.N.V.*.*.*.*.*.*.*")  #has A,F,R,V Feb 2015

#z <- getSDMX('ECB', "BSI.Q.U2.N.*.*.*.*.*.*.*.*") works 
#length(names(z) #927

nm <- names(z)
length(nm) # for A=530  V=92
nm

sum(grepl('A21',nm))  # 0
sum(grepl('.B', names(z)))  # 0
sum(grepl('2250', names(z))) # 22, Feb 2015

# next is ok Feb 2015
if(! any("BSI.Q.U2.N.V.M30.X.1.U2.2250.Z01.E" %in% nm ))
    stop("available series has changed")  



require("TSsdmx")

# RJSDMX::sdmxHelp()  # can be useful for finding series identifiers, etc

oecd <- TSconnect("sdmx", dbname="OECD")

#  environmental indicators. single point, not series?
#  tts <- getSDMX('OECD', 'ENV_KEI.*.NOX_GDP')   # date problems
#  names(tts) # 34 countries
# z <- TSget('ENV_KEI.CAN.NOX_GDP', oecd)
# z <- TSget('ENV_KEI.*.NOX_GDP', oecd)


#### monthly data ####

z <- TSget('G20_PRICES.CAN.CPALTT01.IXOB.M',oecd) 
if(! all(c(1949,1) == start(z))) stop('monthly test 1 start date is changed.')
if(12 != frequency(z))           stop('monthly test 1 frequency error.')

tframe::seriesNames(z)

z <- TSget('G20_PRICES.CAN.CPALTT01.IXOB.M', start=c(1990,6), end =c(2012,12), oecd)  
if(! all(c(1990,6)  == start(z))) stop('monthly test 2 start date error.')
if(! all(c(2012,12) ==   end(z))) stop('monthly test 2 end date error.')


#### quarterly data ####

# quarterly national accounts
#CARSA: national currency, nominal, SAAR (level)

if (FALSE) { # failing with  "Premature end of file."  Dec 11, 2015

z <- TSget('QNA.CAN.PPPGDP.CARSA.Q', oecd)
if(! all(c(1960,1) == start(z))) stop('quarterly test 1 start date is changed.')
if(4 != frequency(z)) stop('quarterly test 1 frequency error.')

z <- TSget('QNA.CAN.PPPGDP.CARSA.Q', start=c(1990,2), end =c(2012,4), oecd)
if(! all(c(1990,2) == start(z))) stop('quarterly test 2 start date error.')
if(! all(c(2012,4) ==   end(z))) stop('quarterly test 2 end date error.')

#   test "+" and "|" in query and test setting names
z <- TSget('QNA.CAN+USA|MEX.PPPGDP.CARSA.Q', 
         names=c("Canada", "United States", "Mexico"),   oecd)

# SDMX + and | queries do not determine the return order, TSget fixes by
# reordering data. (Above query  was not returned in order in Dec 2014.)
# The PPPGDP numbers are all relative to the US so USA numbers are 1.0 and next
#  test checks that as a confirmation that re-order was done. 
# This was BUG #22 which was closed with work around in RJSDMX by using ; to
#  separate queries and maintain order.

if(max(abs(1 - z[,2])) > 1e-16)
          stop('quarterly test reorder series to apply names not working.')

if(! all(c("Canada", "United States", "Mexico") == tframe::seriesNames(z)))
          stop('quarterly test setting series names not working.')

if(! all(c(1955,1) == start(z))) 
          stop('quarterly mulivariate test start date is changed.')

if(4 != frequency(z)) 
          stop('quarterly mulivariate test frequency error.')

# tfplot::tfplot(z, graphs.per.page=3)
# tfplot::tfOnePlot(z, start=c(1990,1))


} # end if (FALSE)

# Annual only ??
z <- TSget('BSI.NAT.EQU.TOT.DIR.CAN', oecd)  
if(! all(c(2009,1) == start(z))) stop('annual test 0 start date is changed.')
if(1 != frequency(z)) stop('annual test 0 frequency error.')

# tts <- getSDMX('OECD', 'BSI.NAT.*.*.*.CAN')
# names(tts)


#####  annual #####


#  USING sdmxHelp()
#>OECD
# 	G20_PRICES > 	>LOCATION:CAN
#			>SUBJECT : CP   (CPI)
#			>MEASURE : IXOB (INDEX)
#			>FREQUENCY: M

z <- TSget('G20_PRICES.CAN.CPALTT01.IXOB.A',oecd)  
if(! all(c(1949,1) == start(z))) stop('annual test 1 start date is changed.')
if(1 != frequency(z)) stop('annual test 1 frequency error.')


#>OECD  Household ... assets and liabilities ..
# 	7HA_A_Q > 	>LOCATION:CAN
#			>TRANSACTION : AF411LI (cons. cr No 7HAL2=liabilities)
#			>ACTIVITY : ST (stocks)
#			>MEASURE : C  (nominal CDN)
#			>FREQUENCY:  A only  *** Q does not exist ***


#z <- TSget('7HA_A_Q.CAN.*.*.*.*', oecd)
#tframe::seriesNames(z)
#tframe::nseries(z) # 44

z <- TSget('7HA_A_Q.CAN.AF411LI.ST.C.A', start=c(1999,1), end =c(2009,1), oecd)
if(! all(c(1999,1) == start(z))) stop('annual test 2 start date error.')
if(! all(c(2009,1) ==   end(z))) stop('annual test 2 end date error.')

z <- TSget('7HA_A_Q.CAN.*.ST.C.A', oecd)
tframe::seriesNames(z)

if(! all(c(1995,1) ==  start(z))) stop('annual mulivariate test start date is changed.')
if(1 != frequency(z)) stop('annual mulivariate test frequency error.')

#silent=TRUE only works properly if levels are set OFF in 
#   the Sdmx configuration file
#z <- try(TSget('G20_PRICES.CAB.CP.IXOB.M',oecd), silent=TRUE)
#z <-     TSget('G20_PRICES.CAB.CP.IXOB.M',oecd)


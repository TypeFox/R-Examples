############################ "ILO" ############################
# http://www.ilo.org/
# The server process which provides information to sdmxHelp() is having problems
#  (Actually, not just the sdmxHelp. The server is not available.)

# The ILO endpoint does not respond to the getFlows call (probably too many flows) and 
# this makes it impossible to use the sdmxHelper for it.

# There is a nice API doc on the ILOSTAT site that helps building queries:

#http://www.ilo.org/ilostat/content/conn/ILOSTATContentServer/path/Contribution%20Folders/statistics/web_pages/static_pages/technical_page/ilostat_appl/SDMX_User_Guide.pdf

# A definitive fix is coming, but in the meantime the getflows already works.

require("RJSDMX")

# z <- getFlows('ILO')  # very slow but works Feb 13, 2015. Previously failed (see above)
# length(z)  #[1] 710

########### annual ########### 

tts <- getTimeSeries("ILO", "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB/YI.MEX.A.463.EMP_TEMP_NB.SEX_F.AGE_10YRBANDS_TOTAL")

if ("1988" != start(tts[[1]]))  stop("ILO test 1 start date changed.")


tts <- getSDMX("ILO",
   "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB/YI.MEX.A.463.EMP_TEMP_NB.SEX_F.AGE_10YRBANDS_TOTAL",
   start="1995-01-01", end="2012-12-31")
#   start="1995", end="2012")

if ("1995" != start(tts[[1]]))  stop("ILO test 1 start date error.")

if ("2012" !=   end(tts[[1]]))  stop("ILO test 1 end date error.") 

#17 provider BUG
# ILOSTAT support said, Feb 11, 2015, the ILOSTAT web services only support 
# data format YYYY-MM-DD. Thus, in order to have the observation for 
# year 2012 you'll have to change the call  from start="1995", end="2012")
# to  use start="1995-01-01", end="2012-12-31")

# works but slow
# tts <- getTimeSeries("ILO", "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB/YI.MEX+ESP.....")
# 
# nm <- names(tts)
# nm[grepl('SEX_T.AGE_AGGREGATE_TOTAL', nm )] 
 
 
tts <- getTimeSeries("ILO", "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB/YI.DEU+FRA+GBR+ITA.A..EMP_TEMP_NB.SEX_T.AGE_AGGREGATE_TOTAL")

length(tts) # 10

if ("2001" != start(tts[[1]]))  stop("ILO test 3 start date changed.")


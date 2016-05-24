# http://www.ilo.org/
# http://www.ilo.org/global/statistics-and-databases/lang--en/index.htm

# There is a nice API doc on the ILOSTAT site that helps building queries:

#http://www.ilo.org/ilostat/content/conn/ILOSTATContentServer/path/Contribution%20Folders/statistics/web_pages/static_pages/technical_page/ilostat_appl/SDMX_User_Guide.pdf

# http://laborsta.ilo.org/
# > browse by country
# > Canada
# >Unemployment
# >by sex and age


require("TSsdmx")

ilo <- TSconnect("sdmx", dbname="ILO")

z <- TSget("DF_YI_ALL_EMP_TEMP_SEX_AGE_NB/YI.MEX.A.463.EMP_TEMP_NB.SEX_F.AGE_10YRBANDS_TOTAL",
        ilo)

if (! all(c(1988,1) == start(z)))  stop("ILO test 1 start date changed.")


z <- TSget("DF_YI_ALL_EMP_TEMP_SEX_AGE_NB/YI.MEX.A.463.EMP_TEMP_NB.SEX_F.AGE_10YRBANDS_TOTAL",
        start=c(1995,1), end=c(2012,1), ilo)

if (! all(c(1995,1) == start(z)))  stop("ILO test 2 start date error.")


#if (!all(c(2012,1) == end(z))) stop("ILO test 2 end date error.") #17 provider BUG must be done as next (see notes in 0ServiceCheck_ILO.R)

z <- TSget("DF_YI_ALL_EMP_TEMP_SEX_AGE_NB/YI.MEX.A.463.EMP_TEMP_NB.SEX_F.AGE_10YRBANDS_TOTAL",
        start="1995-01-01", end="2012-12-31", ilo)
if (!all(c(2012,1) == end(z))) stop("ILO test 2 end date error.")  


z <- TSget( "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB/YI.DEU+FRA+GBR+ITA.A..EMP_TEMP_NB.SEX_T.AGE_AGGREGATE_TOTAL", ilo)

# was 10 circa spring 2015
if (11 !=  tframe::nseries(z))
    stop("ILO test 3 number of series changed (again).")
if (! all(c(1969,1) == start(z)))  stop("ILO test 3 start date error.")
 
tframe::seriesNames(z)
# [1] "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB.YI.DEU.A.1067.EMP_TEMP_NB.SEX_T.AGE_AGGREGATE_TOTAL"
# [2] "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB.YI.DEU.A.1068.EMP_TEMP_NB.SEX_T.AGE_AGGREGATE_TOTAL"
# [3] "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB.YI.DEU.A.2242.EMP_TEMP_NB.SEX_T.AGE_AGGREGATE_TOTAL"
# [4] "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB.YI.FRA.A.47.EMP_TEMP_NB.SEX_T.AGE_AGGREGATE_TOTAL"  
# [5] "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB.YI.FRA.A.1142.EMP_TEMP_NB.SEX_T.AGE_AGGREGATE_TOTAL"
# [6] "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB.YI.FRA.A.2260.EMP_TEMP_NB.SEX_T.AGE_AGGREGATE_TOTAL"

# new as of checking in Sept 2015 
# [7] "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB.YI.GBR.A.666.EMP_TEMP_NB.SEX_T.AGE_AGGREGATE_TOTAL" 

# [7] "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB.YI.GBR.A.1155.EMP_TEMP_NB.SEX_T.AGE_AGGREGATE_TOTAL"
# [8] "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB.YI.GBR.A.2247.EMP_TEMP_NB.SEX_T.AGE_AGGREGATE_TOTAL"
# [9] "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB.YI.ITA.A.325.EMP_TEMP_NB.SEX_T.AGE_AGGREGATE_TOTAL" 
#[10] "DF_YI_ALL_EMP_TEMP_SEX_AGE_NB.YI.ITA.A.2238.EMP_TEMP_NB.SEX_T.AGE_AGGREGATE_TOTAL"


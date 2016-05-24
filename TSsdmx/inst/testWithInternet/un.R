require("TSsdmx")

# note REST yet

# un <- TSconnect("sdmx", dbname == "UN")

# Status:  Not nearly working
#    Preliminary investigation notes. Web site may require account and
#    authentication. 
 #http://unstats.un.org/unsd/tradekb/Knowledgebase/Comtrade-SDMX-Web-Services-and-Data-Exchange

#You can get the SDMX format by clicking the SDMX download in any of query result #pages (for an example at #http://comtrade.un.org/db/dqBasicQueryResults.aspx?cc=62&px=BE&r=251&p=705&rg=2&so=6)

#The SDMX key families and simple power point presentation about Comtrade SDMX #are attached to this article: UN Comtrade SDMX.
#http://unstats.un.org/unsd/tradekb/Knowledgebase/How-to-use-Batch-Processing-feature-in-UN-Comtrade
 
#from #http://unstats.un.org/unsd/tradekb/Knowledgebase/Comtrade-SDMX-Web-Services-and-Data-Exchange?Keywords=SDMX

#URL Example

#Country code : 381 = Italy
#Classification : H1 = HS1996
#Years : 2002 and 2003 
#http://comtrade.un.org/ws/getSdmxV1.aspx?px=H1&r=381&y=2003,2002&cc=TOTAL&p=0&comp=false
#Country code : 36 = Australia
#Classification : S2 = SITC Rev.2
#Year : 2004
#http://comtrade.un.org/ws/getSdmxV1.aspx?px=s2&r=36&y=2004&cc=TOTAL&p=0&comp=false

#(requires permission. access denied)

############################ UN ############################

##R example using json at
##http://comtrade.un.org/data/Doc/api/ex/r


# no SDMX yet but coming
# http://comtrade.un.org/data/doc/api/
# http://comtrade.un.org/data/doc/api/#Future

# UN Comtrade data request takes the following form:
# http://comtrade.un.org/api/get?parameters
# API call: 
# http://comtrade.un.org/api/get?max=50000&type=C&freq=A&px=HS&ps=2013&r=826&p=0&rg=all&cc=AG2&fmt=json

# http://comtrade.un.org/api/get?max=50000&type=C&freq=A&px=HS&ps=2013&r=826&p=0&rg=all&cc=AG2&fmt=sdmx

#Old
#http://unstats.un.org/unsd/tradekb/Knowledgebase/Comtrade-SDMX-Web-Services-and-Data-Exchange
#http://unstats.un.org/unsd/tradekb/Knowledgebase/Comtrade-SDMX-Web-Services-and-Data-Exchange?Keywords=SDMX

#
#  Copyright 2014 jSonar Inc
#  All Rights Reserved.
#
#  Licensed under the GNU AFFERO GENERAL PUBLIC LICENSE version 3
#  See http://www.r-project.org/Licenses/AGPL-3
#

exampleURL <- 'https://example.com/Gateway?host=localhost&port=27017&db=test'
exampleData <- {}

exampleData[paste(exampleURL, 'output=json&type=find&name=delayed_flights&col=ExampleFlights&limit=5', sep='&')] <-
'{"output": [{ "_id" : { "$oid" : "53c6b753754c39f2a9ec5520"} , "FlightNum" : 3231 , "TaxiIn" : 5 , "SecurityDelay" : 0 , "DepTime" : 754 , "DepDelay" : 19 , "WeatherDelay" : 0 , "CRSArrTime" : 1000 , "DayofMonth" : 3 , "DayOfWeek" : 4 , "TaxiOut" : 10 , "Dest" : { "city" : "Tampa" , "country" : "USA" , "iata" : "TPA" , "airport" : "Tampa International " , "long" : -82.53325 , "state" : "FL" , "lat" : 27.97547222} , "CRSElapsedTime" : 145 , "ArrDelay" : 2 , "AirTime" : 113 , "CarrierDelay" : 0 , "CRSDepTime" : 735 , "Diverted" : 0 , "Distance" : 810 , "UniqueCarrier" : { "Code" : "WN" , "Description" : "Southwest Airlines Co."} , "NASDelay" : 0 , "Cancelled" : 0 , "TailNum" : { "status" : "Valid" , "issue_date" : "08/07/2000" , "aircraft_type" : "Fixed Wing Multi-Engine" , "year" : 2000 , "tailnum" : "N772SW" , "model" : "737-7H4" , "type" : "Corporation" , "engine_type" : "Turbo-Fan" , "manufacturer" : "BOEING"} , "Origin" : { "city" : "Chantilly" , "country" : "USA" , "iata" : "IAD" , "airport" : "Washington Dulles International" , "long" : -77.45580972 , "state" : "VA" , "lat" : 38.94453194} , "LateAircraftDelay" : 0 , "Month" : 1 , "ActualElapsedTime" : 128 , "Year" : 2008 , "ArrTime" : 1002 , "CancellationCode" : ""},{ "_id" : { "$oid" : "53c6b753754c39f2a9ec5521"} , "FlightNum" : 448 , "TaxiIn" : 3 , "SecurityDelay" : 0 , "DepTime" : 628 , "DepDelay" : 8 , "WeatherDelay" : 0 , "CRSArrTime" : 750 , "DayofMonth" : 3 , "DayOfWeek" : 4 , "TaxiOut" : 17 , "Dest" : { "city" : "Baltimore" , "country" : "USA" , "iata" : "BWI" , "airport" : "Baltimore-Washington International" , "long" : -76.66819833 , "state" : "MD" , "lat" : 39.17540167} , "CRSElapsedTime" : 90 , "ArrDelay" : 14 , "AirTime" : 76 , "CarrierDelay" : 0 , "CRSDepTime" : 620 , "Diverted" : 0 , "Distance" : 515 , "UniqueCarrier" : { "Code" : "WN" , "Description" : "Southwest Airlines Co."} , "NASDelay" : 0 , "Cancelled" : 0 , "TailNum" : { "status" : "Valid" , "issue_date" : "12/20/2002" , "aircraft_type" : "Fixed Wing Multi-Engine" , "year" : 2002 , "tailnum" : "N428WN" , "model" : "737-7H4" , "type" : "Corporation" , "engine_type" : "Turbo-Fan" , "manufacturer" : "BOEING"} , "Origin" : { "city" : "Indianapolis" , "country" : "USA" , "iata" : "IND" , "airport" : "Indianapolis International" , "long" : -86.29438417 , "state" : "IN" , "lat" : 39.71732917} , "LateAircraftDelay" : 0 , "Month" : 1 , "ActualElapsedTime" : 96 , "Year" : 2008 , "ArrTime" : 804 , "CancellationCode" : ""},{ "_id" : { "$oid" : "53c6b753754c39f2a9ec5524"} , "FlightNum" : 378 , "TaxiIn" : 4 , "SecurityDelay" : 0 , "DepTime" : 1940 , "DepDelay" : 25 , "WeatherDelay" : 0 , "CRSArrTime" : 2110 , "DayofMonth" : 3 , "DayOfWeek" : 4 , "TaxiOut" : 10 , "Dest" : { "city" : "Jacksonville" , "country" : "USA" , "iata" : "JAX" , "airport" : "Jacksonville International" , "long" : -81.68786111 , "state" : "FL" , "lat" : 30.49405556} , "CRSElapsedTime" : 115 , "ArrDelay" : 11 , "AirTime" : 87 , "CarrierDelay" : 0 , "CRSDepTime" : 1915 , "Diverted" : 0 , "Distance" : 688 , "UniqueCarrier" : { "Code" : "WN" , "Description" : "Southwest Airlines Co."} , "NASDelay" : 0 , "Cancelled" : 0 , "TailNum" : { "status" : "Valid" , "issue_date" : "03/31/1999" , "aircraft_type" : "Fixed Wing Multi-Engine" , "year" : 1999 , "tailnum" : "N726SW" , "model" : "737-7H4" , "type" : "Corporation" , "engine_type" : "Turbo-Fan" , "manufacturer" : "BOEING"} , "Origin" : { "city" : "Indianapolis" , "country" : "USA" , "iata" : "IND" , "airport" : "Indianapolis International" , "long" : -86.29438417 , "state" : "IN" , "lat" : 39.71732917} , "LateAircraftDelay" : 0 , "Month" : 1 , "ActualElapsedTime" : 101 , "Year" : 2008 , "ArrTime" : 2121 , "CancellationCode" : ""},{ "_id" : { "$oid" : "53c6b753754c39f2a9ec5523"} , "FlightNum" : 3920 , "TaxiIn" : 3 , "SecurityDelay" : 0 , "DepTime" : 1829 , "DepDelay" : 34 , "WeatherDelay" : 0 , "CRSArrTime" : 1925 , "DayofMonth" : 3 , "DayOfWeek" : 4 , "TaxiOut" : 10 , "Dest" : { "city" : "Baltimore" , "country" : "USA" , "iata" : "BWI" , "airport" : "Baltimore-Washington International" , "long" : -76.66819833 , "state" : "MD" , "lat" : 39.17540167} , "CRSElapsedTime" : 90 , "ArrDelay" : 34 , "AirTime" : 77 , "CarrierDelay" : 2 , "CRSDepTime" : 1755 , "Diverted" : 0 , "Distance" : 515 , "UniqueCarrier" : { "Code" : "WN" , "Description" : "Southwest Airlines Co."} , "NASDelay" : 0 , "Cancelled" : 0 , "TailNum" : { "status" : "Valid" , "issue_date" : "07/23/2004" , "aircraft_type" : "Fixed Wing Multi-Engine" , "year" : 2004 , "tailnum" : "N464WN" , "model" : "737-7H4" , "type" : "Corporation" , "engine_type" : "Turbo-Fan" , "manufacturer" : "BOEING"} , "Origin" : { "city" : "Indianapolis" , "country" : "USA" , "iata" : "IND" , "airport" : "Indianapolis International" , "long" : -86.29438417 , "state" : "IN" , "lat" : 39.71732917} , "LateAircraftDelay" : 32 , "Month" : 1 , "ActualElapsedTime" : 90 , "Year" : 2008 , "ArrTime" : 1959 , "CancellationCode" : ""},{ "_id" : { "$oid" : "53c6b753754c39f2a9ec5525"} , "FlightNum" : 509 , "TaxiIn" : 3 , "SecurityDelay" : 0 , "DepTime" : 1937 , "DepDelay" : 67 , "WeatherDelay" : 0 , "CRSArrTime" : 1940 , "DayofMonth" : 3 , "DayOfWeek" : 4 , "TaxiOut" : 7 , "Dest" : { "city" : "Las Vegas" , "country" : "USA" , "iata" : "LAS" , "airport" : "McCarran International" , "long" : -115.1523333 , "state" : "NV" , "lat" : 36.08036111} , "CRSElapsedTime" : 250 , "ArrDelay" : 57 , "AirTime" : 230 , "CarrierDelay" : 10 , "CRSDepTime" : 1830 , "Diverted" : 0 , "Distance" : 1591 , "UniqueCarrier" : { "Code" : "WN" , "Description" : "Southwest Airlines Co."} , "NASDelay" : 0 , "Cancelled" : 0 , "TailNum" : { "status" : "Valid" , "issue_date" : "05/24/2000" , "aircraft_type" : "Fixed Wing Multi-Engine" , "year" : 2000 , "tailnum" : "N763SW" , "model" : "737-7H4" , "type" : "Corporation" , "engine_type" : "Turbo-Fan" , "manufacturer" : "BOEING"} , "Origin" : { "city" : "Indianapolis" , "country" : "USA" , "iata" : "IND" , "airport" : "Indianapolis International" , "long" : -86.29438417 , "state" : "IN" , "lat" : 39.71732917} , "LateAircraftDelay" : 47 , "Month" : 1 , "ActualElapsedTime" : 240 , "Year" : 2008 , "ArrTime" : 2037 , "CancellationCode" : ""},{ "_id" : { "$oid" : "53c6b753754c39f2a9ec5527"} , "FlightNum" : 11 , "TaxiIn" : 6 , "SecurityDelay" : 0 , "DepTime" : 617 , "DepDelay" : 2 , "WeatherDelay" : 0 , "CRSArrTime" : 650 , "DayofMonth" : 3 , "DayOfWeek" : 4 , "TaxiOut" : 19 , "Dest" : { "city" : "Kansas City" , "country" : "USA" , "iata" : "MCI" , "airport" : "Kansas City International" , "long" : -94.71390556 , "state" : "MO" , "lat" : 39.29760528} , "CRSElapsedTime" : 95 , "ArrDelay" : 2 , "AirTime" : 70 , "CarrierDelay" : 0 , "CRSDepTime" : 615 , "Diverted" : 0 , "Distance" : 451 , "UniqueCarrier" : { "Code" : "WN" , "Description" : "Southwest Airlines Co."} , "NASDelay" : 0 , "Cancelled" : 0 , "TailNum" : { "status" : "Valid" , "issue_date" : "03/11/1997" , "aircraft_type" : "Fixed Wing Multi-Engine" , "year" : 1985 , "tailnum" : "N689SW" , "model" : "737-3Q8" , "type" : "Corporation" , "engine_type" : "Turbo-Fan" , "manufacturer" : "BOEING"} , "Origin" : { "city" : "Indianapolis" , "country" : "USA" , "iata" : "IND" , "airport" : "Indianapolis International" , "long" : -86.29438417 , "state" : "IN" , "lat" : 39.71732917} , "LateAircraftDelay" : 0 , "Month" : 1 , "ActualElapsedTime" : 95 , "Year" : 2008 , "ArrTime" : 652 , "CancellationCode" : ""}]}'

exampleData[paste(exampleURL, 'output=csv&type=agg&name=delays_by_day&col=NYCFlights', sep='&')] <-
'_id.Month,_id.DayofMonth,_id.Origin,count,_sum_AirTime,_sum_ArrDelay,_sum_DepDelay,_avg_AirTime,_avg_ArrDelay,_avg_DepDelay,month,day,origin
1,1,EWR,820,126606,21202,20904,155.91871921182266,26.110837438423644,25.74384236453202,1,1,EWR
1,1,JFK,692,122424,12554,12970,180.56637168141592,18.51622418879056,19.073529411764707,1,1,JFK
1,1,LGA,586,74220,8380,6980,129.30313588850174,14.599303135888501,12.160278745644598,1,1,LGA
1,2,JFK,726,120500,21900,18710,166.89750692520775,30.332409972299168,25.914127423822716,1,2,JFK
1,2,EWR,904,121634,35432,42044,139.80919540229885,40.726436781609195,48.326436781609196,1,2,EWR
1,2,LGA,660,72390,6882,10516,112.40683229813665,10.686335403726709,16.32919254658385,1,2,LGA
1,3,LGA,684,68700,-5138,4830,102.53731343283582,-7.66865671641791,7.208955223880597,1,3,LGA
1,3,EWR,896,116956,18944,24924,133.81693363844394,21.675057208237988,28.45205479452055,1,3,EWR
1,3,JFK,716,115432,8874,15530,164.43304843304844,12.64102564102564,22.122507122507123,1,3,JFK
1,4,JFK,714,114286,-1348,6810,164.20402298850576,-1.9367816091954022,9.645892351274787,1,4,JFK
1,4,LGA,686,70974,-8654,316,106.56756756756756,-12.993993993993994,0.47305389221556887,1,4,LGA
1,4,EWR,894,117838,4450,11240,137.02093023255813,5.174418604651163,13.00925925925926,1,4,EWR
1,5,EWR,782,110434,4680,8022,142.31185567010309,6.030927835051546,10.284615384615385,1,5,EWR
1,5,LGA,480,55058,-3086,1852,118.65948275862068,-6.650862068965517,3.9742489270386265,1,5,LGA
1,5,JFK,686,111176,-1888,5630,165.93432835820894,-2.817910447761194,8.37797619047619,1,5,JFK
1,6,EWR,804,109298,-474,7682,137.30904522613065,-0.5954773869346733,9.650753768844222,1,6,EWR
1,6,JFK,680,113856,-1484,7244,168.92581602373886,-2.201780415430267,10.747774480712167,1,6,JFK
1,6,LGA,600,61252,-5604,2798,106.71080139372822,-9.763066202090592,4.857638888888889,1,6,LGA
1,7,EWR,840,112738,-2920,5736,135.17745803357315,-3.501199040767386,6.844868735083533,1,7,EWR
1,7,JFK,674,113730,-3232,5420,169.7462686567164,-4.823880597014925,8.08955223880597,1,7,JFK
1,7,LGA,732,75768,-390,2482,107.3201133144476,-0.5524079320113314,3.505649717514124,1,7,LGA
1,8,JFK,632,106798,-3212,2406,171.15064102564102,-5.147435897435898,3.8312101910828025,1,8,JFK
1,8,EWR,798,107222,1254,6086,135.38131313131314,1.5833333333333333,7.664987405541562,1,8,EWR'

exampleData[paste(exampleURL, 'output=csv&type=find&name=delayed_flights&col=WAFlights', sep='&')] <-
'_id,YEAR,MONTH,DAY_OF_MONTH,ORIGIN_AIRPORT_ID,ORIGIN_AIRPORT_SEQ_ID,ORIGIN_CITY_MARKET_ID,ORIGIN_CITY_NAME,ORIGIN_STATE_ABR,DEST_AIRPORT_ID,DEST_AIRPORT_SEQ_ID,DEST_CITY_MARKET_ID,DEST_CITY_NAME,DEST_STATE_ABR,ACTUAL_ELAPSED_TIME,WEATHER_DELAY
53750f8b2b860941a2ef5036,2014,1,4,11292,1129202,30325,"Denver, CO",CO,14747,1474703,30559,"Seattle, WA",WA,206.0,2.0
53750f8b2b860941a2ef507e,2014,1,4,11292,1129202,30325,"Denver, CO",CO,11884,1188402,31884,"Spokane, WA",WA,157.0,30.0
53750f8b2b860941a2ef50d4,2014,1,5,11292,1129202,30325,"Denver, CO",CO,14747,1474703,30559,"Seattle, WA",WA,196.0,3.0
53750f8b2b860941a2ef5667,2014,1,1,11292,1129202,30325,"Denver, CO",CO,11884,1188402,31884,"Spokane, WA",WA,182.0,23.0
53750f8b2b860941a2ef6e55,2014,1,4,11109,1110902,30189,"Colorado Springs, CO",CO,14747,1474703,30559,"Seattle, WA",WA,243.0,237.0
53750f8b2b860941a2ef6e62,2014,1,27,11109,1110902,30189,"Colorado Springs, CO",CO,14747,1474703,30559,"Seattle, WA",WA,208.0,6.0
53750f8b2b860941a2ef7f33,2014,1,19,11292,1129202,30325,"Denver, CO",CO,14252,1425202,34252,"Pasco/Kennewick/Richland, WA",WA,137.0,17.0
53750f8b2b860941a2ef8080,2014,1,30,11292,1129202,30325,"Denver, CO",CO,14252,1425202,34252,"Pasco/Kennewick/Richland, WA",WA,173.0,55.0
53750f8b2b860941a2ef81d5,2014,1,5,11292,1129202,30325,"Denver, CO",CO,14252,1425202,34252,"Pasco/Kennewick/Richland, WA",WA,182.0,55.0
53750f8b2b860941a2ef81db,2014,1,28,11292,1129202,30325,"Denver, CO",CO,14252,1425202,34252,"Pasco/Kennewick/Richland, WA",WA,172.0,9.0
53750f8c2b860941a2ef9087,2014,1,27,11292,1129202,30325,"Denver, CO",CO,11884,1188402,31884,"Spokane, WA",WA,173.0,18.0
53750f8c2b860941a2ef9222,2014,1,27,11292,1129202,30325,"Denver, CO",CO,14747,1474703,30559,"Seattle, WA",WA,209.0,33.0
53750f8c2b860941a2efa208,2014,1,5,11292,1129202,30325,"Denver, CO",CO,14747,1474703,30559,"Seattle, WA",WA,185.0,34.0
53750f8d2b860941a2efb1dc,2014,1,4,11292,1129202,30325,"Denver, CO",CO,14747,1474703,30559,"Seattle, WA",WA,197.0,15.0
53750f8d2b860941a2efb3ea,2014,1,6,11292,1129202,30325,"Denver, CO",CO,14747,1474703,30559,"Seattle, WA",WA,167.0,33.0'

exampleData[paste(exampleURL, 'output=csv&type=find&name=flights_to&col=TXFlights&bind.state="CO"', sep='&')] <-
'_id,YEAR,MONTH,DAY_OF_MONTH,ORIGIN_AIRPORT_ID,ORIGIN_AIRPORT_SEQ_ID,ORIGIN_CITY_MARKET_ID,ORIGIN_CITY_NAME,ORIGIN_STATE_ABR,DEST_AIRPORT_ID,DEST_AIRPORT_SEQ_ID,DEST_CITY_MARKET_ID,DEST_CITY_NAME,DEST_STATE_ABR,ACTUAL_ELAPSED_TIME,WEATHER_DELAY
53750f8a2b860941a2ef368d,2014,1,1,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,127.0,0.0
53750f8a2b860941a2ef368e,2014,1,2,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,127.0,0.0
53750f8a2b860941a2ef368f,2014,1,3,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,126.0,0.0
53750f8a2b860941a2ef3690,2014,1,4,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,141.0,0.0
53750f8a2b860941a2ef3691,2014,1,5,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,125.0,0.0
53750f8a2b860941a2ef3692,2014,1,6,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,127.0,0.0
53750f8a2b860941a2ef3693,2014,1,7,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,123.0,null
53750f8a2b860941a2ef3694,2014,1,8,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,120.0,null
53750f8a2b860941a2ef3695,2014,1,9,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,136.0,null
53750f8a2b860941a2ef3696,2014,1,10,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,115.0,0.0
53750f8a2b860941a2ef3697,2014,1,11,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,124.0,null
53750f8a2b860941a2ef3698,2014,1,12,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,122.0,null
53750f8a2b860941a2ef3699,2014,1,13,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,130.0,null
53750f8a2b860941a2ef369a,2014,1,14,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,130.0,null
53750f8a2b860941a2ef369b,2014,1,15,11298,1129803,30194,"Dallas/Fort Worth, TX",TX,11292,1129202,30325,"Denver, CO",CO,124.0,null'

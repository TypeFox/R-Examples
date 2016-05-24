require("chron")
DAILY_DIR <- "DailySourceData"
HOURLY_DIR <- "HourlySourceData"

colnamesHourly <- c( "WBANNO","UTC_DATE","UTC_TIME","LST_DATE","LST_TIME",
                   "CRX_VN"," LONGITUDE","LATITUDE","T_CALC","T_HR_AVG",
                   "T_MAX","T_MIN","P_CALC","SOLARAD","SOLARAD_FLAG",
                   "SOLARAD_MAX","SOLARAD_MAX_FLAG","SOLARAD_MIN","SOLARAD_MIN_FLAG", 
                   "SUR_TEMP","SUR_TEMP_FLAG","SUR_TEMP_MAX","SUR_TEMP_MAX_FLAG", 
                   "SUR_TEMP_MIN","SUR_TEMP_MIN_FLAG","RH_HR_AVG","RH_HR_AVG_FLAG", 
                   "SOIL_MOISTURE_5","SOIL_MOISTURE_10","SOIL_MOISTURE_20", 
                   "SOIL_MOISTURE_50","SOIL_MOISTURE_100", 
                    "SOIL_TEMP_5","SOIL_TEMP_10","SOIL_TEMP_20", 
                    "SOIL_TEMP_50","SOIL_TEMP_100")

colnamesDaily <- c( "WBANNO","LST_DATE",
                   "CRX_VN"," LONGITUDE","LATITUDE","T_DAILY_MAX","T_DAILY_MIN",
                   "T_DAILY_MEAN","T_DAILY_AVE","P_DAILY_CALC","SOLARAD_DAILY", 
                   "SUR_TEMP_DAILY_MAX","SUR_TEMP_DAILY_MIN","SUR_TEMP_DAILY_AVG",  
                     "RH_DAILY_MAX","RH_DAILY_MIN", "RH_DAILY_AVE", "SOIL_MOISTURE_5_DAILY", 
                    "SOIL_MOISTURE_10_DAILY","SOIL_MOISTURE_20_DAILY", 
                    "SOIL_MOISTURE_50_DAILY","SOIL_MOISTURE_100_DAILY", 
                    "SOIL_TEMP_5_DAILY","SOIL_TEMP_10_DAILY","SOIL_TEMP_20_DAILY", 
                    "SOIL_TEMP_50_DAILY","SOIL_TEMP_100_DAILY")

HOURS <- chron(times. = c("00:00:00","01:00:00","02:00:00","03:00:00","04:00:00","05:00:00",
           "06:00:00","07:00:00","08:00:00","09:00:00","10:00:00","11:00:00",
           "12:00:00","13:00:00","14:00:00","15:00:00","16:00:00","17:00:00",
           "18:00:00","19:00:00","20:00:00","21:00:00","22:00:00","23:00:00"), format = "h:m:s")

 

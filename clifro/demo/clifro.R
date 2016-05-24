# Create a public user ----------------------------------------------------

public.user = cf_user()
public.user

# Select datatypes --------------------------------------------------------

# 9am Surface wind (m/s)
wind.dt = cf_datatype(2, 1, 4, 1)

# Daily Rain
rain.dt = cf_datatype(3, 1, 1)

# Daily temperature extremes
temp.dt = cf_datatype(4, 2, 2)

# Combine them together
all.dts = wind.dt + rain.dt + temp.dt
all.dts 

# Select the Reefton Ews station ------------------------------------------

reefton.st = cf_station()
reefton.st

# Submit the query --------------------------------------------------------

# Retrieve all data from ~ six months ago at 9am
reefton.data = cf_query(public.user, all.dts, reefton.st, paste(as.Date(Sys.time()) - 182, "9"))
reefton.data


# Plot the data -----------------------------------------------------------

# Plot the 9am surface wind data (first dataframe in the list) ---
reefton.data[1]

# all identical - although passed to different methods
plot(reefton.data)    #plot,cfDataList,missing-method
plot(reefton.data, 1) #plot,cfDataList,numeric-method
plot(reefton.data[1]) #plot,cfData,missing-method --> plot,cfWind,missing-method

speed_plot(reefton.data)
direction_plot(reefton.data)

# Plot the daily rain data (second dataframe in the list) ---
reefton.data[2]

# With runoff and soil deficit
plot(reefton.data, 2)

# Just plot amount of rain (mm)
plot(reefton.data, 2, include_runoff = FALSE)

# Plot the hourly temperature data (third dataframe in the list) ---
plot(reefton.data, 3)

# Pass an argument to ggplot2::theme
library(ggplot2) # for element_text()
plot(reefton.data, 3, text = element_text(size = 18))

## ---- echo=FALSE---------------------------------------------------------
library(clifro)

## ---- eval = FALSE-------------------------------------------------------
#  me = cf_user("username", "password")

## ---- eval = FALSE-------------------------------------------------------
#  my.dts = cf_datatype(select_1 =     c(7,  4,  3,  2),
#                       select_2 =     c(1,  2,  1,  1),
#                       check_box = list(3,  1,  1,  4),
#                       combo_box =    c(NA, NA, NA, 1))
#  my.dts

## ---- eval = FALSE-------------------------------------------------------
#  my.stations = cf_station(5814, 4241, 2112, 1962)
#  my.stations[, 1:5]

## ---- eval = FALSE-------------------------------------------------------
#  cf.datalist = cf_query(user = me,
#                         datatype = my.dts,
#                         station = my.stations,
#                         start_date = "2012-01-01 00",
#                         end_date = "2014-01-01 00")
#  cf.datalist

## ---- eval = FALSE-------------------------------------------------------
#  # Load the ggplot2 library for element_text() and geom_smooth() functions
#  library(ggplot2)
#  
#  # Increase the text size to 16pt and add a loess smoother with a span equal to a
#  # quarter of the window
#  plot(cf.datalist, ggtheme = "bw", text = element_text(size = 16)) +
#    geom_smooth(method = "loess", span = 1/4)

## ---- eval = FALSE-------------------------------------------------------
#  # Try a different ggtheme
#  plot(cf.datalist, 2, ggtheme = "linedraw")

## ---- eval = FALSE-------------------------------------------------------
#  # Try yet another ggtheme
#  plot(cf.datalist, 3, ggtheme = "light")
#  
#  # Or only plot the rainfall data
#  # plot(cf.datalist, 3, ggtheme = "light", include_runoff = FALSE)

## ---- eval = FALSE-------------------------------------------------------
#  # Defaults to windrose
#  plot(cf.datalist, 4, n_col = 2)

## ---- eval = FALSE-------------------------------------------------------
#  # Plot the wind speeds through time, choose the 'classic' ggtheme and
#  # allow the y-axis scales to differ for each station
#  speed_plot(cf.datalist, 4, ggtheme = "classic", scales = "free_y")
#  
#  # Plot wind direction contours through time
#  direction_plot(cf.datalist, 4, n_col = 2)

## ---- eval = FALSE-------------------------------------------------------
#  # Export the data as separate CSV files to the current working directory
#  for (i in seq_along(cf.datalist))
#    write.csv(cf.datalist[i],
#              file = tempfile(paste0(cf.datalist[i]@dt_name, "_"),
#                              tmpdir = normalizePath("."),
#                              fileext = ".csv"),
#              na = "", row.names = FALSE)
#  
#  # Each dataset is saved separately here:
#  getwd()


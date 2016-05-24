Aide4 <- function () {

aide4 <- tktoplevel()
tktitle(aide4) <- "Help3 - TimeSeries_building"
AIDE4 <- tktext(aide4, height=50, width=120)
tkpack(AIDE4)
tkinsert(AIDE4, "end", paste("This is the forth panel of the interface" ,"\n\n",

"Here, you can observe some of the data characteristics or proceed to temporal","\n",
"trend analysis.","\n\n",

"      _The <Diagnostics> box, propose different options that give some extra informations about your time series.","\n\n",

"       <Spectrum analysis> estimates the spectral density of the time series. Select the pick value in the figure", "\n",
"       and the interface will display the corresponding cycle. If your time step","\n",
"       is monthly and the major frequency is 1 then the time series have a cycle of 12 months." ,"\n",
"       A time series could have more than one identified frequency.","\n\n",

"       <Autocorrelation> compute an autocorrelation function on the time series and display","\n",
"       a figure of the result.","\n\n",

"       <Shapiro normality test>, test the normality of the time series distribution and display a Q-Q Plot.","\n\n",

"       <Anomaly (color.plot)> draw a color plot by year each time.step (months or weeks)","\n",
"       minus the mean of the time.step of all years. Red colors show positive anomalies and","\n",
"       blue colors negative anomalies. ","\n\n",

"       <Anomaly (bar.plot)> display an anomaly barplot (years to weeks time step).","\n",
"       Red colors show positive anomalies and blue colors negative anomalies.","\n\n",

"       <Seasonal decomposition> decompose a time series into seasonal, trend and irregular components using","\n",
"       loess and display the results in a figure (with remainder=residuals).","\n\n\n",

"      _In the <Trend analysis> box you can perform a temporal trend analysis on the time series","\n",
"       you compute through the interface.","\n\n",

"       <Cumulative sum> let you show the cusum curve of your time series ('pastecs' package),","\n",
"       select interesting periods and perform Kendall family test on these sub-series.","\n\n",

"       <Seasonal Trend> performs a Seasonal Mann-Kendall test ('wq' package) on your time series.","\n",
"       Mann-Kendall is therefore compute on each season or time step you choose (i.e. weeks,","\n",
"       months...)","\n\n",

"       <Global Trend> perform a Seasonal Mann Kendall on the entire time series, taking into account","\n",
"       the seasonality of your data.","\n\n",

"       <Trend based on LOESS> fit a polynomial surface determined by one or more numerical predictors,","\n",
"       using local fitting. A Global Trend with Sen's slope is perform on this fitting.","\n\n",

"       <Trend based on mixing diagram> used normalized value of nutrient concentration at salinity npsu for each year.","\n",
"       This is obtain by doing a linear regression between salinity and nutrient concentration of this year.","\n",
"       The nutrient concentration predicted by the regression at salinity npsu is kept for the year.","\n",
"       The final curve is obtain by plotting the predicted nutrient concentration for each year.","\n",
"       This is only relevant for nutrients.", sep="" ))
       
tkpack(tk2button(aide4, text="Ok", command=function(){ tkdestroy(aide4) }, width=20 ), side="bottom")
tkconfigure(AIDE4, state="disabled", background="white")

}
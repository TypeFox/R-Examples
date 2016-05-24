Aide3 <- function () {

aide3 <- tktoplevel()
tktitle(aide3) <- "Help3 - TimeSeries_building"
AIDE3 <- tktext(aide3, height=40, width=120)
tkpack(AIDE3)
tkinsert(AIDE3, "end", paste("This is the third panel of the interface." ,"\n\n",

"Here, you can interact with your data and choose the methods of aggregation to build","\n",
"a regulated time series, essential to perform temporal trend analysis." ,"\n\n",

"    _In the <Data interaction> box you can choose to replace the missing values present" ,"\n",
"     in your data. This can offer you the possibility to perform more analysis or diagnostics","\n",
"     (the ones labelled with an *). It works in two times, first when more than 3 values are"   ,"\n",
"     successively missed, they are replace by the median of equivalent cycle (i.e. months or weeks)."    ,"\n",
"     For small interval missed (less than 3 values), data are replace by a linear regression of"    ,"\n",
"     the data before and after (data are generally link in time series due to autocorrelation)."  ,"\n",
"     However when missing values are numerous (here higher than 1/20 of the data quantity)"    ,"\n",
"     this is not a good idea to replace them and a warning message is display. Replacing to much"  ,"\n",
"     missing values can create aberrant values and modify to much the original data series." ,"\n\n",

"     You can also remove the outliers present in you data, the 'Show boxplot' button display"   ,"\n",
"     box and whiskers by years with outliers as open circle. Outliers are calculated as"         ,"\n",
"     [quantile 0.75 + 1.5(quantile 0.75 - quantile 0.25)] and" ,"\n",
"     [quantile 0.25 - 1.5(quantile 0.75 - quantile 0.25)]"    ,"\n",
"     Be careful, if years are missed in your database, they will not appear in the boxplot." ,"\n\n",

"    _The second box let you choose the frequency at which you want to aggregate your data."    ,"\n",
"     The shorter the time step is, the most information you will probably missed but the less"   ,"\n",
"     missing values you will obtain."  ,"\n",
"     The option 'Guidance to choose the time step' compute the mean time (in days) that" ,"\n",
"     separate two measurement in your database. As a function of the result it suggest a frequency" ,"\n",
"     to use (daily if the result <=5, semi-fortnightly <=10 days; fortnightly <=24 days;"  ,"\n",
"     monthly <=60 days and yearly >60 days)"  ,"\n",
"     Monthly - Climato aggregate the respective months of all years. "  ,"\n",
"     The auto method use the result of the guidance to automatically perform the aggregation,"  ,"\n",
"     this is the default method." ,"\n\n",

"    _The third box do the same things as the second one but with the mathematical method of" ,"\n",
"     aggregation instead of frequency."  ,"\n",
"     <Guidance to choose the method> perform an Anova with a Dunnett's post-hoc test between the"  ,"\n",
"     rawdata and the dataset obtain by each method and select the method with the highest p.value" ,"\n",
"     (less significant differences)." ,"\n\n",

"    _The <Show regularised time series> buttons let you show the table, the figure or the descriptive statistics"  ,"\n",
"     of your data passed through the regulation methods you choose. Figure and tables are automatically saved.", sep=""))
     
tkpack(tk2button(aide3, text="Ok", command=function(){ tkdestroy(aide3) } , width=20 ), side="bottom")
tkconfigure(AIDE3, state="disabled", background="white")

}
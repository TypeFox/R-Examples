Aide2 <- function () {

aide2 <- tktoplevel()
tktitle(aide2) <- "Help2 - Parameters_selection"
AIDE2 <- tktext(aide2, height=37, width=120)
tkpack(AIDE2)
tkinsert(AIDE2, "end", paste("This is the second panel of the interface", "\n\n",

"In this panel you will be able to select the parameter, the station(s) and the dates","\n",
"you want to test.","\n\n",

"     _In the first selection panel, the first list show you the categories (Stations, Taxa,...) present in your file.","\n",
"      You can select them one by one or by group in the first list by using the top blue arrow button.","\n",
"      Select < -All- > to take all categories into account.","\n",
"      The selected categories appear in the second list.",  "\n",
"      To remove category from your selection, select them in the second list and use the bottom blue arrow button.","\n\n",

"     _In the second selection panel, the first list show you all parameters present in your database.", "\n",
"      You can select only one parameter; select it and clicking on the arrow button.",   "\n",
"      Selecting a second one and clicking on the arrow button will replace the previously selected", "\n",
"      parameter.","\n\n",

"     _You can select also the salinity and depth intervals at which you want analyse your parameter.", "\n",
"      To proceed, just slide the bars, if you let the bars at minimal and maximal values,", "\n",
"      parameter values at all salinities or depths will be take into account (including at missing values", "\n",
"      of salinity or depth).",  "\n",
"      If you select different values, values of the parameter at missing salinity or depth (NA) will be removed","\n",
"      from the treatment.","\n\n",

"     _On the right side you can select years and months of treatment. By default, minimal and maximal years","\n",
"      of the database are display, you can change them by clicking on the arrow on the right side of the box.", "\n",
"      Only the months present in your database are display by default, you can delete or add months but each", "\n",
"      month value must be separated by a unique space.","\n\n",

"     _The <Summary> button display descriptive statistics of the selected parameters in the results window", "\n",
"      (in right side of the interface) and in the front of the screen (in a data frame).", sep=""))
      
tkpack(tk2button(aide2, text="Ok", command=function(){ tkdestroy(aide2) } , width=20), side="bottom")
tkconfigure(AIDE2, state="disabled", background="white")

}
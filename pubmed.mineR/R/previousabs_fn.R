previousabs_fn = function(yrs_to_exclude, theme,parentabs){
check = parentabs@Journal
indices=NULL;for (i in 1:length(yrs_to_exclude)){formatyr = paste(". ",yrs_to_exclude[i], sep ="")
check2 = regexpr(formatyr, check, fixed=T)
check3 = which(check2 != -1)
indices = c(indices, check3)}
new_abs = subsetabs(parentabs,-indices)
new_abs_theme = searchabsL(new_abs, include=theme)
return(new_abs_theme)}

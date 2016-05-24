currentabs_fn = function(yr_to_include, theme,parentabs){
check = parentabs@Journal
formatyr = paste(". ",yr_to_include[1], sep ="")
check2 = regexpr(formatyr, check, fixed=T)
check3 = which(check2 != -1)
new_abs = subsetabs(parentabs,check3)
new_abs_theme = searchabsL(new_abs, include=theme)
return(new_abs_theme)}

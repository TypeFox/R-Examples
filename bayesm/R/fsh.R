fsh=function() 
{
# 
# P. Rossi
# revision history: 3/27/05
#
# Purpose:
#  function to flush console (needed only under windows)
#
if (Sys.info()[1] == "Windows") flush.console()
return()
}


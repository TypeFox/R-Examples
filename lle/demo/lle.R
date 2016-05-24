#small demo to show some lle examples 

Pause... <- function(x) invisible(readline())
#example: scurve
#embedd a three-dimensional s-curve area
#correctly into two-dimensional space
lle_scurve()
Pause...("Press enter to continue")
dev.off()

#example: spiral
#vary the number of twists of a 
#three-dimensional spiral and monitor
#the change of intrinsic dimension
lle_spiral()
Pause...("Press enter to continue")
dev.off()

#example: soundfile
#read a soundfile window-wise and
#embedd the windows into lower-
#dimensional space to determine the
#number of modes
lle_sound()
Pause...("Press enter to continue")
dev.off()
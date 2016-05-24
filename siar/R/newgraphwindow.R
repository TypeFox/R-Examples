newgraphwindow <-
function() {
# This function aims to create graphics windows on different OS's where possible
graphoutcome <- try(windows(),silent=TRUE)
if(length(graphoutcome) > 0) graphoutcome <- try(quartz(),silent=TRUE)
if(length(graphoutcome) > 0) graphoutcome <- try(x11(),silent=TRUE)
if(length(graphoutcome) > 0) print("Graphics unavailable on this workstation")
}

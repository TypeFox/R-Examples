summary.BchronCalibratedDates <-
function(object, prob=c(50,95,99), ..., digits = max(3, getOption("digits")-3)) {
  x = object
  for(i in 1:length(x)) {
    currdate = list(x=x[[i]]$ageGrid,y=x[[i]]$densities)
    cat('Highest density regions for',names(x)[i],'\n')
    
    # hdr sometimes fails to this while loop removes that fail step
    bad = TRUE
    while(bad) {
      my_hdr = try(hdrcde::hdr(den=currdate,prob=prob)$hdr,silent=TRUE)
      if(class(my_hdr)!='try-error') bad = FALSE
    }
    print(round(my_hdr,1))
    cat('\n')
  }
}

# rp.checkbox demo code

panel <- rp.control(checkboxvar=list(TRUE,FALSE,TRUE))
rp.checkbox(panel, 
  checkboxvar,
  function(panel)
  { 
    rp.messagebox(panel$checkboxvar); print(panel$checkboxvar); print(panel$checkboxvar[1]); panel
  },
  title="Choose", 
  labels=list("File","Option","Exit")#, 
#  initval=list(FALSE,TRUE,FALSE) # this can be done either here or in the initial creation of the control
)


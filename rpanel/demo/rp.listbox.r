panel <- rp.control(option="Option 1")
rp.listbox(panel, option, title="Choose", labels=list("File","Option 1","Option 2","Option 3","Exit"), rows=3, action=function(panel){ rp.messagebox(panel$option); print(panel$option); panel }) # can also use initval="Option 1", 

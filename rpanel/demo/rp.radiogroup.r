panel <- rp.control(option1="Option", option2="Exit")
rp.radiogroup(panel, option1, pos=list(column=0, row=0), title="Choose", labels=c("File","Option","Exit"), action=function(panel){ print(panel$option1); print(panel$option2); panel })
rp.radiogroup(panel, option2, pos=list(column=0, row=1), background="Green", title="Choose", labels=c("File","Option","Exit"), action=function(panel){ print(panel$option1); print(panel$option2); panel })

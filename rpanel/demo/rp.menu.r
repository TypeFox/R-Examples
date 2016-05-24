panel <- rp.control(title="Menu example", size=c(200,50))
rp.menu(panel, option, labels=list(list("File","Quit"), list("Edit","Copy","Cut","Paste")), "None", function(panel) { rp.messagebox(panel$option); panel }, font="Arial", foreground="Navy", background="Green")


data(USArrests)

# demo(dg.prcomp)

# demo(defaultObjects)

# demo(usermenus)

source(paste(system.file(package = "dynamicGraph"), 
             "demo.source/dg.prcomp.R", sep = "/"))

source(paste(system.file(package = "dynamicGraph"), 
             "demo.source/defaultObjects.R", sep = "/"))

source(paste(system.file(package = "dynamicGraph"), 
             "demo.source/usermenus.R", sep = "/"))

dg.prcomp(data = USArrests, check.plot = TRUE, UserMenus = Menus)

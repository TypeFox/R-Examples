
data(USArrests)

# demo(dg.prcomp)

# demo(defaultObjects)

# demo(usermenus)

source(paste(system.file(package = "dynamicGraph"), 
             "Demo/dg.prcomp.R", sep = "/"))

source(paste(system.file(package = "dynamicGraph"), 
             "Demo/defaultObjects.R", sep = "/"))

source(paste(system.file(package = "dynamicGraph"), 
             "Demo/usermenus.R", sep = "/"))

dg.prcomp(data = USArrests, check.plot = TRUE, UserMenus = Menus)

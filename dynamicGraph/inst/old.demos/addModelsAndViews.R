
# Test of "link" and "add":

# demo("startup", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "Demo/startup.R", sep = "/"))

From <- c(1, 2, 3)
To   <- c(2, 3, 1)

Z <- DynamicGraph(V.Names[1:3], V.Types[1:3], 
                  from = From, to = To, object = Object, UserMenus = Menus, 
                  returnLink = TRUE, width = 150, height = 200, margin = 200,
		  title = "Z")

From <- c(2, 3)
To   <- c(3, 1)

W <- DynamicGraph(from = From, to = To, object = Object, UserMenus = Menus, 
                  returnLink = TRUE, width = 150, height = 200, margin = 200, 
                  title = "W", frameModels = Z, addModel = TRUE)

V <- DynamicGraph(from = From, to = To, object = Object, UserMenus = Menus, 
                  returnLink = TRUE, width = 150, height = 200, margin = 200, 
                  frameModels = W, frameViews = W@models[[2]], 
                  title = "V", addView = TRUE, viewType = "Factor")

From <- 1
To   <- 2

U <- DynamicGraph(from = From, to = To, object = Object, UserMenus = Menus, 
                  returnLink = TRUE, width = 150, height = 200, margin = 200, 
                  title = "U", frameModels = Z, frameViews = Z@models[[1]], 
                  graphWindow =  Z@models[[1]]@graphs[[1]], addModel = TRUE, 
                  overwrite = TRUE)

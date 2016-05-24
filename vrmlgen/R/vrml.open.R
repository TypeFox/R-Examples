
# create a new VRML environment in which multiple
# plotting functions can be combined
`vrml.open` <-
function (filename = "out.wrl", col = "white", navigation = NULL, 
    scale = 1, fov = 0.785, pos = rep(scale + 8, 3), dir = c(0.19, 
        0.45, 0.87, 2.45), html.embed = "out.html", hwidth = 1200, 
    hheight = 800) 
{
    
    # set the global parameters of the 3D scene
    .vrmlgenEnv <<- new.env()
    .vrmlgenEnv$type <- "vrml"
    .vrmlgenEnv$VRMLDir <- tempdir()
    .vrmlgenEnv$filename <- filename
    .vrmlgenEnv$col <- col
    .vrmlgenEnv$html <- html.embed
    .vrmlgenEnv$hheight <- hheight
    .vrmlgenEnv$hwidth <- hwidth
    .vrmlgenEnv$navigation <- navigation
    .vrmlgenEnv$scale <- 1 * scale
    curdir <- getwd()
    setwd(.vrmlgenEnv$VRMLDir)
    
    # write the VRML header and the main 3D scene settings
    write("#VRML V2.0 utf8\n", file = filename, append = FALSE)
    
    # set the viewpoint
    write(paste("\nViewpoint {\n\tfieldOfView", fov, "\n\tposition", 
        pos[1], pos[2], pos[3], "\n\torientation", dir[1], dir[2], 
        dir[3], dir[4], "\n\tjump TRUE\n\tdescription \"viewpoint1\"\n}\n", 
        sep = " "), file = filename, append = TRUE)
    
    # set the navigation type
    if (!is.null(navigation)) 
        write(paste("\nNavigationInfo { type \"", navigation, 
            "\" }\n", sep = ""), file = filename, append = TRUE)
            
    # set the background color
    bg_rcol <- (col2rgb(col)/255)[1]
    bg_gcol <- (col2rgb(col)/255)[2]
    bg_bcol <- (col2rgb(col)/255)[3]
    
    write(paste("Background {\n\t skyColor [\n\t\t ", bg_rcol, 
        bg_gcol, bg_bcol, " \n\t]\n}", sep = " "), file = filename, 
        append = TRUE)
    
    # return to the current directory    
    setwd(curdir)
}


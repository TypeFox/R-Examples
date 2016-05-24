
# create a new Livegraphics3D environment in which multiple
# plotting functions can be combined
`lg3d.open` <-
function (filename = "out.m", col = "white", scale = 1, html.embed = "out.html", 
    hwidth = 1200, hheight = 800, ambientlight = 0.5) 
{
    
    # set the global parameters of the 3D scene
    .vrmlgenEnv <<- new.env()
    .vrmlgenEnv$type <- "lg3d"
    .vrmlgenEnv$VRMLDir <- tempdir()
    .vrmlgenEnv$filename <- filename
    .vrmlgenEnv$col <- col
    .vrmlgenEnv$html <- html.embed
    .vrmlgenEnv$hheight <- hheight
    .vrmlgenEnv$hwidth <- hwidth
    .vrmlgenEnv$scale <- scale
    .vrmlgenEnv$ambientlight <- ambientlight

		# write the Livegraphics 3D header        
    curdir <- getwd()
    setwd(.vrmlgenEnv$VRMLDir)
    write("Graphics3D[\n{\n", file = filename, append = FALSE)
    setwd(curdir)
}


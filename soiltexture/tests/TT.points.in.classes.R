require( "soiltexture" ) 

# Create a dummy data frame of soil textures:
my.text <- data.frame( 
    "CLAY"  = c(05,60,15,05,25,05,25,45,65,75,13,47), 
    "SILT"  = c(05,08,15,25,55,85,65,45,15,15,17,43), 
    "SAND"  = c(90,32,70,70,20,10,10,10,20,10,70,10), 
    "OC"    = c(20,14,15,05,12,15,07,21,25,30,05,28)  
)   #

# Display the table:
my.text

# Classify according to the HYPRES / European Soil Map classification
TT.points.in.classes( 
    tri.data    = my.text[1:5,], 
    class.sys   = "HYPRES.TT"  
)   #

# Classify according to the USDA classification
TT.points.in.classes( 
    tri.data    = my.text[1:5,], 
    class.sys   = "USDA.TT"  
)   #

# Classify according to the HYPRES / European Soil Map classification, 
#   returns logical values
TT.points.in.classes( 
    tri.data    = my.text[1:5,], 
    class.sys   = "HYPRES.TT", 
    PiC.type    = "l" 
)   #

# Classify according to the HYPRES / European Soil Map classification, 
#   returns text
TT.points.in.classes( 
    tri.data    = my.text[1:5,], 
    class.sys   = "HYPRES.TT", 
    PiC.type    = "t" 
)   #

# Classify according to the HYPRES / European Soil Map classification, 
#   returns text, 
#   custom class separator in case of points belonging to 
#   several classes.
TT.points.in.classes( 
    tri.data    = my.text[1:5,], 
    class.sys   = "HYPRES.TT", 
    PiC.type    = "t", 
    collapse    = ";"
)   #

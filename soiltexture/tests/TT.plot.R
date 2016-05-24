require( soiltexture ) 

# ::: Texture triangles without data

# :: Base plot (HYPRES / European Soil Map triangle) 
TT.plot() 

# same as
TT.plot( class.sys = "HYPRES.TT" ) 

# :: Same plot, but with USDA texture triangle 
TT.plot( class.sys = "USDA.TT" ) 

# :: Same plot, but with a color gradient 
TT.plot( 
    class.sys       = "USDA.TT", 
    class.p.bg.col  = TRUE
)   #

# :: No texture classification system
TT.plot( class.sys = "none" ) 

# ::: Texture triangles with texture data 

# :: 1st create a dummy texture dataset 
my.text <- data.frame( 
    "CLAY"  = c(05,60,15,05,25,05,25,45,65,75,13,47), 
    "SILT"  = c(05,08,15,25,55,85,65,45,15,15,17,43), 
    "SAND"  = c(90,32,70,70,20,10,10,10,20,10,70,10), 
    "OC"    = c(20,14,15,05,12,15,07,21,25,30,05,28)  
)   #

# :: And plot it on a French Aisne texture triangle
#    with a title
TT.plot( 
    class.sys   = "FR.AISNE.TT", 
    tri.data    = my.text, 
    main        = "Soil texture data" 
)   #

# ::: Bubble plots (4th variable) 

# :: 1st generate a dummy texture dataset with a 4th variable 
#    with TT.dataset() 
rand.text   <- TT.dataset( n = 100, seed.val = 1980042401 ) 

# :: Plot the dummy dataset as a bubble plot
TT.plot( 
    class.sys   = "none", 
    tri.data    = rand.text, 
    z.name      = "Z", 
    main        = "Soil texture triangle and Z bubble plot" 
)   #

# ::: Test all the texture triangles
TT.plot( class.sys = "none" )           # no classification 
TT.plot( class.sys = "HYPRES.TT" )      # HYPRES / European Soil Map 
TT.plot( class.sys = "USDA.TT" )        # USDA 
TT.plot( class.sys = "FR.AISNE.TT" )    # French Aisne 
TT.plot( class.sys = "FR.GEPPA.TT" )    # French GEPPA 
TT.plot( class.sys = "DE.BK94.TT" )     # Germany 
TT.plot( class.sys = "DE.SEA74.TT" )    # German SEA 1974 
TT.plot( class.sys = "DE.TGL85.TT" )    # German TGL 1985 
TT.plot( class.sys = "UK.SSEW.TT" )     # UK 
TT.plot( class.sys = "BE.TT" )          # Belgium 
TT.plot( class.sys = "CA.FR.TT" )       # Canada (fr) 
TT.plot( class.sys = "CA.EN.TT" )       # Canada (en) 
TT.plot( class.sys = "AU2.TT" )         # Australian 
TT.plot( class.sys = "ISSS.TT" )        # ISSS 
TT.plot( class.sys = "ROM.TT" )         # Romanian 
TT.plot( class.sys = "USDA1911" )       # USDA 1911 (M. Whitney, 1911)
TT.plot( class.sys = "BRASIL.TT" )      # Brasil (Lemos & Santos 1996)
TT.plot( class.sys = "SiBCS13.TT" )     # Brasil (Lemos & Santos 1996)

#   Triangles with special characters
#   (may not work on all platforms + some accents can be missing)
try( TT.plot( class.sys = "PL.TT" ) )   # Polish 


# ::: Test all the languages:
TT.plot( class.sys = "USDA.TT", lang = "en" )  # English, default 
TT.plot( class.sys = "USDA.TT", lang = "fr" )  # French 
TT.plot( class.sys = "USDA.TT", lang = "de" )  # German 
TT.plot( class.sys = "USDA.TT", lang = "es" )  # Spanish 
TT.plot( class.sys = "USDA.TT", lang = "it" )  # Italian 
TT.plot( class.sys = "USDA.TT", lang = "nl" )  # Dutch 
TT.plot( class.sys = "USDA.TT", lang = "fl" )  # Dutch (Belgium) / Flemish 
TT.plot( class.sys = "USDA.TT", lang = "se" )  # Swedish 
TT.plot( class.sys = "USDA.TT", lang = "ro" )  # Romanian 

#   Languages with special characters
#   (may not work on all platforms + some accents can be missing)
try( TT.plot( class.sys = "USDA.TT", lang = "pl"  ) ) # Polish 
try( TT.plot( class.sys = "USDA.TT", lang = "pt"  ) ) # Portuguese 
try( TT.plot( class.sys = "USDA.TT", lang = "es2" ) ) # Spanish
try( TT.plot( class.sys = "USDA.TT", lang = "ro2" ) ) # Romanian


# Environment for storing, hiding and protecting internal variables and functions.
TT.env <- new.env() 

# assign( 
#     x       = "TT.env", 
#     value   = new.env() 
# )   #






# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# | LIST: TT.par                        |
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# [ TT.par :: A list of default parameters, contained in a list, stored in the TT.env environment (~ invisible)
assign( 
    envir   = TT.env,
    x       = "TT.par", 
    value   = list( 
        # +-------------------------------------------------------------------------+
        # | SCRIPT PARAMETERS SPECIFICATION:                                        |
        # +-------------------------------------------------------------------------+
        #
        # +---------------------------------+
        # | TRIANGLE GEOMETRY               |
        # +---------------------------------+
        #
        # Geometry of the texture triangle, in cases where no specific texture triangle 
        #   is used. Be aware that some of these default parameters will be overwritten 
        #   by triangle specific parameters.
        blr.clock       = rep(T,3), # clockwise axes (T) ? or counterclock ? (F)
                                    #   clock wise mean bottom axis is oriented from right to left, 
                                    #   left axis is oriented bottom to top and right axis is 
                                    #   oriented from top to bottom.
                                    #   Can also be c(F,T,NA) or c(T,NA,F) or c(NA,F,T)
        # Top, Left and Right ANGLES of the triangle, in DEGREES
        tlr.an          = c(60,60,60), 
                                    #   For nicer triangle, higher angle must be on the LEFT in the case
                                    #   of a clock orientation, and higher angle must be on the RIGHT 
                                    #   in the case of a counter-clock orientation
                                    # Bottom, Left and Right TEXTURES of the triangles
        # 
        # NB: in the default options, coordinates are always expressed as a "fraction" (0 to 1 values)
        # 
        # Limits of the triangle to be drawn (get rid of this???)
        base.frame      = data.frame( 
            # Rows (here vertical)      = triangle submits
            # Coulmns (here horizontal) = coordinates of the submits for bootom, left and right variables 
            "b" = c( 1, 0, 0 ), 
            "l" = c( 0, 1, 0 ), 
            "r" = c( 0, 0, 1 ) 
        ),  #
        # 
        # Limits (min,max values) of the different variables
        b.lim           = c(0,1),   # Bootom variable
        l.lim           = c(0,1),   # Left   variable
        r.lim           = c(0,1),   # Right  variable
        # 
        # Limits "tolerance", as a fraction of the triangle maximum range (0 to 1)
        lim.tol         = 0.1, 
        # 
        # +---------------------------------+
        # | TRIANGLE BASE CONTENT           |
        # +---------------------------------+
        # See also pre-formatted triangles with classe below
        #
        # Bottom, Left and Right TeXtures (= which texture 
        #   displaying on B, L, R axis). MUST belong to CLAY SILT 
        #   and SAND (the order is free).
        blr.tx          = c("SAND","CLAY","SILT"), 
        # 
        # Clay, Silt and Sand (columns) NAMES in the input soil 
        #   texture data table (tri.data). 
        #   example: c("ARGILE","LIMON","SABLE") 
        css.names       = c("CLAY","SILT","SAND"), 
        #
        # A vector of 3 character strings, or expressions, for the 
        #   LABELS of Clay, Silt and Sand, respectively. If 
        #   non-null, will overwrite any default label (in any lang).
        css.lab         = NULL, 
        #
        #blr.lab         = NULL,     # a vector of 3 character strings 
        #                           # or expressions, with the labels 
        #                           # of bottom, left and right axis.
        #
        lang            = "en", 
        #
        # blr.psize.lim = data.frame( 
        #   "SAND"  = c(50,2000), 
        #   "CLAY"  = c(0,2), 
        #   "SILT"  = c(2,50) 
        # ),    #
        # 
        unit.ps         = quote(bold(mu) * bold('m')), 
        unit.tx         = quote(bold('%')), 
        #
        # Input data:
        text.sum        = 100,      # Value of the SUM OF CLAY SILT SAND TEXTURE: either 1 or 100 (or fancy)
        text.tol        = 1/1000,   # Error tolerance on the sum of the 3 particle size classes
        tri.sum.tst     = TRUE,     # Perform a sum test on the tri-variable data (sums == text.sum) ??
        tri.pos.tst     = TRUE,     # Test if all tri-variable data are positive ??
        # 
        # +---------------------------------+
        # | INTERNATIONALISATION            |
        # +---------------------------------+
        lang.par    = data.frame( 
            "lang"  = c(    "en",                           "fr",                       "it", 
                            "es",                           "de",                       "nl",
                            "se",                           "fl",                       "ro" ),   
            #
            "CLAY"  = c(    "\"Clay\"",                     "\"Argile\"",               "\"Argilla\"", 
                            "\"Arcilla\"",                  "\"Ton\"",                  "\"Lutum\"", 
                            "\"Ler\"",                      "\"Klei\"",                 "\"Argila\"" ),   
            #
            "SILT"  = c(    "\"Silt\"",                     "\"Limon\"",                "\"Limo\"", 
                            "\"Limo\"",                     "\"Schluff\"",              "\"Silt\"", 
                            "\"Silt\"",                     "\"Leem\"",                 "\"Praf\"" ),   
            #
            "SAND"  = c(    "\"Sand\"",                     "\"Sable\"",                "\"Sabbia\"", 
                            "\"Arena\"",                    "\"Sand\"",                 "\"Zand\"", 
                            "\"Sand\"",                     "\"Zand\"",                 "\"Nisip\"" ),   
            #
            "TT"    = c(    "\"Texture triangle\"",         "\"Triangle de texture\"",  "\"Triangolo della tessitura\"", 
                            "\"Triangulo de textura\"",     "\"Bodenartendiagramm\"",   "\"Textuurdriehoek\"", 
                            "\"Texturtriangel\"",           "\"Textuurdriehoek\"",      "\"Diagrama triunghiulara a texturii\"" ), 
            
            # NOTE: accents removed!
            
            stringsAsFactors    = FALSE  
        ),  #
        # Acknowledgments: Rosca Bogdan, from the Romanian Academy, Iasi Branch, Geography team, provided the Romanian translation (thanks!).
        # 
        # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
        # | TRIANGLE CUSTOMISATION          |
        # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
        # 
        # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
        # | General graphical parameters
        # 
        # FONTS:
        font            = NULL,     # for plotted points and text
        font.axis       = 2,        # for ticks labels
        font.lab        = 2,        # for axis labels (arrows) and texture classes labels
        font.main       = 2,        # main title 
        # 
        # COLORS:
        bg              = NULL,     # plot and plotted points background
        fg              = NULL,     # for plotted points and text (foreground) 
        col             = NULL,     # for plotted points and text 
        col.axis        = NULL,     # for ticks, triangle "axis/frame" (but not texture classes or grid)
        col.lab         = NULL,     # for axis labels and arrows (but not texture classes labels)
        col.main        = NULL,     # main title 
        # NOTICE: grid lines and "frame" background colors are set with other options (non generic)
        # 
        # CEX:
        cex             = 1.5,      # for plotted points and text 
        cex.axis        = 1.5,      # for ticks labels
        cex.lab         = 1.5,      # for axis labels and texture classes labels
        cex.main        = 1.5,      # main title 
        # 
        # LWD: 
        lwd             = 3,        # for plotted lines (not implemented yet?) 
        #
        #
        #
        # warning. belw, not in par()
        lwd.axis        = 3,        # for ticks, triangle "axis/frame" and grid
        lwd.lab         = 3,        # for axis arrows and texture classes polygons
        # 
        # FAMILY: 
        family.op       = NULL, 
        #
        # graph margins (as in par("mar")):
        new.mar         = NULL, 
        #
        # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
        # | Specific graphical parameters
        #
        # TRIANGLE FRAME parameters:
        frame.bg.col    = NULL,             # triangle frame background color
        #
        # Texture values _at_ which starting grid or thicks points (values between 0 and 1!!!)
        at              = seq(from=.1,to=.9,by=.1), 
        # 
        # CLASSES SYSTEM (polygons / Texture Triangle) used by default:
        class.sys       = "HYPRES.TT", 
        class.lab.col   = NULL,             # Color of classes names (abreviation)
        class.p.bg.col  = FALSE,            # Fill classes polygon with color gradient ???
        class.p.bg.hue  = 0.04,             # Hue (unique) of the classes polygon color gradient
        class.line.col  = NULL,             # Classes lines (foreground) color
        #                                   # 0.04 is salmon-pink, 0.21 is olive green, 0.57 is sky blue, 0.72 a nice purple-blue
        class.lty       = NULL, 
        class.lab.show  = "abr",            # Show
        #
        # GRID LINES parameters:
        grid.show       = TRUE, 
        grid.col        = NULL,             # grid line colors (fg) # default in fact gray(.5)
        grid.lty        = 1,                # Grid line type
        #
        # TICKS MARK parameters:
        ticks.shift         = 02.5/100,     # Value of the ticks shift from the border of the triangle: ]0;1[ (but 1 extra big)
        ticks.lab.shift     = 05.0/100,     # Value of the ticks label shift from the border of the triangle: ]0;1[ (but 1 extra big) 
        #
        # ARROWS parameters:
        arrows.show         = TRUE, 
        arrows.lims         = c(0.150,0.450), 
        arrows.head.shift   = 05.0/100,     # Distance between the arrow head and the triangle frame
        arrows.base.shift   = 11.0/100,     # Distance between the arrow base and the triangle frame
        arrows.text.shift   = 11.0/100,     # Distance between the arrow label/text and the triangle frame
        arrows.text.shift2  = 10.0/100,     # Distance between the arrow "turn" and the text (parallel to axis)
        arrows.lty          = 1, 
        #
        # POINTS-DATA parameters:
        points.type         = "p", 
        pch                 = NULL,  # Added 20090617 
        # 
        # TEXT DATA parameters:
        adj             = c(0.5,0.5), 
        pos             = NULL, 
        offset          = 0, 
        #
        # +---------------------------------+
        # | TT.points.in.classes parameters |
        # +---------------------------------+
        PiC.type            = c("n","l","t")[1], 
        #                                   # Type of output for TT.points.in.classes
        #                                   # "n" stands for numeric (0 = out, 1 = in, 3 = border) 
        #                                   # "l" stands for logical (F = out, T = in) 
        #                                   # "t" stands for text (concatenated texture class symbol)
        # 
        # +---------------------------------+
        # | z: the 4th variable             |
        # +---------------------------------+
        # [ z may be plotted as a bubble chart or as an interpolated map
        z.type              = c("bubble","map")[1], 
        z.col.hue           = 0.21, 
        z.cex.range         = c(1.00,4.00), 
        # z.cex.min         = 0.50, 
        z.pch               = c(16,1),  # Respectively for color fill and border
        #
        # +---------------------------------+
        # | Texture particle size           |
        # | transformation:                 |
        # +---------------------------------+ 
        # Base (plot, geometry, reference) Clay Silt Sand 
        #   particle size limits (in micro-meters), in the form of 
        #   c(clay min,clay max = silt min, silt max = sand min, sand.max) 
        #   If plotted texture triangles or soil point data have 
        #   a different css.ps.lim, and if transformation is allowed, 
        #   then all Clay Silt Sand data will be coverted TO that 
        #   system.
        base.css.ps.lim = c(0,2,50,2000), 
        #
        # Same, but for the plotted texture triangle
        #   particle size limits (in micro-meters) (see above)
        tri.css.ps.lim  = NULL, 
        #
        # Same, but for the plotted soil point data (tri.data)
        #   particle size limits (in micro-meters) (see above)
        dat.css.ps.lim  = NULL, 
        # 
        # Should any texture triangle or soil point data, with 
        #   different Clay Silt Sand particle size limits from the 
        #   base plot, be transformed into the base plot system?
        css.transf      = FALSE, 
        #
        # Default function to be used when transforming soil 
        #   texture data or triangles (particle sizes) 
        #   alternatives functions should either have the same 
        #   options (even if some are not used), or a "..." option 
        #   that will be used as a "dump" for unused otions.
        text.transf.fun     = "TT.text.transf",  # CHARACTER vector
        #
        # In case an alternative to "TT.text.transf" is used in 
        #   text.transf.fun, it is possible to provide 2 additional 
        #   options to that new functions (free), unused in 
        #   "TT.text.transf".
        trsf.add.opt1       = NA,   # Additionnal option 1 
        trsf.add.opt2       = NA,   # Additionnal option 2 
        #
        # +-------------------------------------------------+
        # | HYPRES TEXTURE TRIANGLE -- ORIGINAL, WRONG NAME |
        # +-------------------------------------------------+
        #
        FAO50.TT  = list( # FAO TRIANGLE PARAMETERS :
            
            main            = "HYPRES (please use 'HYPRES.TT' instead)", 
            
            #                 The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  =-P01-   P02    P03    P04    P05    P06    P07    P08   -P09-   P10    P11   -P12-    
            "tt.points"     = data.frame( 
                "CLAY"      = c( 1.000, 0.600, 0.600, 0.350, 0.350, 0.350, 0.180, 0.180, 0.000, 0.000, 0.000, 0.000), 
                "SILT"      = c( 0.000, 0.000, 0.400, 0.000, 0.500, 0.650, 0.000, 0.170, 0.000, 0.350, 0.850, 1.000), 
                "SAND"      = c( 0.000, 0.400, 0.000, 0.650, 0.150, 0.000, 0.820, 0.650, 1.000, 0.650, 0.150, 0.000)  
            ),  #
            #
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "VF"        = list( "name" = "Very fine",       "points" = c(02,01,03)          ), 
                "F"         = list( "name" = "Fine",            "points" = c(04,02,03,06)       ), 
                "M"         = list( "name" = "Medium",          "points" = c(07,04,05,11,10,08) ), 
                "MF"        = list( "name" = "Medium fine",     "points" = c(11,05,06,12)       ), 
                "C"         = list( "name" = "Coarse",          "points" = c(09,07,08,10)       )  
            ),  #
            #
            # Traingle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = rep(T,3), 
            tlr.an          = c(60,60,60), 
            #
            blr.tx      = c("SAND","CLAY","SILT"), 
            # 
            base.css.ps.lim = c(0,2,50,2000), 
            tri.css.ps.lim  = c(0,2,50,2000), 
            #
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
            #
            # In fact it is the FAO soil texture classes: Info from SysCan
            # http://sis.agr.gc.ca/cansis/nsdb/lpdb/faotext.html
            # FAO Soil Texture
            # Texture is the relative proportion of sand, silt and clay of the dominant 
            # soil for each soil map polygon. Texture classes are:
            #
            # Coarse texture: sands, loamy sand and sandy loams with less than 18 % clay, 
            # and more than 65 % sand.
            #
            # Medium texture: sandy loams, loams, sandy clay loams, silt loams with less 
            # than 35 % clay and less than 65 % sand; the sand fractions may be as high as 82 % if a minimum of 18 % clay is present.
            #
            # Fine texture: clays, silty clays, sandy clays, clay loams and silty clay loams 
            # with more than 35 % clay.
            #
            # Where two or three texture names appear, this means that all named textures 
            # are present in the map unit.
            # 
            # Texture Codeset
            # COARSE
            # FINE
            # FINE-COARSE
            # FINE-MED-CRS
            # FINE-MEDIUM
            # MEDIUM
            # MEDIUM-COARSE
            #
        ),  #
        
        # +----------------------------------+
        # | HYPRES TEXTURE TRIANGLE, RENAMED |
        # +----------------------------------+
        
        HYPRES.TT  = list( # EU SOIL MAP TRIANGLE PARAMETERS :
            
            main            = "HYPRES / European Soil Map", 
            
            #                 The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  =-P01-   P02    P03    P04    P05    P06    P07    P08   -P09-   P10    P11   -P12-    
            "tt.points"     = data.frame( 
                "CLAY"      = c( 1.000, 0.600, 0.600, 0.350, 0.350, 0.350, 0.180, 0.180, 0.000, 0.000, 0.000, 0.000), 
                "SILT"      = c( 0.000, 0.000, 0.400, 0.000, 0.500, 0.650, 0.000, 0.170, 0.000, 0.350, 0.850, 1.000), 
                "SAND"      = c( 0.000, 0.400, 0.000, 0.650, 0.150, 0.000, 0.820, 0.650, 1.000, 0.650, 0.150, 0.000)  
            ),  #
            
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "VF"        = list( "name" = "Very fine",       "points" = c(02,01,03)          ), 
                "F"         = list( "name" = "Fine",            "points" = c(04,02,03,06)       ), 
                "M"         = list( "name" = "Medium",          "points" = c(07,04,05,11,10,08) ), 
                "MF"        = list( "name" = "Medium fine",     "points" = c(11,05,06,12)       ), 
                "C"         = list( "name" = "Coarse",          "points" = c(09,07,08,10)       )  
            ),  #
            
            # Traingle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = rep(T,3), 
            tlr.an          = c(60,60,60), 
            #
            blr.tx      = c("SAND","CLAY","SILT"), 
            # 
            base.css.ps.lim = c(0,2,50,2000), 
            tri.css.ps.lim  = c(0,2,50,2000), 
            #
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
            #
            # In fact it is the FAO soil texture classes: Info from SysCan
            # http://sis.agr.gc.ca/cansis/nsdb/lpdb/faotext.html
            # FAO Soil Texture
            # Texture is the relative proportion of sand, silt and clay of the dominant 
            # soil for each soil map polygon. Texture classes are:
            #
            # Coarse texture: sands, loamy sand and sandy loams with less than 18 % clay, 
            # and more than 65 % sand.
            #
            # Medium texture: sandy loams, loams, sandy clay loams, silt loams with less 
            # than 35 % clay and less than 65 % sand; the sand fractions may be as high as 82 % if a minimum of 18 % clay is present.
            #
            # Fine texture: clays, silty clays, sandy clays, clay loams and silty clay loams 
            # with more than 35 % clay.
            #
            # Where two or three texture names appear, this means that all named textures 
            # are present in the map unit.
            # 
            # Texture Codeset
            # COARSE
            # FINE
            # FINE-COARSE
            # FINE-MED-CRS
            # FINE-MEDIUM
            # MEDIUM
            # MEDIUM-COARSE
            #
        ),  #
        
        # +-----------------+
        # | OTHER TRIANGLES |
        # +-----------------+
        
        USDA.TT = list(  #  USDA Triangle parameters
            #
            main            = "USDA", 
            # 
            #                The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  = P01    P02    P03    P04    P05    P06    P07    P08    P09    P10    P11    P12
            #                  = P13    P14    P15    P16    P17    P18    P19    P20    P21    P22    P23   
            #                  = P24    P25    P26 (submits)
            "tt.points"     = data.frame( 
                "CLAY"      = c( 0.550, 0.600, 0.350, 0.350, 0.400, 0.400, 0.400, 0.200, 0.200, 0.275, 0.275, 0.275,  
                                 0.275, 0.150, 0.100, 0.075, 0.075, 0.125, 0.125, 0.000, 0.000, 0.000, 0.000,         
                                 1.000, 0.000, 0.000  ),  
                            #
                "SILT"      = c( 0.000, 0.400, 0.000, 0.200, 0.150, 0.400, 0.600, 0.000, 0.275, 0.275, 0.500, 0.525,  
                                 0.725, 0.000, 0.000, 0.400, 0.500, 0.800, 0.875, 0.150, 0.300, 0.500, 0.800,         
                                 0.000, 0.000, 1.000  ),  
                            #
                "SAND"      = c( 0.450, 0.000, 0.650, 0.450, 0.450, 0.200, 0.000, 0.800, 0.525, 0.450, 0.225, 0.200,  
                                 0.000, 0.850, 0.900, 0.525, 0.425, 0.075, 0.000, 0.850, 0.700, 0.500, 0.200,         
                                 0.000, 1.000, 0.000  )  
            ),  #
            #
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "Cl"        = list( "name" = "clay",            "points" = c(24,01,05,06,02)            ), 
                "SiCl"      = list( "name" = "silty clay",      "points" = c(02,06,07)                  ), 
                "SaCl"      = list( "name" = "sandy clay",      "points" = c(01,03,04,05)               ), 
                "ClLo"      = list( "name" = "clay loam",       "points" = c(05,04,10,11,12,06)         ), 
                "SiClLo"    = list( "name" = "silty clay loam", "points" = c(06,12,13,07)               ), 
                "SaClLo"    = list( "name" = "sandy clay loam", "points" = c(03,08,09,10,04)            ), 
                "Lo"        = list( "name" = "loam",            "points" = c(10,09,16,17,11)            ), 
                "SiLo"      = list( "name" = "silty loam",      "points" = c(11,17,22,23,18,19,13,12)   ), 
                "SaLo"      = list( "name" = "sandy loam",      "points" = c(08,14,21,22,17,16,09)      ), 
                "Si"        = list( "name" = "silt",            "points" = c(18,23,26,19)               ), 
                "LoSa"      = list( "name" = "loamy sand",      "points" = c(14,15,20,21)               ), 
                "Sa"        = list( "name" = "sand",            "points" = c(15,25,20)                  )  
            ),  #
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = rep(T,3), 
            tlr.an          = c(60,60,60), 
            #
            blr.tx      = c("SAND","CLAY","SILT"), 
            # 
            base.css.ps.lim = c(0,2,50,2000), 
            tri.css.ps.lim  = c(0,2,50,2000), 
            #
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  #
        # 
        FR.AISNE.TT = list( # AISNE/FRENCH TRIANGLE PARAMETERS :
            #
            main            = "Aisne (FR)", 
            #
            #                 The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  = P01    P02    P03    P04    P05    P06    P07    P08    P09    P10    P11    
            #                    P12    P13    P14    P15    P16    P17    P18   -P19-   P20    P21    P22    
            #                    P23    P24    P25    P26   -P27-  -P28-   P29    
            "tt.points"     = data.frame( 
                "CLAY"      = c( 0.450, 0.450, 0.450, 0.450, 0.250, 0.250, 0.300, 0.300, 0.300, 0.300, 0.100, 
                                 0.100, 0.125, 0.125, 0.175, 0.175, 0.175, 0.175, 0.000, 0.000, 0.075, 0.075, 
                                 0.075, 0.075, 0.000, 0.000, 0.000, 1.000, 0.300  ),   
                "SILT"      = c( 0.000, 0.100, 0.350, 0.550, 0.000, 0.200, 0.250, 0.350, 0.500, 0.700, 0.000, 
                                 0.100, 0.100, 0.325, 0.275, 0.475, 0.675, 0.825, 0.000, 0.200, 0.375, 0.575, 
                                 0.775, 0.925, 0.450, 0.850, 1.000, 0.000, 0.550  ),   
                "SAND"      = c( 0.550, 0.450, 0.200, 0.000, 0.750, 0.550, 0.450, 0.350, 0.200, 0.000, 0.900, 
                                 0.800, 0.775, 0.550, 0.550, 0.350, 0.150, 0.000, 1.000, 0.800, 0.550, 0.350, 
                                 0.150, 0.000, 0.550, 0.150, 0.000, 0.000, 0.150  )    
            ),  #
            # 
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "ALO"       = list( "name" = "Argile lourde",           "points" = c(28, 01, 04                 ) ), 
                "A"         = list( "name" = "Argile",                  "points" = c(02, 07, 09, 03             ) ), 
                "AL"        = list( "name" = "Argile limoneuse",        "points" = c(03, 09, 10, 04             ) ), 
                "AS"        = list( "name" = "Argile sableuse",         "points" = c(01, 05, 06, 07, 02         ) ), 
                "LA"        = list( "name" = "Limon argileux",          "points" = c(29, 17, 18, 10             ) ), 
                "LAS"       = list( "name" = "Limon argilo-sableux",    "points" = c(08, 16, 17, 29, 09         ) ), 
                "LSA"       = list( "name" = "Limon sablo-argileux",    "points" = c(07, 06, 15, 16, 08         ) ), 
                "SA"        = list( "name" = "Sable argileux",          "points" = c(05, 11, 12, 13, 14, 15, 06 ) ), 
                "LM"        = list( "name" = "Limon moyen",             "points" = c(17, 23, 24, 18             ) ), 
                "LMS"       = list( "name" = "Limon moyen sableux",     "points" = c(16, 22, 23, 17             ) ), 
                "LS"        = list( "name" = "Limon sableux",           "points" = c(15, 14, 21, 22, 16         ) ), 
                "SL"        = list( "name" = "Sable limoneux",          "points" = c(13, 12, 20, 25, 21, 14     ) ), 
                "S"         = list( "name" = "Sable",                   "points" = c(11, 19, 20, 12             ) ), 
                "LL"        = list( "name" = "Limon leger",             "points" = c(23, 26, 27, 24             ) ), 
                "LLS"       = list( "name" = "Limon leger sableux",     "points" = c(21, 25, 26, 23, 22         ) )  
                #
            ),  #
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = rep(T,3), 
            tlr.an          = c(60,60,60), 
            #
            blr.tx      = c("SAND","CLAY","SILT"), 
            # 
            base.css.ps.lim = c(0,2,50,2000), 
            tri.css.ps.lim  = c(0,2,50,2000), 
            # 
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  #
        # 
        FR.GEPPA.TT = list( # GEPPA/FRENCH TRIANGLE PARAMETERS :
            #
            main            = "GEPPA (FR)", 
            #
            #                 The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  P01          P02          P03          P04          P05          P06          P07          
            #                  P08          P09          P10          P11          P12          P13          P14          
            #                  P15          P16          P17          P18          P19          P20          P21          
            #                  P22          P23          P24          P25          P26          P27          P28          
            "tt.points"     = data.frame( 
                "CLAY"  =   c( 1.000000000, 0.600000000, 0.550000000, 0.450000000, 0.426000000, 0.394160600, 0.375000000, 
                               0.325000000, 0.307758600, 0.287600800, 0.275000000, 0.225000000, 0.209848500, 0.203698200, 
                               0.187764100, 0.175000000, 0.125000000, 0.111486500, 0.103378400, 0.087890630, 0.075000000, 
                               0.075000000, 0.033333330, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000),  
                "SILT"  =   c( 0.000000000, 0.000000000, 0.450000000, 0.000000000, 0.200000000, 0.465328500, 0.625000000, 
                               0.000000000, 0.250000000, 0.542288300, 0.725000000, 0.000000000, 0.250000000, 0.351479300, 
                               0.614392600, 0.825000000, 0.000000000, 0.250000000, 0.400000000, 0.686523440, 0.925000000, 
                               0.000000000, 0.250000000, 0.000000000, 0.250000000, 0.450000000, 0.750000000, 1.000000000),  
                "SAND"  =   c( 0.000000000, 0.400000000, 0.000000000, 0.550000000, 0.374000000, 0.140510900, 0.000000000, 
                               0.675000000, 0.442241400, 0.170110900, 0.000000000, 0.775000000, 0.540151500, 0.444822500, 
                               0.197843300, 0.000000000, 0.875000000, 0.638513500, 0.496621600, 0.225585930, 0.000000000, 
                               0.925000000, 0.716666670, 1.000000000, 0.750000000, 0.550000000, 0.250000000, 0.000000000)  
            ),  #
            # 
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "AA"    = list( "name" = "Argile lourde",           "points" = c(01, 02, 03             ) ),  
                "A"     = list( "name" = "Argileux",                "points" = c(02, 04, 05, 06, 07, 03 ) ),  
                "As"    = list( "name" = "Argile sableuse",         "points" = c(04, 08, 09, 05         ) ),  
                "Als"   = list( "name" = "Argile limono-sableuse",  "points" = c(05, 09, 10, 06         ) ),  
                "Al"    = list( "name" = "Argile limoneuse",        "points" = c(06, 10, 11, 07         ) ),  
                "AS"    = list( "name" = "Argilo-sableux",          "points" = c(08, 12, 13, 09         ) ),  
                "LAS"   = list( "name" = "Limon argilo-sableux",    "points" = c(09, 13, 14, 15, 10     ) ),  
                "La"    = list( "name" = "Limon argileux",          "points" = c(10, 15, 16, 11         ) ),  
                "Sa"    = list( "name" = "Sable argileux",          "points" = c(12, 17, 18, 13         ) ),  
                "Sal"   = list( "name" = "Sable argilo-limoneux",   "points" = c(13, 18, 19, 14         ) ),  
                "Lsa"   = list( "name" = "Limon sablo-argileux",    "points" = c(14, 19, 20, 15         ) ),  
                "L"     = list( "name" = "Limon",                   "points" = c(15, 20, 21, 16         ) ),  
                "S"     = list( "name" = "Sableux",                 "points" = c(17, 22, 23, 18         ) ),  
                "SS"    = list( "name" = "Sable",                   "points" = c(22, 24, 25, 23         ) ),  
                "Sl"    = list( "name" = "Sable limoneux",          "points" = c(18, 23, 25, 26, 19     ) ),  
                "Ls"    = list( "name" = "Limon sableux",           "points" = c(19, 26, 27, 20         ) ),  
                "LL"    = list( "name" = "Limon pur",               "points" = c(20, 27, 28, 21         ) )   
                #
            ),  #
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = c(F,T,NA), 
            tlr.an          = c(45,90,45), 
            #
            blr.tx      = c("SILT","CLAY","SAND"), 
            # 
            base.css.ps.lim = c(0,2,50,2000), 
            tri.css.ps.lim  = c(0,2,50,2000), 
            #
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  #
        #
        DE.BK94.TT  = list( # GERMAN TRIANGLE PARAMETERS :
            #
            main            = "Bodenkundliche Kartieranleitung 1994 (DE)", 
            #
            #                 The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  P01    P02    P03    P04    P05    P06    P07    P08    P09    P10    P11    P12    P13    
            #                  P14    P15    P16    P17    P18    P19    P20    P21    P22    P23    P24    P25    P26    
            #                  P27    P28    P29    P30    P31    P32    P33    P34    P35    P36    P37    P38    P39    
            #                  P40    P41    P42    P43    P44    P45    P46    P47    P48    P49    P50    P51    P52    
            #                  P53    
            "tt.points"     = data.frame( 
                "CLAY"  =   c( 0.000, 0.080, 0.120, 0.170, 0.000, 0.080, 0.250, 0.080, 0.120, 0.170, 0.250, 0.300, 0.350, 
                               0.450, 0.000, 0.080, 0.170, 0.250, 0.300, 0.350, 0.450, 0.000, 0.080, 0.120, 0.170, 0.250, 
                               0.650, 0.170, 0.250, 0.350, 0.450, 0.650, 0.000, 0.050, 0.080, 0.170, 0.250, 0.350, 0.450, 
                               0.650, 0.000, 0.050, 0.080, 0.120, 0.170, 0.000, 0.050, 0.170, 0.250, 0.350, 0.450, 0.650, 
                               1.000),  
                "SILT"  =   c( 1.000, 0.920, 0.880, 0.830, 0.800, 0.800, 0.750, 0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 
                               0.550, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.400, 0.400, 0.400, 0.400, 0.400, 
                               0.350, 0.300, 0.300, 0.300, 0.300, 0.300, 0.250, 0.250, 0.250, 0.150, 0.150, 0.150, 0.150, 
                               0.150, 0.100, 0.100, 0.100, 0.100, 0.100, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
                               0.000),  
                "SAND"  =   c( 0.000, 0.000, 0.000, 0.000, 0.200, 0.120, 0.000, 0.270, 0.230, 0.180, 0.100, 0.050, 0.000, 
                               0.000, 0.500, 0.420, 0.330, 0.250, 0.200, 0.150, 0.050, 0.600, 0.520, 0.480, 0.430, 0.350, 
                               0.000, 0.530, 0.450, 0.350, 0.250, 0.050, 0.750, 0.700, 0.670, 0.680, 0.600, 0.500, 0.400, 
                               0.200, 0.900, 0.850, 0.820, 0.780, 0.730, 1.000, 0.950, 0.830, 0.750, 0.650, 0.550, 0.350, 
                               0.000)  
            ),  #
            # 
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "Ss"    = list( "name" = "reiner Sand",             "points" = c(41, 46, 47, 42         ) ), 
                "Su2"   = list( "name" = "Schwach schluffiger Sand","points" = c(33, 41, 42, 34         ) ), 
                "Sl2"   = list( "name" = "Schwach lehmiger Sand",   "points" = c(34, 42, 43, 35         ) ), 
                "Sl3"   = list( "name" = "Mittel lehmiger Sand",    "points" = c(23, 35, 43, 44, 24     ) ), 
                "St2"   = list( "name" = "Schwach toniger Sand",    "points" = c(42, 47, 48, 45, 44, 43 ) ), 
                "Su3"   = list( "name" = "Mittel schluffiger Sand", "points" = c(22, 33, 34, 35, 23     ) ), 
                "Su4"   = list( "name" = "Stark schluffiger Sand",  "points" = c(15, 22, 23, 16         ) ), 
                "Slu"   = list( "name" = "Schluffig-lehmiger Sand", "points" = c(16, 23, 24, 25, 17     ) ), 
                "Sl4"   = list( "name" = "Stark lehmiger Sand",     "points" = c(24, 44, 45, 36, 28, 25 ) ), 
                "St3"   = list( "name" = "Mittel toniger Sand",     "points" = c(36, 45, 48, 49, 37     ) ), 
                "Ls2"   = list( "name" = "Schwach sandiger Lehm",   "points" = c(17, 25, 26, 18         ) ), 
                "Ls3"   = list( "name" = "Mittel sandiger Lehm",    "points" = c(25, 28, 29, 26         ) ), 
                "Ls4"   = list( "name" = "Stark sandiger Lehm",     "points" = c(28, 36, 37, 29         ) ), 
                "Lt2"   = list( "name" = "Schwach toniger Lehm",    "points" = c(18, 26, 29, 30, 20, 19 ) ), 
                "Lts"   = list( "name" = "Sandig-toniger Lehm",     "points" = c(29, 37, 38, 39, 31, 30 ) ), 
                "Ts4"   = list( "name" = "Stark sandiger Ton",      "points" = c(37, 49, 50, 38         ) ), 
                "Ts3"   = list( "name" = "Mittel sandiger Ton",     "points" = c(38, 50, 51, 39         ) ), 
                "Uu"    = list( "name" = "Reiner Schluff",          "points" = c(01, 05, 06, 02         ) ), 
                "Us"    = list( "name" = "Sandiger Schluff",        "points" = c(05, 15, 16, 08, 06     ) ), 
                "Ut2"   = list( "name" = "Schwach toniger Schluff", "points" = c(02, 06, 08, 09, 03     ) ), 
                "Ut3"   = list( "name" = "Mittel toniger Schluff",  "points" = c(03, 09, 10, 04         ) ), 
                "Uls"   = list( "name" = "Sandig-lehmiger Schluff", "points" = c(08, 16, 17, 10, 09     ) ), 
                "Ut4"   = list( "name" = "Stark toniger Schluff",   "points" = c(04, 10, 11, 07         ) ), 
                "Lu"    = list( "name" = "Schluffiger Lehm",        "points" = c(10, 17, 18, 19, 12, 11 ) ), 
                "Lt3"   = list( "name" = "Mittel toniger Lehm",     "points" = c(20, 30, 31, 21         ) ), 
                "Tu3"   = list( "name" = "Mittel schluffiger Ton",  "points" = c(12, 19, 20, 21, 14, 13 ) ), 
                "Tu4"   = list( "name" = "Stark schluffiger Ton",   "points" = c(07, 11, 12, 13         ) ), 
                "Ts2"   = list( "name" = "Schwach sandiger Ton",    "points" = c(39, 51, 52, 40         ) ), 
                "Tl"    = list( "name" = "Lehmiger Ton",            "points" = c(31, 39, 40, 32         ) ), 
                "Tu2"   = list( "name" = "Schwach schluffiger Ton", "points" = c(14, 21, 31, 32, 27     ) ), 
                "Tt"    = list( "name" = "Reiner Ton",              "points" = c(27, 32, 40, 52, 53     ) ) 
                #
            ),  #
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = c(F,T,NA), 
            tlr.an          = c(45,90,45), 
            #
            blr.tx      = c("CLAY","SILT","SAND"), 
            # 
            base.css.ps.lim = c(0,2,63,2000), 
            tri.css.ps.lim  = c(0,2,63,2000), 
            # 
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  #
        #
        UK.SSEW.TT  = list( # SSEW-DEFRA TRIANGLE PARAMETERS :
            #
            main            = "Soil Survey of England and Wales (UK)", 
            # 
            #                The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  = P01   P02   P03   P04   P05   P06   P07   P08   P09   P10   P11   
            #                    P12   P13   P14   P15   P16   P17   P18   P19   P20
            "tt.points"     = data.frame( 
                "CLAY"      = c( 1.00, 0.55, 0.55, 0.35, 0.35, 0.35, 0.30, 0.30, 0.18, 0.18, 0.18, 
                                 0.18, 0.15, 0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 ),  
                            #                 
                "SILT"      = c( 0.00, 0.00, 0.45, 0.20, 0.45, 0.65, 0.00, 0.20, 0.00, 0.32, 0.62, 
                                 0.82, 0.00, 0.00, 0.00, 0.15, 0.30, 0.50, 0.80, 1.00 ),  
                            #                 
                "SAND"      = c( 0.00, 0.45, 0.00, 0.45, 0.20, 0.00, 0.70, 0.50, 0.82, 0.50, 0.20, 
                                 0.00, 0.85, 0.90, 1.00, 0.85, 0.70, 0.50, 0.20, 0.00 )  
            ),  #
            #
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "Cl"        = list( "name" = "Clay",            "points" = c(01,02,04,05,03)    ), 
                "SaCl"      = list( "name" = "Sandy clay",      "points" = c(02,07,08,04)       ), 
                "SiCl"      = list( "name" = "Silty clay",      "points" = c(03,05,06)          ), 
                "ClLo"      = list( "name" = "Clay loam",       "points" = c(04,08,10,11,05)    ), 
                "SiClLo"    = list( "name" = "Silty clay loam", "points" = c(05,11,12,06)       ), 
                "SaClLo"    = list( "name" = "Sandy clay loam", "points" = c(07,09,10,08)       ), 
                "SaLo"      = list( "name" = "Sandy loam",      "points" = c(09,13,17,18,10)    ), 
                "SaSiLo"    = list( "name" = "Sandy silt loam", "points" = c(10,18,19,11)       ), 
                "SiLo"      = list( "name" = "Silt loam",       "points" = c(11,19,20,12)       ), 
                "LoSa"      = list( "name" = "Loamy sand",      "points" = c(13,14,16,17)       ), 
                "Sa"        = list( "name" = "Sand",            "points" = c(14,15,16)          )  
            ),  #
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = rep(T,3), 
            tlr.an          = c(60,60,60), 
            #
            blr.tx      = c("SAND","CLAY","SILT"), 
            # 
            base.css.ps.lim = c(0,2,60,2000), 
            tri.css.ps.lim  = c(0,2,60,2000), 
            # 
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  #
        # 
        "AU.TT" = list( # Australian TRIANGLE PARAMETERS :
            #
            main            = "Autralia (AU)", 
            #
            #                 The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  P01  P02  P03  P04  P05  P06  P07  
            #                  P08  P09  P10  P11  P12  P13  P14  
            #                  P15  P16  P17  P18  P19  P20  P21  
            #                  P22  
            "tt.points"     = data.frame( 
                "CLAY"  =   c( 1.00,  0.74,  0.40,  0.26,  0.00,  0.40, 
                               0.26,  0.50,  0.12,  0.30,  0.205, 0.00, 
                               0.265, 0.105, 0.17,  0.09,  0.00,  0.085, 
                               0.00 ),  
                "SILT"  =   c( 0.00,  0.26,  0.60,  0.74,  1.00,  0.26, 
                               0.26,  0.00,  0.26,  0.07,  0.10,  0.26, 
                               0.00,  0.1325,0.00,  0.04,  0.08,  0.00, 
                               0.00 ), 
                "SAND"  =   c( 0.00,  0.00,  0.00,  0.00,  0.00,  0.34, 
                               0.48,  0.50,  0.62,  0.63,  0.695, 0.74, 
                               0.735, 0.7625,0.83,  0.87,  0.92, 0.915, 
                               1.00 )   
            ),  #
            # 
            #
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "Cl"        = list( "name" = "Clay",            "points" = c(01, 02, 06, 10, 08) ),  
                "SiCl"      = list( "name" = "Silty clay",      "points" = c(02, 03, 06) ),  
                "SiClLo"    = list( "name" = "Silt clay loam",  "points" = c(03, 04, 07, 06) ),  
                "SiLo"      = list( "name" = "Silty loam",      "points" = c(04, 05, 12, 09, 07) ),  
                "ClLo"      = list( "name" = "Clay loam",       "points" = c(06, 07, 11, 10) ),  
                "Lo"        = list( "name" = "Loam",            "points" = c(07, 09, 14, 11) ),  
                "LoSa"      = list( "name" = "Loamy sand",      "points" = c(09, 12, 17, 19, 18, 16, 14) ),  
                "SaCl"      = list( "name" = "Sandy clay",      "points" = c(08, 10, 13) ),  
                "SaClLo"    = list( "name" = "Sandy clay loam", "points" = c(10, 11, 15, 13) ),  
                "SaLo"      = list( "name" = "Sandy loam",      "points" = c(11, 14, 16, 18, 15) ),  
                "Sa"        = list( "name" = "Sand",            "points" = c(16, 17, 19, 18) )   
            ),  #
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = c(F,T,NA), 
            tlr.an          = c(45,90,45), 
            #
            blr.tx      = c("SAND","CLAY","SILT"), 
            # 
            base.css.ps.lim = c(0,2,20,2000), 
            tri.css.ps.lim  = c(0,2,20,2000), 
            #
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  #
        # 
        "AU2.TT" = list( # Australian TRIANGLE PARAMETERS :
            main            = "Autralia (AU)", 
            #
            # The coordinates of this triangle were kindly provided by Budiman Minasni
            #    as a replacement for the original implementation of the Australian 
            #    texture triangle.
            #
            #                 The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  P01    P02    P03    P04    P05    P06    P07  
            #                  P08    P09    P10    P11    P12    P13    P14  
            #                  P15    P16    P17    P18    P19  
            "tt.points"     = data.frame( 
                "CLAY"  =   c( 0.510, 0.300, 0.400, 0.750, 1.000, 0.400, 0.260, 
                               0.260, 0.000, 0.000, 0.120, 0.085, 0.000, 0.000, 
                               0.080, 0.100, 0.210, 0.170, 0.260 ), 
                "SILT"  =   c( 0.000, 0.070, 0.250, 0.250, 0.000, 0.600, 0.740, 
                               0.250, 0.250, 1.000, 0.250, 0.040, 0.075, 0.000, 
                               0.000, 0.130, 0.100, 0.000, 0.000 ), 
                "SAND"  =   c( 0.490, 0.630, 0.350, 0.000, 0.000, 0.000, 0.000, 
                               0.490, 0.750, 0.000, 0.630, 0.875, 0.925, 1.000, 
                               0.920, 0.770, 0.690, 0.830, 0.740 ) 
            ),  #
            # 
            #
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "C"   = list( "name" = "Clay",            "points" = c(01, 02, 03, 04, 05, 01) ),  
                "ZC"  = list( "name" = "Silty clay",      "points" = c(03, 06, 04, 03) ),  
                "ZCL" = list( "name" = "Silty clay loam", "points" = c(08, 07, 06, 03, 08) ),  
                "ZL"  = list( "name" = "Silty loam",      "points" = c(09, 10, 07, 08, 11, 09) ),  
                "LS"  = list( "name" = "Loamy sand",      "points" = c(13, 09, 11, 12) ),  
                "S"   = list( "name" = "Sand",            "points" = c(14, 13, 12, 15, 14) ),  
                "SL"  = list( "name" = "Sandy loam",      "points" = c(15, 16, 17, 18, 15) ),  
                "L"   = list( "name" = "Loam",            "points" = c(16, 11, 08, 17, 16) ),  
                "SCL" = list( "name" = "Sandy clay loam", "points" = c(18, 17, 02, 19, 18) ),  
                "CL"  = list( "name" = "Clay loam",       "points" = c(17, 08, 03, 02) ),  
                "SC"  = list( "name" = "Sandy clay",      "points" = c(19, 02, 01, 19) )   
            ),  #
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = c(F,T,NA), 
            tlr.an          = c(45,90,45), 
            #
            blr.tx      = c("SAND","CLAY","SILT"), 
            # 
            base.css.ps.lim = c(0,2,20,2000), 
            tri.css.ps.lim  = c(0,2,20,2000), 
            #
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  #
        #
        "BE.TT" = list( # Belgian TRIANGLE PARAMETERS :
            #
            main            = "Belgium (BE)", 
            #
            #                 The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  =-P01-   P02    P03    P04    P05    P06    P07    P08      P09     P10    P11    
            #                    P12    P13    P14    P15   -P16-   P17    P18    P19      P20   -P21-   
            "tt.points"     = data.frame( 
                "CLAY"      = c( 1.000, 0.450, 0.350, 0.350, 0.175, 0.175, 0.175, 0.18125, 0.200, 0.225, 0.300,   
                                 0.0875,0.0875,0.1125,0.1125,0.000, 0.000, 0.000, 0.000,   0.000, 0.000 ),   
                "SILT"      = c( 0.000, 0.550, 0.000, 0.550, 0.000, 0.150, 0.525, 0.56875, 0.600, 0.625, 0.700,   
                                 0.000, 0.0875,0.2125,0.3875,0.000, 0.175, 0.325, 0.500,   0.850, 1.000 ),   
                "SAND"      = c( 0.000, 0.000, 0.650, 0.100, 0.825, 0.675, 0.300, 0.250,   0.200, 0.150, 0.000,   
                                 0.9125,0.825, 0.675, 0.500, 1.000, 0.825, 0.675, 0.500,   0.150, 0.000 )    
            ),  #
            # 
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "U" = list( "name" = "Argile lourde | Zware klei",              "points" = c(01, 03, 04, 02                     ) ), 
                "E" = list( "name" = "Argile | Klei",                           "points" = c(03, 04, 02, 11, 10, 09, 08, 07, 05 ) ), 
                "A" = list( "name" = "Limon | Leem",                            "points" = c(10, 11, 21, 20                     ) ), 
                "L" = list( "name" = "Limon sableux | Zandleem",                "points" = c(06, 07, 08, 09, 10, 20, 19, 15, 14 ) ), 
                "P" = list( "name" = "Limon sableux leger | Licht zandleem",    "points" = c(14, 15, 19, 18                     ) ), 
                "S" = list( "name" = "Sable limoneux | Lemig zand",             "points" = c(05, 06, 14, 18, 17, 13, 12         ) ), 
                "Z" = list( "name" = "Sable | Zand",                            "points" = c(12, 13, 17, 16                     ) )  
                #
            ),  #
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = rep(F,3), 
            tlr.an          = c(60,60,60), 
            #
            blr.tx      = c("SILT","SAND","CLAY"), 
            # 
            base.css.ps.lim = c(0,2,50,2000), 
            tri.css.ps.lim  = c(0,2,50,2000), 
            # 
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  #
        # 
        "CA.FR.TT" = list( # Canadian TRIANGLE PARAMETERS : (added 2010-04-16) 
            #
            main            = "Canada (CA)", 
            #
            #                 The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  =-P01-   P02    P03    P04    P05    P06    P07    P08    P09    P10    P11    
            #                    P12    P13    P14    P15    P16    P17    P18    P19    P20    P21   -P22-   
            #                    P23    P24    P25    P26   -P27-   
            "tt.points"     = data.frame( 
                "CLAY"      = c( 1.000, 0.600, 0.600, 0.550, 0.400, 0.400, 0.400, 0.350, 0.350, 0.270, 0.270, 
                                 0.270, 0.270, 0.200, 0.200, 0.120, 0.120, 0.150, 0.100, 0.070, 0.070, 0.000, 
                                 0.000, 0.000, 0.000, 0.000, 0.000 ),   
                "SILT"      = c( 0.000, 0.400, 0.000, 0.000, 0.600, 0.400, 0.150, 0.200, 0.000, 0.730, 0.530, 
                                 0.500, 0.280, 0.280, 0.000, 0.880, 0.800, 0.000, 0.000, 0.500, 0.410, 1.000, 
                                 0.800, 0.500, 0.300, 0.150, 0.000 ),   
                "SAND"      = c( 0.000, 0.000, 0.400, 0.450, 0.000, 0.200, 0.450, 0.450, 0.650, 0.000, 0.200, 
                                 0.230, 0.450, 0.520, 0.800, 0.000, 0.080, 0.850, 0.900, 0.430, 0.520, 0.000, 
                                 0.200, 0.500, 0.700, 0.850, 1.000 )    
                #                http://sis.agr.gc.ca/cansis/glossary/texture,_soil.html
            ),  #
            # 
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "ALo"  = list( "name" = "Argile lourde",         "points" = c( 01, 02, 03 ) ), 
                "ALi"  = list( "name" = "Argile limoneuse",      "points" = c( 02, 05, 06 ) ), 
                "A"    = list( "name" = "Argile",                "points" = c( 02, 03, 04, 07, 06 ) ), 
                "AS"   = list( "name" = "Argile sableuse",       "points" = c( 04, 07, 08, 09 ) ), 
                "LLiA" = list( "name" = "Loam limono-argileux",  "points" = c( 05, 06, 11, 10 ) ), 
                "LA"   = list( "name" = "Loam argileux",         "points" = c( 06, 07, 08, 13, 12, 11 ) ), 
                "LSA"  = list( "name" = "Loam sablo-argileux",   "points" = c( 08, 09, 15, 14, 13 ) ), 
                "LLi"  = list( "name" = "Loam limoneux",         "points" = c( 10, 11, 12, 20, 24, 23, 17, 16 ) ), 
                "L"    = list( "name" = "Loam",                  "points" = c( 12, 13, 14, 21, 20 ) ), 
                "LS"   = list( "name" = "Loam sableux",          "points" = c( 14, 15, 18, 25, 24, 20, 21 ) ), 
                "SL"   = list( "name" = "Sable loameux",         "points" = c( 18, 19, 26, 25 ) ), 
                "Li"   = list( "name" = "Limon",                 "points" = c( 16, 17, 23, 22 ) ), 
                "S"    = list( "name" = "Sable",                 "points" = c( 19, 27, 26 ) )  
                #
            ),  #
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = c(F,T,NA), 
            tlr.an          = c(45,90,45), 
            #
            blr.tx      = c("SAND","CLAY","SILT"), 
            # 
            base.css.ps.lim = c(0,2,50,2000), # http://sis.agr.gc.ca/cansis/glossary/separates,_soil.html
            tri.css.ps.lim  = c(0,2,50,2000), 
            # 
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  #
        #
        "CA.EN.TT" = list( # Canadian TRIANGLE PARAMETERS : (added 2010-04-16) 
            #
            main            = "Canada (CA)", 
            #
            #                 The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  =-P01-   P02    P03    P04    P05    P06    P07    P08    P09    P10    P11    
            #                    P12    P13    P14    P15    P16    P17    P18    P19    P20    P21   -P22-   
            #                    P23    P24    P25    P26   -P27-   
            "tt.points"     = data.frame( 
                "CLAY"      = c( 1.000, 0.600, 0.600, 0.550, 0.400, 0.400, 0.400, 0.350, 0.350, 0.270, 0.270, 
                                 0.270, 0.270, 0.200, 0.200, 0.120, 0.120, 0.150, 0.100, 0.070, 0.070, 0.000, 
                                 0.000, 0.000, 0.000, 0.000, 0.000 ),   
                "SILT"      = c( 0.000, 0.400, 0.000, 0.000, 0.600, 0.400, 0.150, 0.200, 0.000, 0.730, 0.530, 
                                 0.500, 0.280, 0.280, 0.000, 0.880, 0.800, 0.000, 0.000, 0.500, 0.410, 1.000, 
                                 0.800, 0.500, 0.300, 0.150, 0.000 ),   
                "SAND"      = c( 0.000, 0.000, 0.400, 0.450, 0.000, 0.200, 0.450, 0.450, 0.650, 0.000, 0.200, 
                                 0.230, 0.450, 0.520, 0.800, 0.000, 0.080, 0.850, 0.900, 0.430, 0.520, 0.000, 
                                 0.200, 0.500, 0.700, 0.850, 1.000 )    
                #                http://sis.agr.gc.ca/cansis/glossary/texture,_soil.html
            ),  #
            # 
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "HCl"    = list( "name" = "Heavy clay",       "points" = c( 01, 02, 03 ) ), 
                "SiCl"   = list( "name" = "Silty clay",       "points" = c( 02, 05, 06 ) ), 
                "Cl"     = list( "name" = "Clay",             "points" = c( 02, 03, 04, 07, 06 ) ), 
                "SaCl"   = list( "name" = "Sandy clay",       "points" = c( 04, 07, 08, 09 ) ), 
                "SiClLo" = list( "name" = "Silty clay loam",  "points" = c( 05, 06, 11, 10 ) ), 
                "ClLo"   = list( "name" = "Clay loam",        "points" = c( 06, 07, 08, 13, 12, 11 ) ), 
                "SaClLo" = list( "name" = "Sandy clay loam",  "points" = c( 08, 09, 15, 14, 13 ) ), 
                "SiLo"   = list( "name" = "Silty loam",       "points" = c( 10, 11, 12, 20, 24, 23, 17, 16 ) ), 
                "L"      = list( "name" = "Loam",             "points" = c( 12, 13, 14, 21, 20 ) ), 
                "SaLo"   = list( "name" = "Sandy loam",       "points" = c( 14, 15, 18, 25, 24, 20, 21 ) ), 
                "LoSa"   = list( "name" = "Loamy sand",       "points" = c( 18, 19, 26, 25 ) ), 
                "Si"     = list( "name" = "Silt",             "points" = c( 16, 17, 23, 22 ) ), 
                "Sa"     = list( "name" = "Sand",             "points" = c( 19, 27, 26 ) )  
                #
            ),  #
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = c(F,T,NA), 
            tlr.an          = c(45,90,45), 
            #
            blr.tx      = c("SAND","CLAY","SILT"), 
            # 
            base.css.ps.lim = c(0,2,50,2000), # http://sis.agr.gc.ca/cansis/glossary/separates,_soil.html
            tri.css.ps.lim  = c(0,2,50,2000), 
            # 
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  #
        #
        "ISSS.TT" = list(  #  ISSS Triangle parameters
            #                 by Wei Shangguan, School of geography, Beijing normal university 
            #                 and after Verheye, W., and J. Ameryckx. 1984. Mineral fractions 
            #                 and classificaton of soil texture. Pedologie, 2, 215-225.
            main            = "ISSS", 
            # 
            #                The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  = P01    P02    P03    P04    P05    P06    P07    P08    P09    P10    P11    P12
            #                  = P13    P14    P15    P16    P17    P18      
            "tt.points"     = data.frame( 
                "CLAY"      = c( 1.000, 0.450, 0.450, 0.450, 0.250, 0.250, 0.250, 0.250, 0.150, 0.150, 0.150, 0.150,  
                                 0.050, 0.000, 0.000, 0.000, 0.000, 0.000),  
                            #
                "SILT"      = c( 0.000, 0.000, 0.450, 0.550, 0.000, 0.200, 0.450, 0.750, 0.000, 0.200, 0.450, 0.850,  
                                 0.000, 0.000, 0.150, 0.350, 0.450, 1.000),  
                            #
                "SAND"      = c( 0.000, 0.550, 0.100, 0.000, 0.750, 0.550, 0.300, 0.000, 0.850, 0.650, 0.400, 0.000,  
                                 0.950, 1.000, 0.850, 0.650, 0.550, 0.000)  
            ),  #
            #
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "HCl"        = list( "name" = "heavy clay",      "points" = c(01,02,03,04)        ), 
                "SaCl"       = list( "name" = "sandy clay",      "points" = c(02,05,06)           ), 
                "LCl"        = list( "name" = "light clay",      "points" = c(02,06,07,03)        ), 
                "SiCl"       = list( "name" = "silty clay",      "points" = c(03,07,08,04)        ), 
                "SaClLo"     = list( "name" = "sandy clay loam", "points" = c(05,09,10,06)        ), 
                "ClLo"       = list( "name" = "clay loam",       "points" = c(06,10,11,07)        ), 
                "SiClLo"     = list( "name" = "silty clay loam", "points" = c(07,11,12,08)        ), 
                "LoSa"       = list( "name" = "loamy sand",      "points" = c(09,13,15)           ), 
                "Sa"         = list( "name" = "sand",            "points" = c(13,14,15)           ), 
                "SaLo"       = list( "name" = "sandy loam",      "points" = c(09,15,16,10)        ), 
                "Lo"         = list( "name" = "loam",            "points" = c(10,16,17,11)        ), 
                "SiLo"       = list( "name" = "silt loam",       "points" = c(11,17,18,12)        )  
            ),  #
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = rep(T,3), 
            tlr.an          = c(60,60,60), 
            #
            blr.tx      = c("SAND","CLAY","SILT"), 
            # 
            base.css.ps.lim = c(0,2,20,2000), 
            tri.css.ps.lim  = c(0,2,20,2000), 
            #
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  #
        #
        "ROM.TT" = list(# ROM TRIANGLE PARAMETERS: Added 2010/06/07 
            #                  by Rosca Bogdan, Romanian Academy 
            #                  Iasi Branch, Geography team
            #
            main = "SRTS 2003", 
            # 
            #                The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  = P01    P02    P03    P04    P05    P06    P07    P08    P09    P10    P11    P12
            #                  = P13    P14    P15    P16    P17    P18    P19    P20    P21    P22    P23   
            #                  = P24    P25    P26    P27 (submits)
            "tt.points"     = data.frame( 
                "CLAY"      = c( 0.700, 0.700, 0.600, 0.600, 0.600, 0.450, 0.450, 0.450, 0.450, 0.330, 0.330, 0.330,  
                                 0.330, 0.200, 0.200, 0.200, 0.200, 0.200, 0.120, 0.120, 0.050, 0.050, 0.000,         
                                 0.000, 1.000, 0.000, 0.000 ),  
                            #
                "SILT"      = c( 0.000, 0.300, 0.000, 0.330, 0.400, 0.000, 0.140, 0.330, 0.550, 0.000, 0.140, 0.330,  
                                 0.670, 0.000, 0.140, 0.330, 0.500, 0.800, 0.000, 0.330, 0.000, 0.330, 0.330,         
                                 0.500, 0.000, 0.000, 1.000 ),  
                            #
                "SAND"      = c( 0.300, 0.000, 0.400, 0.070, 0.000, 0.550, 0.410, 0.220, 0.000, 0.670, 0.530, 0.340,  
                                 0.000, 0.800, 0.660, 0.470, 0.300, 0.000, 0.880, 0.550, 0.950, 0.620, 0.670,         
                                 0.500, 0.000, 1.000, 0.000 )  
            ),  #
            #
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "AF"        = list( "name" = "argila fina",          "points" = c(01,25,02                 ) ), 
                "AA"        = list( "name" = "argila medie",         "points" = c(01,03,04,05,02           ) ), 
                "AP"        = list( "name" = "argila prafoasa",      "points" = c(04,08,09,05              ) ), 
                "AL"        = list( "name" = "argila lutoasa",       "points" = c(03,06,07,08,04           ) ), 
                "TP"        = list( "name" = "lut argilo-prafos",    "points" = c(08,12,13,09              ) ), 
                "TT"        = list( "name" = "lut argilos mediu",    "points" = c(07,11,12,08              ) ), 
                "TN"        = list( "name" = "argila nisipoasa",     "points" = c(06,10,11,07              ) ), 
                "LP"        = list( "name" = "lut prafos",           "points" = c(12,16,17,18,13           ) ), 
                "LL"        = list( "name" = "lut mediu",            "points" = c(11,15,16,12              ) ), 
                "LN"        = list( "name" = "lut nisipo-argilos",   "points" = c(10,14,15,11              ) ), 
                "SP"        = list( "name" = "praf",                 "points" = c(17,24,27,18              ) ), 
                "SS"        = list( "name" = "lut nisipos prafos",   "points" = c(16,20,22,23,24,17        ) ),
                "SG+SM+SF"  = list( "name" = "lut nisipos",          "points" = c(14,19,20,16,15           ) ),
                "UG+UM+UF"  = list( "name" = "nisip lutos",          "points" = c(19,21,22,20              ) ),
                "NG+NM+NF"  = list( "name" = "nisip",                "points" = c(21,26,23,22              ) )            
            ),  #
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = rep(T,3), 
            tlr.an          = c(60,60,60), 
            #
            blr.tx      = c("SAND","CLAY","SILT"), 
            # 
            base.css.ps.lim = c(0,2,20,2000), 
            tri.css.ps.lim  = c(0,2,20,2000), 
            #
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  #
        #
        "DE.SEA74.TT"  = list( # GDR Forest soils TRIANGLE PARAMETERS :
            #
            main            = "Standortserkundungsanweisung SEA 1974 (DE)", 
            #
            #                 The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes. The clay definition of points 08, 09, 10 closely follows the
            #                   triangle, plotted in SEA 1974, conforms to the GDR standard (TGL 24300-05:1985-06), but
            #                   is 1 Percent larger than the actual German texture triangle in DE.BK94.TT :
            #                  P01    P02    P03    P04    P05    P06    P07    P08    P09    P10    P11    P12    P13    
            #                  P14    P15    P16    P17    P18    P19    P20    P21    P22    P23    P24    P25    P26    
            "tt.points"     = data.frame( 
                "CLAY"  =   c( 1.000, 0.450, 0.450, 0.300, 0.300, 0.300, 0.300, 0.180, 0.180, 0.180, 0.150, 0.120, 0.050, 
                               0.000, 0.000, 0.000, 0.000, 0.050, 0.100, 0.025, 0.050, 0.080, 0.080, 0.000, 0.000, 0.450 ),
                "SILT"  =   c( 0.000, 0.550, 0.000, 0.700, 0.550, 0.150, 0.000, 0.820, 0.550, 0.150, 0.000, 0.150, 0.550, 
                               0.550, 0.200, 0.100, 0.000, 0.000, 0.000, 0.075, 0.150, 0.800, 0.920, 1.000, 0.800, 0.400 ),  
                "SAND"  =   c( 0.000, 0.000, 0.550, 0.000, 0.150, 0.550, 0.700, 0.000, 0.270, 0.670, 0.850, 0.730, 0.400, 
                               0.450, 0.800, 0.900, 1.000, 0.950, 0.900, 0.900, 0.800, 0.120, 0.000, 0.000, 0.200, 0.150 )  
            ),  #
            # 
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "L"   = list( "name" = "Lehm",                "points" = c(10, 06, 05, 09 ) ), 
                "stL" = list( "name" = "sandig-toniger Lehm", "points" = c(11, 07, 06, 10, 12 ) ), 
                "sL"  = list( "name" = "sandiger Lehm",       "points" = c(12, 10, 09, 13 ) ), 
                "S"   = list( "name" = "Sand",                "points" = c(17, 18, 20, 16 ) ), 
                "alS" = list( "name" = "anlehmiger Sand",     "points" = c(18, 19, 21, 15, 16, 20 ) ), 
                "lS"  = list( "name" = "lehmiger Sand",       "points" = c(19, 11, 12, 13, 14, 15, 21 ) ), 
                "T"   = list( "name" = "Ton",                 "points" = c(03, 01, 02 ) ), 
                "uT"  = list( "name" = "schluffiger Ton",     "points" = c(05, 26, 02, 04 ) ), 
                "lT"  = list( "name" = "lehmiger Ton",        "points" = c(03, 26, 05, 06 ) ), 
                "sT"  = list( "name" = "sandiger Ton",        "points" = c(07, 03, 06 ) ), 
                "U"   = list( "name" = "Schluff",             "points" = c(25, 22, 23, 24 ) ), 
                "UL"  = list( "name" = "Schlufflehm",         "points" = c(09, 05, 04, 08 ) ), 
                "lU"  = list( "name" = "lehmiger Schluff",    "points" = c(14, 13, 09, 08, 23, 22, 25 ) ) 
                #
            ),  #
            
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = c(F,T,NA), 
            tlr.an          = c(45,90,45), 
            
            blr.tx      = c("CLAY","SILT","SAND"), 
            
            base.css.ps.lim = c(0,2,63,2000), 
            tri.css.ps.lim  = c(0,2,63,2000), 
            
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            
            text.sum        = 100 
        ),  #
        #
        "DE.TGL85.TT" = list( # GDR Arable soils TRIANGLE PARAMETERS :
            #
            main = "TGL 24300-05, landwirtschaftliche Boeden (DE)", 
            #
            #                 The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes. Definitions follow the GDR Standard for arable soils
            #                   (TGL 24300-05:1985-06; Aufnahme landwirtschaftlich genutzter Standorte):
            #                  P01    P02    P03    P04    P05    P06    P07    P08    P09    P10    P11    P12    P13    
            #                  P14    P15    P16    P17    P18    P19    P20    P21    P22    P23    P24    P25    P26
            "tt.points" = data.frame( 
                "CLAY"  =   c( 1.000, 0.500, 0.500, 0.500, 0.300, 0.300, 0.300, 0.300, 0.180, 0.180, 0.180, 0.080, 0.080, 
                               0.000, 0.000, 0.140, 0.110, 0.080, 0.050, 0.000, 0.000, 0.000, 0.050, 0.085, 0.000, 0.050),
                "SILT"  =   c( 0.000, 0.500, 0.300, 0.000, 0.700, 0.500, 0.200, 0.000, 0.820, 0.500, 0.000, 0.920, 0.720, 
                               0.800, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.150, 0.300, 0.300, 0.300, 0.500, 0.500),  
                "SAND"  =   c( 0.000, 0.000, 0.200, 0.500, 0.000, 0.200, 0.500, 0.700, 0.000, 0.320, 0.820, 0.000, 0.200, 
                               0.200, 0.000, 0.860, 0.890, 0.920, 0.950, 1.000, 0.850, 0.700, 0.650, 0.615, 0.500, 0.450)  
            ),  #
            # 
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "rS"    = list( "name" = "reiner Sand",                "points" = c(19, 20, 21 ) ), 
                "l''S"  = list( "name" = "sehr schwach lehmiger Sand", "points" = c(18, 19, 21, 22 ) ), 
                "l'S"   = list( "name" = "schwach lehmiger Sand",      "points" = c(17, 18, 22, 23 ) ), 
                #"lS"   = list( "name" = "stark lehmiger Sand",        "points" = c(17, 18, 20, 16 ) ), 
                "lS"    = list( "name" = "stark lehmiger Sand",        "points" = c(16, 17, 23, 24 ) ), 
                "uS"    = list( "name" = "schluffiger Sand",           "points" = c(24, 23, 22, 25, 26 ) ), 
                "U"     = list( "name" = "Schluff",                    "points" = c(13, 14, 15, 12 ) ), 
                "lU"    = list( "name" = "lehmiger Schluff",           "points" = c(10, 26, 25, 14, 13, 12, 09 ) ), 
                "sL"    = list( "name" = "sandiger Lehm",              "points" = c(11, 16, 24, 26,  10 ) ), 
                "L"     = list( "name" = "Lehm",                       "points" = c(08, 11, 10, 06, 07 ) ), 
                "UL"    = list( "name" = "Schlufflehm",                "points" = c(06, 10, 09, 05 ) ), 
                "uT"    = list( "name" = "schluffiger Ton",            "points" = c(03, 06, 05, 02 ) ), 
                "lT"    = list( "name" = "lehmiger Ton",               "points" = c(04, 07, 06, 03 ) ), 
                "sT"    = list( "name" = "sandiger Ton",               "points" = c(04, 08, 07 ) ),
                "T"     = list( "name" = "Ton",                        "points" = c(01, 04, 03, 02 ) )
                #
            ),  #
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = c(F,T,NA), 
            tlr.an          = c(45,90,45), 
            #
            blr.tx          = c("CLAY","SILT","SAND"), 
            # 
            base.css.ps.lim = c(0,2,63,2000), 
            tri.css.ps.lim  = c(0,2,63,2000), 
            # 
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  
        
        "USDA1911" = list(
            # USDA 1911 (M. Whitney, 1911)
            # Courtesy of Nic Jelinski University of Minnesota, USA
            # 2014-04-23
            
            "main"  = "USDA 1911 (M. Whitney, 1911) - Ternary Plot", 
            
            "tt.points" = read.table( sep = "", header = TRUE, 
                text = "CLAY SILT SAND
                        0.0  0.0  1.0       # 1
                        0.2  0.0  0.8       # 2
                        0.0  0.2  0.8       # 3
                        0.5  0.0  0.5       # 4
                        0.3  0.2  0.5       # 5
                        0.2  0.3  0.5       # 6
                        0.0  0.5  0.5       # 7
                        0.2  0.5  0.3       # 8
                        0.3  0.5  0.2       # 9
                        0.3  0.7  0.0       # 10
                        0.2  0.8  0.0       # 11
                        0.0  1.0  0.0       # 12
                        1.0  0.0  0.0" ),   # 13
            
            "tt.polygons"   = list( 
                "C" = list( 
                    "name"      = "Clay", 
                    "points"    = c( 4, 13, 10,  5 ) 
                ),  
                "SaC" = list(   # INSTEAD OF SC (can be mistaken for Silty Clay)
                    "name"      = "Sandy Clay", 
                    "points"    = c( 2, 4, 5 ) 
                ),  
                "CL" = list( 
                    "name"      = "Clay Loam", 
                    "points"    = c( 5, 9, 8, 6 ) 
                ),  
                "SaCL" = list(  # INSTEAD OF SiCL
                    "name"      = "Sandy Clay Loam", 
                    "points"    = c( 8, 9, 10, 11 ) 
                ),  
                "L" = list( 
                    "name"      = "Loam", 
                    "points"    = c( 7, 6, 8 ) 
                ),  
                "SiL" = list( 
                    "name"      = "Silt Loam", 
                    "points"    = c( 7, 8, 11, 12 ) 
                ),  
                "SaL" = list(   # INSTEAD OF SL (can be mistaken for Silty Loam)
                    "name"      = "Sandy Loam", 
                    "points"    = c( 3, 2, 5, 7 ) 
                ),  
                "Sa" = list(    # INSTEAD OF S (Can be mistaken for Silt)
                    "name"      = "Sand", 
                    "points"    = c( 1, 2, 3 ) 
                )   
            ),  
            "blr.clock"         = c( FALSE, TRUE, NA ), 
            "tlr.an"            = c( 45, 90, 45 ), 
            "blr.tx"            = c( "SILT", "CLAY", "SAND" ), 
            "base.css.ps.lim"   = c( 0, 5, 50, 2000 ), 
            "tri.css.ps.lim"    = c( 0, 5, 50, 2000 ), 
            "unit.ps"           = quote( bold( mu ) * bold( "m" ) ), 
            "unit.tx"           = quote( bold( "%" ) ), 
            "text.sum"          = 100
        ),  
        
        BRASIL.TT = list(  #  Brazilian Triangle parameters (Lemos and Santos 1996)
            # Lemos, R. C. & Santos, R. D. Manual de descricao e 
            # coleta de solo no campo. 3a ed. Campinas, Sociedade 
            # Brasileira de Ciencia do solo, 1996.
            
            # Information is a courtesy of Rodolfo Marcondes Silva 
            # Souza, UFPE, Brasil (base-triangle is USDA triangle, 
            # modified for Brasil)
            
            main            = "Brasil - Lemos & Santos (1996)", 
            # 
            #                The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  = P01    P02    P03    P04    P05    P06    P07    P08    P09    P10    P11    P12
            #                  = P13    P14    P15    P16    P17    P18    P19    P20    P21    P22    P23   
            #                  = P24    P25    P26    P27 (submits)
            "tt.points"     = data.frame( 
                "CLAY"      = c( 0.550, 0.600, 0.350, 0.350, 0.400, 0.400, 0.400, 0.200, 0.200, 0.275, 0.275, 0.275,  
                                 0.275, 0.150, 0.100, 0.075, 0.075, 0.125, 0.125, 0.000, 0.000, 0.000, 0.000,         
                                 1.000, 0.000, 0.000, 0.600 ),  
                            #
                "SILT"      = c( 0.000, 0.400, 0.000, 0.200, 0.150, 0.400, 0.600, 0.000, 0.275, 0.275, 0.500, 0.525,  
                                 0.725, 0.000, 0.000, 0.400, 0.500, 0.800, 0.875, 0.150, 0.300, 0.500, 0.800,         
                                 0.000, 0.000, 1.000, 0.000 ),  
                            #
                "SAND"      = c( 0.450, 0.000, 0.650, 0.450, 0.450, 0.200, 0.000, 0.800, 0.525, 0.450, 0.225, 0.200,  
                                 0.000, 0.850, 0.900, 0.525, 0.425, 0.075, 0.000, 0.850, 0.700, 0.500, 0.200,         
                                 0.000, 1.000, 0.000, 0.400 )  
            ),  #
            #
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"   = list( 
                "MA"        = list( "name" = "muito argilosa",        "points" = c(02,24,27)                  ), 
                "A"         = list( "name" = "argila",                "points" = c(27,01,05,06,02)            ), 
                "As"        = list( "name" = "argila siltosa",        "points" = c(02,06,07)                  ), 
                "AAr"       = list( "name" = "argila arenosa",        "points" = c(01,03,04,05)               ), 
                "FA"        = list( "name" = "franco argiloso",       "points" = c(05,04,10,11,12,06)         ), 
                "FAS"       = list( "name" = "franco argilo siltoso", "points" = c(06,12,13,07)               ), 
                "FAAr"      = list( "name" = "franco argilo arenoso", "points" = c(03,08,09,10,04)            ), 
                "F"         = list( "name" = "franco",                "points" = c(10,09,16,17,11)            ), 
                "FS"        = list( "name" = "franco siltoso",        "points" = c(11,17,22,23,18,19,13,12)   ), 
                "FAr"       = list( "name" = "franco arenoso",        "points" = c(08,14,21,22,17,16,09)      ), 
                "S"         = list( "name" = "silte",                 "points" = c(18,23,26,19)               ), 
                "ArF"       = list( "name" = "areia franca",          "points" = c(14,15,20,21)               ), 
                "Ar"        = list( "name" = "areia",                 "points" = c(15,25,20)                  )  
            ),  
            #
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = rep(T,3), 
            tlr.an          = c(60,60,60), 
            #
            blr.tx      = c("SAND","CLAY","SILT"), 
            # 
            base.css.ps.lim = c(0,2,50,2000), 
            tri.css.ps.lim  = c(0,2,50,2000), 
            #
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            #
            text.sum        = 100 
        ),  #
        
        "SiBCS13.TT" = list(  
            # Subagrupamento Textural SiBCS 2013 parameters (Embrapa 2013)
            #   Embrapa. Sistema Brasileiro de Classificacao de Solos /
            #   Humberto Golcalves dos Santos ... [et al.]. 3a ed. rev. ampl.
            #   Brasilia, DF: Embrapa, 2013.
            
            # Information is a courtesy of Jose Lucas Safanelli and
            # Alexandre ten Caten, UFSC Curitibanos, Brasil.
            
            # main            = "Subagrupamento textural SiBCS 2013 - Embrapa 2013", 
            main            = "SiBCS 2013 (Embrapa)", # Shorter title and more international?
            # 
            #                The list below specify the CSS coordinates of the different POINTS
            #                   that are used to draw soil texture classes. One points can be 
            #                   used by several classes :
            #                  = P01    P02    P03    P04    P05    P06    P07    P08    P09    P10    P11    P12
            #                  = P13    P14    P15    P16    P17    P18    (submits)
            "tt.points"     = data.frame( 
                "CLAY"      = c( 1.000, 0.600, 0.600, 0.350, 0.350, 0.350, 0.350, 0.200, 0.200, 0.250, 0.150, 0.100,  
                                 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 ),
                            #
                "SILT"      = c( 0.000, 0.000, 0.400, 0.000, 0.175, 0.500, 0.650, 0.000, 0.250, 0.275, 0.000, 0.000,  
                                 0.000, 0.150, 0.300, 0.475, 0.850, 1.000 ),
                            #
                "SAND"      = c( 0.000, 0.400, 0.000, 0.650, 0.475, 0.150, 0.000, 0.800, 0.550, 0.475, 0.850, 0.900,  
                                 1.000, 0.850, 0.700, 0.525, 0.150, 0.000 )
            ),  #
            
            #   Abreviations;       Names of the texture cl;    Points marking the class limits (points specified above)
            "tt.polygons"  = list( 
                "MA"       = list( "name" = "muito argilosa", "points" = c(01,02,03)          ), 
                "A"        = list( "name" = "argilosa",       "points" = c(02,03,07,06,05,04) ), 
                "S"        = list( "name" = "siltosa",        "points" = c(06,07,18,17)       ), 
                "MeS"      = list( "name" = "media siltosa",  "points" = c(05,06,17,16,09,10) ), 
                "MeA"      = list( "name" = "media argilosa", "points" = c(04,05,10,09,08)    ), 
                "MeAr"     = list( "name" = "media arenosa",  "points" = c(08,09,16,15,11)    ), 
                "ArMe"     = list( "name" = "arenosa media",  "points" = c(11,15,14,12)       ), 
                "MAr"      = list( "name" = "muito arenosa",  "points" = c(12,14,13)          )   
            ),  
            
            # Triangle specific parameters for triangle geometry / appearance
            #   See general parameters above for detailed description of them
            blr.clock       = rep(T,3), 
            tlr.an          = c(60,60,60), 
            
            blr.tx          = c("SAND","CLAY","SILT"), 
             
            base.css.ps.lim = c(0,2,50,2000), 
            tri.css.ps.lim  = c(0,2,50,2000), 
            
            unit.ps         = quote(bold(mu) * bold('m')), 
            unit.tx         = quote(bold('%')), 
            
            text.sum        = 100, 
            
            #   New parameter
            class.lab.cex   = 2/3
        )   
        
        # +-------------------------------------------------------------------------+
        # | END(SCRIPT PARAMETERS SPECIFICATION)                                    |
        # +-------------------------------------------------------------------------+
        
    )   #
)   #



# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# | LIST:   TT.par.bkp                  |
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# [ TT.par.bkp :: A backup of "TT.par" parameter list, used for resetting "TT.par"
assign( 
    envir   = TT.env,
    x       = "TT.par.bkp", 
    value   = get( x = "TT.par",  envir = TT.env )
)   #
#   get( "TT.par.bkp", envir = TT.env )






TT.set <- function(# Function to change / set the default package parameters. 
### Function to change / set the default package parameters as they 
### are stored in the list TT.par in the environment TT.env. Use 
### this function to change some deafult parameters for all the 
### current R cession. Many functions of soiltexture take some of 
### their parameter values in TT.par.

...,
### List of parameters and value in the form "par.name1" = par.value1, 
### "par.name2" = par.value2... List of parameters to change.

reset=FALSE,
### Single logical. If set to TRUE the parameter list is reset to default

par.list="TT.par",
### Single character. Name of the list containing the parameters

bkp.par.list="TT.par.bkp",
### Single character. Name of the backuped list containing the default parameters

par.env=TT.env
### An R environment. Name of the environment containing the parameter lists (no quotes)

){  #
    argz <- list(...)
    #
    # Basic checkup:
    if( !reset & (length(argz) == 0) )
    {   stop("No argument specified!")  }
    #
    if( reset & (length(argz) != 0) )
    {   stop("'reset' is not compatible with other parameters specification")   }
    #
    argz.nm <- names(argz)
    #
    if( (length(argz) != 0) & any(is.null(argz.nm)) )
    {   stop("parameterS should not be text names of variables: use TT.get() to get variables") }
    #
    # Reset option:
    if( reset )
    {   #
        assign( 
            envir   = par.env,
            x       = par.list, 
            value   = get( x = bkp.par.list,  envir = par.env )
        )   #
    }else{  # Non reset options:
        # Import the environment variables (function internal)
        argz.list   <- get( x = par.list,  envir = par.env )
        #
        # Basic check:
        if(   any(  !( names(argz) %in% names(argz.list) )  )   )
        {   stop("At least one specified parameter doesn't exist!") }
        #
        # SETTING the parameters value:
        for( i in 1:length(argz) )
        {   #
            if( is.null( argz[[i]] ) )
            {   #
                argz.list[  names(argz)[i]  ] <- list( argz[[i]] ) 
            }else{ 
                argz.list[[ names(argz)[i] ]] <- argz[[i]] 
            }   #
        }   #
        # Re-export the environment variables
        assign( 
            envir   = par.env,
            x       = par.list, 
            value   = argz.list 
        )   #
        #
        return( invisible( argz.list ) )
    }   #
}   #






TT.get <- function(# Function to retrieve / get the default package parameters. 
### Function to retrieve / get the default package parameters. 

    ...,                            # List of parameters to change
    par.list        = "TT.par",     # Name of the list containing the parameters
    bkp.par.list    = "TT.par.bkp", # Name of the backuped list containing the default parameters
    par.env         = TT.env        # name of the environment containing the parameter lists
){  #
    argz <- list(...)
    #
    # Import the environment variables (function internal)
    argz.list   <- get( x = par.list,  envir = par.env )
    #
    # Basic checkup:
    if( length(argz) == 0 )
    {   #
        res <- argz.list
    }else{ 
        #
        argz.nm <- names(argz)
        #
        if( any(!is.null(argz.nm)) ) 
        {  stop("parameterS should be text names of variables: use TT.set() to set variables")  }
        #
        # Flattern argument list (to text vector)
        argz <- unlist(argz) 
        #
        # Basic check:
        if(  any( !(argz %in% names(argz.list)) )  ) 
        {   stop("At least one specified parameter doesn't exist!") }
        #
        # RETRIEVING the parameters value (in the order they were asked)
        res <- lapply(X=argz,FUN=function(X,dat){dat[[X]]},dat=argz.list)
        # res <- argz.list[ argz ]
        if( length(res) == 1 )
        {   res <- res[[1]] } 
    }   #
    return( res )
}   #






TT.add <- function(# Function to add a new default package parameters. 
### Function to add a new default package parameters. Mostly used 
### to add a new texture triangle definition.

    ...,                            # List of parameters to change
    par.list        = "TT.par",     # Name of the list containing the parameters
    bkp.par.list    = "TT.par.bkp", # Name of the backuped list containing the default parameters
    par.env         = TT.env        # name of the environment containing the parameter lists
){  #
    argz    <- list(...) 
    argz.nm <- names(argz) 
    #
    # Import the environment variables (function internal)
    argz.list       <- get( x = par.list,  envir = par.env ) 
    argz.list.nm    <- names( argz.list ) 
    #
    # Basic checkup:
    if( length(argz) == 0 )
    {   #
        stop( 
            "TT.add must have at least one '...' option specified (named list). Now zero length" 
        )   #
    }   #
    #
    if( is.null( argz.nm ) | any( argz.nm == "" ) ) 
    {   #
        stop( 
            paste( 
                sep = "", 
                "TT.add ... option must be a named list (eventually length 1),\n", 
                "with _names_ different from ''.\n", 
                "Now names(...) = ", paste(argz.nm,collapse=" ")  
            )   #
        )   #
    }   #
    #
    if( any( duplicated( argz.nm ) ) ) 
    {   #
        stop( 
            paste( 
                sep = "", 
                "TT.add ... option must be a named list (eventually length 1),\n", 
                "with unique (non duplicated) names.\n", 
                "Now names(...) = ", paste(names(argz),collapse=" ")  
            )   #
        )   #
    }   #
    #
    if( any( argz.nm %in% argz.list.nm ) ) 
    {   #
        stop( 
            paste( 
                sep = "", 
                "TT.add ... option must be a named list (eventually length 1),\n", 
                "with names that are not already in the default argument list (TT.par).\n", 
                "Now names(...) = ", paste(names(argz),collapse=" "), "\n", 
                "If you want to replace/change an existing argument, use TT.set() instead"  
            )   #
        )   #
    }   #
    # 
    argz.list   <- c(argz.list,argz) 
    sel.vec     <- (length(argz.list)-length(argz)+1):length(argz.list) 
    names( argz.list )[ sel.vec ]   <- argz.nm
    # 
    assign( 
        envir   = par.env,
        x       = par.list, 
        value   = argz.list 
    )   #
}   #






TT.str <- function(# Internal. Stretch or reshape the range of value of some data set. 
### Function to 'stretch' or reshape the range of value of some data set. Usefull for cex parameter in plot. 
##keywords<< internal

    x, 
    str.min = 0, 
    str.max = 1
){  #
    aa <- (str.min - str.max)/( min(x) - max(x) )
    bb <- str.min - aa*min(x)
    aa*x + bb
}   #
#   TT.str(50:100,1,0)






# # TT.nightC <- function(# Internal. Inverse RGB values of a vector of colors.
# # ### Inverse RGB values of a vector of colors.

# #     cl,         # vector of colors, html stype "#808080"
# #     ic  = TRUE  # Really inverse colors ?
# # ){  #
# #     if( ic ) 
# #     {   #
# #         cl <- col2rgb(
# #             col     = cl, 
# #             alpha   = FALSE 
# #         )   #
# #         #
# #         cl <- apply( 
# #             X       = cl, 
# #             MARGIN  = 2, 
# #             FUN     = function(X){ 
# #                 rep(255,3) - X
# #             }   #
# #         )   #
# #         #
# #         cl <- cl/255
# #         #
# #         rgb( 
# #             red     = cl["red",], 
# #             green   = cl["green",], 
# #             blue    = cl["blue",], 
# #         )   #
# #     }else{ 
# #         cl
# #     }   #
# # }   #






TT.gen.op.set  <- function(# Internal. Retrieve and set default values from options. 
### Retrieve and set default values from options (that do _not_ superseed par()). 
##keywords<< internal

    param, 
    assign.op   = TRUE, 
    p.env       = parent.frame() 
){  #
    # Get the parameter values:
    param.val   <- lapply( 
        X   = param, 
        FUN = function(X){ 
            get(X,envir=p.env) 
        }   #
    )   #
    names(param.val)    <- param 
    #
    # Find null parameters values:
    null.param  <- unlist(  lapply( 
            X   = param.val, 
            FUN = function(X){ 
                any( is.null( X ) ) 
            }   #
    )   )   #
    #
    if( any(null.param) )
    {   #
        # Get the default value in the options list: 
        param.val[ null.param ] <- lapply( 
            X   = param[ null.param ], 
            FUN = function(X){ 
                TT.get(X) 
            }   #
        )   #
        #
        # Assign the values in the higher level function
        if( assign.op )
        {   #
            silent  <- lapply( 
                X   = 1:length(param[ null.param ]), 
                FUN = function(X){ 
                    assign( 
                        x       = param[ null.param ][X], 
                        value   = param.val[ null.param ][[X]], 
                        envir   = p.env 
                    )   #
                }   #
            )   #
        }   #
    }   #
    return( invisible( param.val ) )
}   #

#   test.fun1   <- function( 
#       cex         = NULL, 
#       cex.lab     = 2, 
#       col.axis    = "blue", 
#       font.axis   = NULL  
#   ){  #
#       invres <- TT.gen.op.set(c("cex","cex.lab","col.axis","font.axis"))
#       #
#       list("cex"=cex,"cex.lab"=cex.lab,"col.axis"=col.axis,"font.axis"=font.axis,invres) 
#   }   #
#   test.fun1()






TT.par.op.set  <- function(# Internal. Retrieve and set default values from options with default in "par()". 
### Retrieve and set default values from options with default in "par()"
##keywords<< internal

    param, 
    assign.op   = TRUE, 
    p.env       = parent.frame() 
){  #
    param.val   <- TT.gen.op.set( 
        param       = param, 
        assign.op   = assign.op, 
        p.env       = p.env  
    )   #
    #
    # Find null parameters values:
    null.param  <- unlist(  lapply( 
            X   = param.val, 
            FUN = function(X){ 
                any( is.null( X ) ) 
            }   #
    )   )   #
    #
    if( any(null.param) )
    {   #
        # Get the default value in par() 
        param[ param == "family.op" ]       <- "family" # for compatibility with family() 
        param.val[ null.param ]             <- par(param[ null.param ]) 
        param[ param == "family" ]          <- "family.op" 
        # 
        # Assign the values in the higher level function
        if( assign.op )
        {   #
            silent  <- lapply( 
                X   = 1:length(param[ null.param ]), 
                FUN = function(X){ 
                    assign( 
                        x       = param[ null.param ][X], 
                        value   = param.val[ null.param ][[X]], 
                        envir   = p.env 
                    )   #
                }   #
            )   #
        }   #
    }   #
    return( invisible( param.val ) )
}   #

#   test.fun    <- function( 
#       cex         = NULL, 
#       cex.lab     = 2, 
#       col.axis    = "blue", 
#       font.axis   = NULL, 
#       family.op   = NULL  
#   ){  #
#       invres <- TT.par.op.set(c("cex","cex.lab","col.axis","font.axis","family.op"))
#       #
#       list("cex"=cex,"cex.lab"=cex.lab,"col.axis"=col.axis,"font.axis"=font.axis,"family.op"=family.op,invres) 
#   }   #
#   test.fun()






TT.auto.set    <- function(# Internal. Retrieve and set default values for parameters (par() or not), when NULL.
### Retrieve and set default values for parameters (par() or not), when NULL.
##keywords<< internal

    fun         = sys.function(which=-1), 
    assign.op   = TRUE, 
    p.env       = parent.frame(), 
    set.par     = TRUE 
){  #
    param   <- names( formals(fun) ) 
    #
    if( set.par )
    {   #
        sel.par <- (param %in% names(par())) 
        #
        l1  <- TT.par.op.set( 
            param       = param[ sel.par ], 
            assign.op   = assign.op, 
            p.env       = p.env
        )   #
    }else{ 
        sel.par <- rep(FALSE,length(param)) 
        #
        l1  <- vector()
    }   #
    #
    sel.TT  <- ((param %in% names(TT.get())) & !sel.par)
    #
    l2  <- TT.gen.op.set( 
        param       = param[ sel.TT ], 
        assign.op   = assign.op, 
        p.env       = p.env  
    )   #
    #
    return( invisible( c(l1,l2) ) ) 
}   #

#   test.fun    <- function( 
#       cex         = NULL, 
#       cex.lab     = 2, 
#       col.axis    = "blue", 
#       font.axis   = NULL, 
#       blr.clock   = NULL, 
#       tlr.an      = c(50,60,70) 
#   ){  #
#       invres <- TT.auto.set(set.par=TRUE) 
#       #
#       list( 
#           invres, 
#           lapply( 
#               X   = names(formals(test.fun)), 
#               FUN = function(X){ 
#                   get(X)
#               }   #
#           )   #
#       )   #
#   }   #
#   test.fun()






# # TT.inv.par <- function(# Internal. Same as the par() function, but reverse colors.
# # ### Same as the par() function, but reverse colors.

# #     ic      = TRUE, # Really inverse colors ? To be used for autonated shift
# #     par.opt = c("bg","col","col.axis","col.lab","col.main","col.sub","fg"), 
# #     ...
# # ){  #
# #     old.par <- par(no.readonly=TRUE) 
# #     #
# #     dots    <- list(...) 
# #     #
# #     par.opt <- par.opt[ !(par.opt %in% names(dots)) ] 
# #     #
# #     if( ic )
# #     {   #
# #         cl          <- TT.nightC( cl = old.par[ par.opt ], ic = ic ) 
# #         cl          <- as.list( cl ) 
# #         names( cl ) <- par.opt 
# #         cl          <- c(dots,cl) 
# #     }else{ cl <- dots }
# #     #
# #     do.call( what = "par", args = cl )
# #     #
# #     return( invisible( old.par ) )
# # }   #






TT.DJ.col <- function(# Internal. A function to obtaine a weight average 'mix' of different colors!
### A function to obtaine a weight average 'mix' of different colors!
##keywords<< internal

    cl,             # vector of colors, html stype "#808080"
    w,              # vector of weight corresponding to the colors
    gray.l  = FALSE # if TRUE Produce a gray level color, instead of a 'colored' color
){  #
    cl  <- col2rgb( cl, alpha = FALSE )
    #
    m.cl    <- apply( 
        X       = cl, 
        MARGIN  = 1, 
        FUN     = function(X){ 
            weighted.mean(x=X,w=w) 
        }   #
    )   #
    #
    if( gray.l ){ m.cl[] <- rep(mean(m.cl),3) }     # 1:3 stands here in case of alpha value...
    #
    rgb( 
        red           = m.cl["red"], 
        green         = m.cl["green"], 
        blue          = m.cl["blue"], 
        maxColorValue = 255  
    )   #
}   #





TT.col2hsv  <- function(# Internal. Convert any colors to hsv. 
### Convert any colors to hsv. Wrapper around rgb2hsv() and col2rgb(). 
##keywords<< internal

    col 
){  #
    t(  #
        rgb2hsv( 
            col2rgb( 
                col     = col, 
                alpha   = FALSE
            ),  #
            #gamma          = 1, 
            maxColorValue   = 255  
        )   #
    )   #
}   #






TT.blr.tx.check <- function(# Internal. Check the consistency between blr.tx and css.names. 
### Check the consistency between blr.tx and css.names. All values 
### in blr.tx should be found in css.names and vice-versa.
##keywords<< internal

    blr.tx, 
    css.names  
){  #
    css <- TT.get("css.names")
    #
    names(css.names)    <- css
    #
    if( !all( blr.tx %in% css ) | !all( css %in% blr.tx ) )
    {   #
        stop( 
            paste( 
                sep = "", 
                "Every blr.tx (", paste(blr.tx,collapse=", "), ")\n", 
                "\t should be one of ", paste(css,collapse=", "), "\n", 
                "\t and vice-versa."  
            )   #
        )   #
    }   #
    #
    blr.tx.nm       <- blr.tx 
    blr.tx          <- css.names[ blr.tx ] 
    names(blr.tx)   <- blr.tx.nm 
    #
    return( blr.tx ) 
}   #

#     TT.blr.tx.check(c("SAND","CLAY","SILT"),c("ARGILE","LIMON","SABLE"))
#     TT.blr.tx.check(c("SAND","CLAY","SILT"),c("CLAY","SILT","SAND"))






TT.blr.ps.lim <- function(# Internal. Create a tabular version of clay silt sand particle size limits. 
### Create a tabular version of clay silt sand particle size limits. 
##keywords<< internal

    blr.tx, 
    css.ps.lim  
){  #
    css.ps.lim  <- do.call( 
        what    = "cbind", 
        args    = lapply( 
            X   = 1:3, 
            FUN = function(X){ 
                c(css.ps.lim[X],css.ps.lim[X+1])  
            }   #
        )   #
    )   #
    #
    colnames(css.ps.lim)    <- TT.get("css.names") 
    rownames(css.ps.lim)    <- c("ps.min","ps.max") 
    #
    css.ps.lim              <- css.ps.lim[, blr.tx ] 
    #
    colnames(css.ps.lim)    <- c("B","L","R") 
    #
    return( css.ps.lim ) 
}   #

#     TT.blr.ps.lim( blr.tx = c("CLAY","SILT","SAND"), css.ps.lim = c(0,2,50,2000) ) 
#     TT.blr.ps.lim( blr.tx = c("SAND","CLAY","SILT"), css.ps.lim = c(0,2,50,2000) ) 






TT.geo.set  <- function(# Internal. Takes "geo" values and assign them individually in the parent function. 
### Takes "geo" values and assign them individually in the parent function. 
##keywords<< internal

    geo, 
    p.env   = parent.frame()  
){  #
    silent  <- lapply( 
        X   = names(geo), 
        FUN = function(X){ 
            assign( 
                x       = X, 
                value   = geo[[X]], 
                envir   = p.env 
            )   #
        }   #
    )   #
}   #

#   rm("blr.clock","tlr.an","blr.tx","text.sum","blr.psize.lim")
#   test.fun    <- function( 
#       geo = list( 
#           "blr.clock"     = rep(T,3), 
#           "tlr.an"        = rep(60,3), 
#           "blr.tx"        = c("CLAY","SILT","SAND"), 
#           "text.sum"      = 100, 
#           "blr.psize.lim" = 50  
#       )   #
#   ){  #
#       TT.geo.set( 
#           geo     = geo  
#           #p.env  = parent.frame()  
#       )   #
#       #
#       return( list( 
#           "blr.clock"     = blr.clock, 
#           "tlr.an"        = tlr.an, 
#           "blr.tx"        = blr.tx, 
#           "text.sum"      = text.sum, 
#           "blr.psize.lim" = blr.psize.lim  
#       )   )   #
#   }   #
#   test.fun() ; blr.clock 

#   rm("blr.clock","tlr.an","blr.tx","text.sum","blr.psize.lim")
#   test.fun    <- function( 
#       geo = list( 
#           "blr.clock"     = rep(T,3), 
#           "tlr.an"        = rep(60,3), 
#           "blr.tx"        = c("CLAY","SILT","SAND"), 
#           "text.sum"      = 100, 
#           "blr.psize.lim" = 50  
#       )   #
#   ){  #
#       TT.geo.set( 
#           geo     = geo  
#           #p.env  = environment()  
#       )   #
#       #
#       return( list( 
#           "blr.clock"     = blr.clock, 
#           "tlr.an"        = tlr.an, 
#           "blr.tx"        = blr.tx, 
#           "text.sum"      = text.sum, 
#           "blr.psize.lim" = blr.psize.lim  
#       )   )   #
#   }   #
#   test.fun() ; blr.clock 
    
#   rm("blr.clock","tlr.an","blr.tx","text.sum","blr.psize.lim")
#   test.fun    <- function( 
#       geo = list( 
#           "blr.clock"     = rep(T,3), 
#           "tlr.an"        = rep(60,3), 
#           "blr.tx"        = c("CLAY","SILT","SAND"), 
#           "text.sum"      = 100, 
#           "blr.psize.lim" = 50  
#       )   #
#   ){  #
#       TT.geo.set( 
#           geo     = geo 
#       )   #
#       #
#       return( list( 
#           "blr.clock"     = blr.clock, 
#           "tlr.an"        = tlr.an, 
#           "blr.tx"        = blr.tx, 
#           "text.sum"      = text.sum, 
#           "blr.psize.lim" = blr.psize.lim  
#       )   )   #
#   }   #
#   test.fun() ; blr.clock 

#   test.fun    <- function( p.env = environment() ){ p.env }
#   test.fun2   <- function(){ test.fun() } 
#   test.fun2() 






TT.geo.get  <- function(# Internal. Retrieve and return the geometrical parameters from a list of parameter values (NULL or not).
### Retrieve and return the geometrical parameters from a list of parameter values (NULL or not).
##keywords<< internal

    class.sys       = NULL,  
    blr.clock       = NULL,  
    tlr.an          = NULL,  
    blr.tx          = NULL,  
    text.sum        = NULL,  
    base.css.ps.lim = NULL   
){  #
    if( is.null(class.sys) ){ class.sys <- TT.get("class.sys") } 
    
    if( class.sys == "FAO50.TT" ){ 
        warning( "class.sys = 'FAO50.TT' must be replaced by class.sys = 'HYPRES.TT'. See the package vignette." )
    }   
    
    geo.par         <- c("blr.clock","tlr.an","blr.tx","text.sum","base.css.ps.lim") 
    # 
    p.env           <- environment() 
    #
    null.geo.par    <- unlist(  lapply( 
            X   = geo.par, 
            FUN = function(X){ 
                is.null( get(x=X,envir=p.env) ) 
            }   #
    )   )   #
    
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Attributes either classes-system or default values 
    # to triangle geometry parameters (and others)
    
    if( any( null.geo.par ) )
    {   #
        geo.par <- geo.par[ null.geo.par ] 
        #
        if( class.sys != "none" )
        {   #
            # Retrieve classes-system (texture triangle) parameters:
            TT.data <- TT.get(class.sys) 
            #
            silent  <- lapply( 
                X   = geo.par, 
                FUN = function(X){ 
                    assign(x=X,value=TT.data[[X]],envir=p.env)
                }   #
            )   #
        }else{ 
            silent  <- lapply( 
                X   = geo.par, 
                FUN = function(X){ 
                    assign(x=X,value=TT.get(X),envir=p.env) 
                }   #
            )   #
        }   #
    }   #
    #
    return( list( 
            "blr.clock"         = blr.clock, 
            "tlr.an"            = tlr.an, 
            "blr.tx"            = blr.tx, 
            "text.sum"          = text.sum, 
            "base.css.ps.lim"   = base.css.ps.lim  
            #"blr.psize.lim"    = blr.psize.lim  
    )   )   #
}   #






TT.data.test <- function(# Test the validity of some soil texture data table (3 particle size classes). 
### Test the validity of some soil texture data table. (1) Test that 
### it is a data.frame or matrix, (2) Test that column names contains 
### 'css.names', (3) Test that there are no missing values, (4) that 
### all values are >= 0, (5) That the sum of the 3 particle size classes 
### is >= 'text.sum'*(1-'text.tol') or <= 'text.sum'*(1+'text.tol'). 
### 'tri.data' may contain other variables than the 3 textuer classes 
### (ignored).

    tri.data, 
    css.names   = NULL, 
    text.sum    = NULL, 
    text.tol    = NULL, 
    #
    tri.sum.tst = NULL, 
    tri.pos.tst = NULL  
){  #
    # Set rest of variables:
    TT.auto.set(set.par=FALSE) 
    #
    # 1. Check if tri.data is a matrix or a data.frame:
    if( !is.data.frame(tri.data) & !is.matrix(tri.data) )
    {   #
        stop("tri.data MUST be a data.frame, a matrix, or NULL") 
    }   #
    # 
    # 2. Check if columns names correspond to the list provided by css.names 
    if( !all( css.names %in% colnames(tri.data) ) )
    {   #
        stop(   paste( 
                sep = "", 
                "tri.data column names (", 
                paste(colnames(tri.data),collapse=", "), 
                ") don't correspond to ", 
                paste(css.names,collapse=", ") 
        )   )   #
    }   #
    # 
    # 3. Sub-select only the interest variables
    tri.data <- tri.data[,css.names] 
    #
    # 3b. Test the presence of missing values (error)
    if( tri.pos.tst ) 
    {   #
        row.na <- apply( 
            X       = tri.data,  
            MARGIN  = 1,  
            FUN     = function(X){ any( is.na( X ) ) } 
        )   #
        #
        if(  any( row.na )  )
        {   #
            print( tri.data[row.na,] ) 
            cat("\n") 
            stop( "No missing values are allowed in tri.data: check the data" ) 
        }   #   
    }   #
    #
    # 4. Test the presence of negative values (error)
    if( tri.pos.tst ) 
    {   #
        row.neg <- apply( 
            X       = tri.data,  
            MARGIN  = 1,  
            FUN     = function(X){ any(X < 0) } 
        )   #
        #
        if(  any( row.neg )  )
        {   #
            print( tri.data[row.neg,] ) 
            cat("\n") 
            stop( "Each of the 3 plotted variables should be >= 0: check the data" ) 
        }   #   
    }   #
    #
    if( tri.sum.tst ) 
    {   #
        # 5. Check if all row sum of variable triplets correspond to total value provided by text.sum
        # -  Compute the sum
        row.sums <- apply( 
            X       = tri.data,  
            MARGIN  = 1,  
            FUN     = sum
        )   #
        # -  Check if sum results are ok (+/- tolerance)
        sums.bad <- (row.sums < text.sum*(1-text.tol)) | (row.sums > text.sum*(1+text.tol))
        # 
        if(  any( sums.bad )  )
        {   #
            print( data.frame( tri.data[sums.bad,], "sum" = row.sums[sums.bad] ) ) 
            cat("\n") 
            stop( 
                paste( 
                    sep="", 
                    "The sum of the 3 plotted variables should be around ", 
                    text.sum, 
                    ": check the data, or change 'text.tol' parameter." 
                )   #
            )   #
        }   #
    }   #
}   #






TT.data.test.X <- function(# Test the validity of some soil texture data table (X particle size classes). 
### Test the validity of some soil texture data table. (1) Test that 
### it is a data.frame or matrix, (3) Test that there are no missing 
### values, (4) that all values are >= 0, (5) That the sum of the 
### X particle size class is >= 'text.sum'*(1-'text.tol') or <= 
### 'text.sum'*(1+'text.tol'). Contrary to TT.data.test() no test 
### are performed for the particle size classes and columns names, so 
### 'tri.data' should only contains texture data, and nothing else.

    tri.data,   # Only texture data here. No additionnal variables
    text.sum    = NULL, 
    text.tol    = NULL, 
    #
    tri.sum.tst = NULL, 
    tri.pos.tst = NULL  
){  #
    # Set rest of variables:
    TT.auto.set( set.par = FALSE ) 
    #
    # 1. Check if tri.data is a matrix or a data.frame:
    if( !is.data.frame(tri.data) & !is.matrix(tri.data) )
    {   #
        stop("tri.data MUST be a data.frame, a matrix, or NULL") 
    }   #
    #
    # 2. Test the presence of negative values (error)
    if( tri.pos.tst ) 
    {   #
        row.neg <- apply( 
            X       = tri.data,  
            MARGIN  = 1,  
            FUN     = function(X){ any(X < 0) } 
        )   #
        #
        if(  any( row.neg )  )
        {   #
            print( tri.data[row.neg,] ) 
            cat("\n") 
            stop( "Each of the n particle size classes must be >= 0: check the data" ) 
        }   #   
    }   #
    #
    if( tri.sum.tst ) 
    {   #
        # 3. Check if all row sum of variable triplets correspond to total value provided by text.sum
        # -  Compute the sum
        row.sums <- apply( 
            X       = tri.data,  
            MARGIN  = 1,  
            FUN     = sum
        )   #
        # -  Check if sum results are ok (+/- tolerance)
        sums.bad <- (row.sums < text.sum*(1-text.tol)) | (row.sums > text.sum*(1+text.tol))
        # 
        if(  any( sums.bad )  )
        {   #
            print( data.frame( tri.data[sums.bad,], "sum" = row.sums[sums.bad] ) ) 
            cat("\n") 
            stop( 
                paste( 
                    sep="", 
                    "The sum of the n particle size classes should be around ", 
                    text.sum, 
                    ": check the data, or change 'text.tol' parameter." 
                )   #
            )   #
        }   #
    }   #
}   #






TT.dia2phi <- function(# Internal. Convert a soil particle diameter dia [micro-meters] into phi = -log2(dia/1000)
### Convert a soil particle diameter dia [micro-meters] into 
### phi = -log2(dia). See also TT.phi2dia().
##keywords<< internal

 dia
### Particle size diameter in micro-meters (will be converted in milli-meters)

){  #
    return( -logb(dia/1000,base=2) ) 
}   #






TT.phi2dia <- function(# Internal. Convert a soil particle phi value into diameter dia [micro-meters]. 
### Convert a soil particle phi value into diameter dia [micro-meters]. 
### See also TT.dia2phi(). dia = (2^-phi)*1000. Not used by the package. 
##keywords<< internal

phi

){  #
    return( (2^-phi)*1000 )
}   #






TT.check.ps.lim <- function(# Internal. Check the consistency between 'base.ps.lim' and 'dat.ps.lim'. 
### Check the consistency between 'base.ps.lim' and 'dat.ps.lim'. 
### 5 tests performed.
##keywords<< internal

    base.ps.lim,  
    dat.ps.lim,  
    ps.lim.length=c(4,4)
### vector of 2 integers. Number of particle size classes + 1. c(base,dat)

){  #
    # if( length( base.ps.lim ) != length( dat.ps.lim ) ) 
    # {   #
    #     stop( paste( 
    #         sep="", 
    #         "The length of the 'base' particle size classes limits must be equal to\n", 
    #         "the length of the 'dat' particle size classes limits.\n", 
    #         "Either check the 'base' particle size classes limits vector,\n", 
    #         "or check number of column in tri.data.\n"  
    #     ) )   #
    # }   #
    #
    if( length( base.ps.lim ) != ps.lim.length[1] ) 
    {   #
        stop( paste( 
            sep="", 
            "The length of the 'base' particle size classes limits must be equal to\n", 
            ps.lim.length[1], " (number of particle size classes+1; from ps min to ps.max)\n", 
            "Actual value: ", length( base.ps.lim ), ".\n", 
            "Either check the 'base' particle size classes limits,\n", 
            "or check number of column in tri.data.\n"  
        ) )   #
    }   #
    #
    if( length( dat.ps.lim ) != ps.lim.length[2] ) 
    {   #
        stop( paste( 
            sep="", 
            "The length of the 'dat' particle size classes limits must be equal to\n", 
            ps.lim.length[2], " (number of particle size classes +1; from ps min to ps.max)\n", 
            "Actual value: ", length( dat.ps.lim ), ".\n", 
            "Either check the 'dat' particle size classes limits,\n", 
            "or check number of column in tri.data.\n"  
        ) )   #
    }   #
    #
    if( base.ps.lim[1] != dat.ps.lim[1] ) 
    {   #
        stop( paste( 
            sep="", 
            "The first value of the 'dat' particle size classes limits must be equal to\n", 
            "the first value of the 'base' particle size classes limits.\n", 
            "Actual value, base: ", base.ps.lim[1], ", dat: ", dat.ps.lim[1]  
        ) )   #
    }   #
    #
    if( base.ps.lim[ps.lim.length[1]] != dat.ps.lim[ps.lim.length[2]] ) 
    {   #
        stop( paste( 
            sep="", 
            "The last value of the 'dat' particle size classes limits must be equal to\n", 
            "the last value of the 'base' particle size classes limits.\n", 
            "Actual value, base: ", base.ps.lim[ps.lim.length[1]], ", dat: ", dat.ps.lim[ps.lim.length[2]]  
        ) )   #
    }   #
    #
    if( base.ps.lim[1] == 0 ) 
    {   #
        if( base.ps.lim[2] < dat.ps.lim[2] )
        stop( paste( 
            sep="", 
            "When the 1st value of 'dat' and 'base' particle size classes limits is 0\n", 
            "The 2nd value of the 'base' particle size classes limits must higher or equal to\n", 
            "the 2nd value of the 'dat' particle size classes limits.\n"  
        ) )   #
    }   #
}   #






TT.text.transf <- function(# Log-linear transformation of a soil texture data table between 2 particle size systems (3 classes).
### Log-linear transformation of a soil texture data table 
### ('tri.data') from one 
### particle size system ('dat.css.ps.lim') into another 
### ('base.css.ps.lim'). Only 3 particle size classes allowed. See 
### TT.text.transf.X for transformation involving more than 3 
### particle classes. 'tri.data' may contain other variables 
### (not in 'css.names'). They are returned unchanged with the 
### transformed texture data.

    tri.data,  
    base.css.ps.lim,  
    dat.css.ps.lim,  
    css.names       = NULL,  
    blr.tx          = NULL,  
    text.sum        = NULL,  
    text.tol        = NULL,  
    tri.sum.tst     = NULL,  
    tri.pos.tst     = NULL,  
    trsf.add.opt1   = NULL,   # unused here (but required) 
    trsf.add.opt2   = NULL    # unused here (but required) 
){  #
    #
    TT.auto.set( set.par = FALSE ) 
    # 
    TT.data.test( 
        tri.data    = tri.data, 
        css.names   = css.names, 
        text.sum    = text.sum, 
        text.tol    = text.tol, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst  
    )   #
    #
    blr.tx <- TT.blr.tx.check( 
        blr.tx      = blr.tx, 
        css.names   = css.names
    )   # 
    # 
    # Saving the other columns of the data.frame:
    old.col.nm  <- colnames(tri.data) 
    other.data  <- as.data.frame( 
        tri.data[,!(old.col.nm %in% css.names)] 
    )   #
    colnames(other.data)    <- old.col.nm[!(old.col.nm %in% css.names)] 
    #
    # Taking only texture data:
    colnames(tri.data[,blr.tx]) <- names(blr.tx) 
    #
    # Giving them international names:
    tri.data <- tri.data[,css.names] 
    #
    tri.data <- t(  apply( 
        X       = tri.data, 
        MARGIN  = 1, 
        FUN     = function(X){ 
            cumsum(X) 
        }   #
    )   )   #
    #
    TT.check.ps.lim( 
        base.ps.lim     = base.css.ps.lim,  
        dat.ps.lim      = dat.css.ps.lim,  
        ps.lim.length   = c(4,4) 
    )   #
    #
    if( base.css.ps.lim[1] != 0 ) 
    {   #
        tri.data <- cbind( 
            "ZERO" = rep(0,dim(tri.data)[1]), 
            tri.data  
        )   #
        #
        ps.start <- 1 
    }else{ 
        ps.start <- 2 
    }   #
    #
    base.css.ps.lim2 <- TT.dia2phi(base.css.ps.lim) 
    dat.css.ps.lim2  <- TT.dia2phi(dat.css.ps.lim) 
    #
    tri.data <- t(  apply( 
        X       = tri.data, 
        MARGIN  = 1, 
        FUN     = function(X,base.css.ps.lim2,dat.css.ps.lim2){ 
            c( X[1], diff( approx( 
                x       = dat.css.ps.lim2[ ps.start:4 ], 
                y       = X, 
                xout    = base.css.ps.lim2[ ps.start:4 ], 
                method  = "linear", 
                rule    = 1, 
                ties    = function(...){ 
                    stop("error in TT.text.transf: Unexpected ties in text.cum = f(phi)")
                }   #
            )$"y"   )   )   #
        },  #
        base.css.ps.lim2, 
        dat.css.ps.lim2  
    )   )   #
    #
    if( base.css.ps.lim[1] != 0 ) 
    {   #
        tri.data <- tri.data[,-1] 
    }   #
    #
    tri.data            <- as.data.frame(tri.data) 
    colnames(tri.data)  <- c("CLAY","SILT","SAND") 
    tri.data            <- tri.data[,names(blr.tx)] 
    colnames(tri.data)  <- blr.tx 
    tri.data            <- cbind( 
        tri.data, 
        other.data  
    )   #
    tri.data            <- tri.data[,old.col.nm] 
    #
    return( tri.data ) 
}   #

#     my.text <- data.frame( 
#         "CLAY"  = c(05,60,15,05,25,05,25,45,65,75,13,47), 
#         "SILT"  = c(05,08,15,25,55,85,65,45,15,15,17,43), 
#         "SAND"  = c(90,32,70,70,20,10,10,10,20,10,70,10), 
#         "OC"    = c(20,14,15,05,12,15,07,21,25,30,05,28)  
#     )   #
#     my.text 
#     TT.text.transf( 
#       tri.data        = my.text, 
#       base.css.ps.lim = c(0,2,50,2000), 
#       dat.css.ps.lim  = c(0,2,50,2000)  
#     )   #
#     TT.text.transf( 
#       tri.data        = my.text, 
#       base.css.ps.lim = c(0,2,50,2000), 
#       dat.css.ps.lim  = c(0,2,60,2000)  
#     )   #
#     tmp <- TT.text.transf( 
#       tri.data        = my.text, 
#       base.css.ps.lim = c(1,2,50,2000), 
#       dat.css.ps.lim  = c(1,1.5,60,2000)  
#     )   #
#     tmp ; all( rowSums(tmp[,1:3]) == 100 ) 


#     my.text2 <- my.text

#     my.text2 <- do.call( 
#         what    = "rbind", 
#         args    = lapply( 
#             X   = 1:8500, 
#             FUN = function(X){ 
#                 my.text
#             }   #
#         )   #
#     )   #

#     dim(my.text2)

#     system.time( 
#         TT.text.transf( 
#             tri.data        = my.text2, 
#             base.css.ps.lim = c(0,2,50,2000), 
#             dat.css.ps.lim  = c(0,2,60,2000)  
#         )   #
#     )   #
#     # 4450 texture transformed per second






TT.text.transf.X <- function(# Log-linear transformation of a soil texture data table between 2 particle size systems (X classes).
### Log-linear transformation of a soil texture data table 
### ('tri.data') from one 
### particle size system ('dat.css.ps.lim') into another 
### ('base.css.ps.lim'). No limit in the number of partile size classes 
### in the inputed and outputed texture tables. See TT.text.transf 
### for transformation involving only 3 particle classes. 'tri.data' 
### can only contain texture data.

    tri.data,  
    base.ps.lim,  
    dat.ps.lim,  
    text.sum        = NULL,  
    text.tol        = NULL,  
    tri.sum.tst     = NULL,  
    tri.pos.tst     = NULL   
){  #
    TT.auto.set( set.par = FALSE ) 
    # 
    TT.data.test.X( 
        tri.data    = tri.data, 
        text.sum    = text.sum, 
        text.tol    = text.tol, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst  
    )   #
    #
    tri.data <- t(  apply( 
        X       = tri.data, 
        MARGIN  = 1, 
        FUN     = function(X){ 
            cumsum(X) 
        }   #
    )   )   #
    #
    ps.end   <- dim( tri.data )[2] + 1
    #
    TT.check.ps.lim( 
        base.ps.lim     = base.ps.lim,  
        dat.ps.lim      = dat.ps.lim,  
        ps.lim.length   = c(length(base.ps.lim),ps.end) 
    )   #
    #
    if( base.ps.lim[1] != 0 ) 
    {   #
        tri.data <- cbind( 
            "ZERO" = rep(0,dim(tri.data)[1]), 
            tri.data  
        )   #
        #
        ps.start    <- 1 
    }else{ 
        ps.start    <- 2 
    }   #
    #
    base.ps.lim2 <- TT.dia2phi(base.ps.lim) 
    dat.ps.lim2  <- TT.dia2phi(dat.ps.lim) 
    #
    old.col.nm   <- colnames( tri.data ) 
    #
    tri.data <- t(  apply( 
        X       = tri.data, 
        MARGIN  = 1, 
        FUN     = function(X,base.ps.lim2,dat.ps.lim2){ 
            c( X[1], diff( approx( 
                x       = dat.ps.lim2[ ps.start:ps.end ], 
                y       = X, 
                xout    = base.ps.lim2[ ps.start:length(base.ps.lim) ], 
                method  = "linear", 
                rule    = 1, 
                ties    = function(...){ 
                    stop("error in TT.text.transf: Unexpected ties in text.cum = f(phi)")
                }   #
            )$"y"   )   )   #
        },  #
        base.ps.lim2, 
        dat.ps.lim2  
    )   )   #
    #
    if( base.ps.lim[1] != 0 ) 
    {   #
        tri.data <- tri.data[,-1] 
    }   #
    #
    tri.data            <- as.data.frame(tri.data) 
    colnames(tri.data)  <- paste(sep="","C",1:dim(tri.data)[2]) 
    #
    return( tri.data )
}   #

#     my.text4 <- data.frame( 
#         "CLAY"  = c(05,60,15,05,25,05,25,45,65,75,13,47), 
#         "FSILT" = c(02,04,10,15,25,40,35,20,10,05,10,20), 
#         "CSILT" = c(03,04,05,10,30,45,30,25,05,10,07,23), 
#         "SAND"  = c(90,32,70,70,20,10,10,10,20,10,70,10)  
#     )   #
#     my.text4 
#     TT.text.transf.X( 
#       tri.data    = my.text4, 
#       base.ps.lim = c(0,2,20,50,2000), 
#       dat.ps.lim  = c(0,2,20,50,2000)  
#     )   #
#     TT.text.transf.X( 
#       tri.data    = my.text4, 
#       base.ps.lim = c(0,2,20,50,2000), 
#       dat.ps.lim  = c(0,2,30,60,2000)  
#     )   #
#     tmp <- TT.text.transf.X( 
#       tri.data    = my.text4, 
#       base.ps.lim = c(0,2,50,2000), 
#       dat.ps.lim  = c(0,2,30,60,2000)  
#     )   #
#     tmp ; all( rowSums( tmp ) == 100 ) 






TT.deg2rad <- function(# Internal. Function to convert angle in degree to angle in radian.
### Function to convert angle in degree to angle in radian.
##keywords<< internal

 A
### Angle in Degrees

){  #
    (pi/180)*A
}   #






TT.ifelse <- function(# Internal. Flexible version of ifelse. 
### Flexible version of ifelse. 
##keywords<< internal

 test,
 yes,
 no
 
){  #
    if(test){ res <- yes }else{ res <- no } 
    return(res)
}   #






TT.switch <- function(# Internal. Used in the plot axis drawings.
### Used in the plot axis drawings.
##keywords<< internal

    blr.clock   = TT.get("blr.clock"), 
    c1          = NA, 
    c2          = NA, 
    c3          = NA, 
    c4          = NA, 
    blr.order   = c(1,3,2) 
){  #
    TT.ifelse( 
        "test"  = blr.clock[ blr.order[1] ], 
        "yes"   = TT.ifelse( 
            "test"  = blr.clock[ blr.order[2] ], 
            "yes"   = c1,   # Side 1 == clock and Side 2 == clock
            "no"    = c2    # Side 1 == clock and Side 2 == Aclock
        ),  # 
        "no"    = TT.ifelse( 
            "test"  = blr.clock[ blr.order[3] ], 
            "yes"   = c3,   # Side 1 == Aclock and Side 3 == clock
            "no"    = c4    # Side 1 == Aclock and Side 3 == Aclock
        )   #
    )   #
}   #






TT.css2xy <- function(# Internal. Converts texture data (3 classes) into x-y coordinates. 
### Converts texture data (3 classes) into x-y coordinates. This 
### function is the 'heart' of most soiltexture plot functions.
##keywords<< internal

    tri.data, 
    geo, 
    css.names       = NULL, 
    #
    text.tol        = NULL, 
    #
    tri.sum.tst     = NULL, 
    tri.pos.tst     = NULL, 
    set.par         = FALSE, 
    text.sum        = NULL, 
    blr.clock       = NULL  
){  #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    # Set rest of variables:
    TT.auto.set(set.par=set.par) 
    #
    blr.tx <- TT.blr.tx.check( 
        blr.tx      = blr.tx, 
        css.names   = css.names
    )   # 
    #
    # Check for tlr.an: data type
    if( any( is.null(tlr.an) | is.na(tlr.an) ) | !is.numeric(tlr.an) | (length(tlr.an) != 3) )
    {   #
        stop( paste( 
            sep = "", 
            "tlr.an (=", 
            paste(tlr.an,collapse=";"), 
            ") must be a numeric, non-null, non-na vector of length 3" 
        )   )  
    }   #
    #
    # Check for tlr.an: angle sum must be 180 degrees
    if( sum(tlr.an) != 180 ) 
    {   #
        stop( paste( 
            sep = "", 
            "sum(tlr.an) (=", 
            paste(tlr.an,collapse=";"), 
            ") must be 180 (degrees)"  
        )   )  
    }   #
    #
    # Test the data provided:
    TT.data.test( 
        tri.data    = tri.data, 
        css.names   = css.names, 
        text.sum    = text.sum, 
        text.tol    = text.tol, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst  
        #set.par     = set.par  
    )   #
    # 
    # Test (anti)clock settings
    ok.clock <- list( 
        #       #    Bottom Left    Right  
        "TTT"   = c( TRUE,  TRUE,   TRUE    ), 
        "TXF"   = c( TRUE,  NA,     FALSE   ), 
        "FTX"   = c( FALSE, TRUE,   NA      ), 
        "FFF"   = c( FALSE, FALSE,  FALSE   )  
    )   #
    #
    ok.clock <- unlist( lapply( 
            X           = ok.clock, 
            FUN         = function(X,blr.clock){ 
                identical(blr.clock,X) 
            },  #
            blr.clock   = blr.clock
    )   )   #
    #
    if( !any(ok.clock) )
    {   #
        stop( paste( 
                sep = "", 
                "blr.clock (=", 
                paste(as.character(blr.clock),collapse=";"), 
                ") MUST be one of: ", 
                paste(names(ok.clock),collapse=";"), 
                "; [with X == NA]. consider revising" 
        )   )   #
    }   #
    #
    # Angle transformation: degree to radian
    tlr.an <- TT.deg2rad(tlr.an)
    #
    # # "reverse" the bottom and right orientation to fit x and y orientation
    rev.dt <- function( 
        i, 
        blr.c   = blr.clock, 
        tri.d   = tri.data, 
        blr.t   = blr.tx, 
        text.s  = text.sum  
    ){  #
        val <- tri.d[  , blr.t[i] ]
        if( !is.na(blr.c[i]) )  # Do not reverse NA sides
        {   #
            if( (blr.c[i] & (i != 2)) | (!blr.c[i] & (i == 2)) ) 
            {   #
                val <- ( text.s - val )
            }   #
        }   #
        val 
    }   #
    #
    for( j in 1:3 ){ tri.data[,blr.tx[j]] <- rev.dt("i"=j) }
    #
    # The x,y coordnates calculation is 1st separated depending on blr.clock[2]
    if( !is.na(blr.clock[2]) ){ cond2 <- blr.clock[2] }else{ cond2 <- FALSE } 
    #
    if( cond2 )
    {   #
        ypos    <- tri.data[  , blr.tx[2] ] * sin(tlr.an[2])
    }else{ 
        ypos    <- tri.data[  , blr.tx[3] ] * sin(tlr.an[3])
    }   #
    #
    if( blr.clock[1] )
    {   # if cond2 this is the TTT case, else (!cond2) this is the TXF case.
        xpos    <- tri.data[  , blr.tx[1] ] - ypos/tan(tlr.an[3]) 
    }else{ 
        # if cond2 this is the FTX case, else (cond2) this is the FFF case.
        xpos    <- tri.data[  , blr.tx[1] ] + ypos/tan(tlr.an[2]) 
    }   #
    #
    return( data.frame( "xpos" = xpos , "ypos" = ypos ) )
}   #






TT.points <- function(# Plot a soil texture data table as points on an existing texture plot. 
### Plot a soil texture data table as points on an existing texture plot. 

    tri.data, 
    geo, 
    css.names       = NULL, 
    z.name          = NULL, 
    base.css.ps.lim = NULL, 
    dat.css.ps.lim  = NULL, 
    css.transf      = NULL, 
    text.transf.fun = NULL, 
    trsf.add.opt1   = NULL,   # Additionnal option 1 
    trsf.add.opt2   = NULL,   # Additionnal option 2 
    text.tol        = NULL,
    #
    pch             = NULL, 
    fg              = NULL, 
    col             = NULL, 
    bg              = NULL, 
    cex             = NULL, 
    lwd             = NULL, 
    points.type     = NULL, 
    #
    tri.sum.tst     = NULL, 
    tri.pos.tst     = NULL, 
    #
    z.type          = NULL, 
    z.col.hue       = NULL, 
    z.cex.range     = NULL, 
    z.pch           = NULL, 
    text.sum        = NULL, 
    blr.clock       = NULL, 
    blr.tx          = NULL  
){  #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    if( any( is.null(dat.css.ps.lim) ) )
    {   #
        dat.css.ps.lim  <- base.css.ps.lim 
    }   #   
    #
    # Set the rest of parameters
    TT.auto.set() 
    #
    # Basic checks
    if( dev.cur() == 1 ) 
    {   #
        stop("Cannot add points unless the TT.plot has been drawn")
    }   #
    # 
    if( css.transf )
    {   #
        text.transf.fun <- get( text.transf.fun ) 
        #
        tri.data <- text.transf.fun( 
            tri.data        = tri.data,  
            base.css.ps.lim = base.css.ps.lim,  
            dat.css.ps.lim  = dat.css.ps.lim,  
            css.names       = css.names,  
            blr.tx          = blr.tx,  
            text.sum        = text.sum,  
            text.tol        = text.tol,  
            tri.sum.tst     = tri.sum.tst,  
            tri.pos.tst     = tri.pos.tst,  
            trsf.add.opt1   = trsf.add.opt1, 
            trsf.add.opt2   = trsf.add.opt2  
        )   #
    }   #   
    #
    xy.coord <- TT.css2xy( 
        tri.data    = tri.data, 
        css.names   = css.names, 
        geo         = geo, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        set.par     = FALSE, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    if( !is.null(z.name) & z.type == "bubble" ) 
    {   #
        z           <- order( tri.data[,z.name], decreasing = TRUE ) 
        #
        xy.coord    <- xy.coord[ z ,  ] 
        #
        z           <- tri.data[ z , z.name ] 
        #
        old.fg.col  <- col 
        #
        night.cols  <-  TT.col2hsv(bg)[,"v"] < 0.5 
        #
        if( night.cols )
        {   #
            points.sat <- TT.str(z,0.25,1.00) 
            points.val <- TT.str(z,0.25,1.00) 
        }else{ 
            points.sat <- TT.str(z,1.00,0.25) 
            points.val <- TT.str(z,1.00,0.25) 
        }   #
        #
        col <- hsv( 
            h   = z.col.hue, 
            s   = points.sat, 
            v   = points.val 
        )   #
        #
        cex <- TT.str(z,z.cex.range[1],z.cex.range[2])
        #
        pch1    <- z.pch[1] # Added 20090617 
        pch2    <- z.pch[2] # Added 20090617 
    }else{ 
        pch1    <- pch 
        pch2    <- pch 
    }   #
    #
    nobs <- dim(xy.coord)[1]
    #
    points( 
        x       = xy.coord$"xpos", 
        y       = xy.coord$"ypos", 
        pch     = pch1, # Added 20090617 
        col     = col, 
        bg      = bg, 
        type    = points.type, 
        cex     = cex, 
        lwd     = lwd  # Added 20090617 
    )   # 
    #
    if( !is.null(z.name) & z.type == "bubble" ) 
    {   #
        points( 
            x       = xy.coord$"xpos", 
            y       = xy.coord$"ypos", 
            pch     = pch2, 
            col     = old.fg.col, 
            bg      = bg, 
            type    = points.type, 
            cex     = cex, 
            lwd     = lwd  # Added 20090617 
        )   # 
    }   #
    #
    invisible( xy.coord ) 
}   # 






TT.text <- function(# Plot text labels for each values of a soil texture data table on an existing texture plot. 
### Plot text labels for each values of a soil texture data table on an existing texture plot. 

    tri.data, 
    geo, 
    labels          = NULL, 
    css.names       = NULL, 
    base.css.ps.lim = NULL, 
    dat.css.ps.lim  = NULL, 
    css.transf      = NULL, 
    text.transf.fun = NULL, 
    trsf.add.opt1   = NULL,   # Additionnal option 1 
    trsf.add.opt2   = NULL,   # Additionnal option 2 
    text.tol        = NULL, 
    text.sum        = NULL, 
    blr.clock       = NULL, 
    #
    fg              = NULL, 
    col             = NULL, 
    cex             = NULL, 
    font            = NULL, 
    family.op       = NULL, 
    adj             = NULL, 
    pos             = NULL, 
    offset          = NULL, 
    #
    tri.sum.tst     = NULL, 
    tri.pos.tst     = NULL, 
    blr.tx          = NULL  
){  #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    if( any( is.null(dat.css.ps.lim) ) )
    {   #
        dat.css.ps.lim  <- base.css.ps.lim 
    }   #   
    #
    if( any( is.null(labels) ) )
    {   #
        labels <- 1:dim(tri.data)[1] 
    }   #
    #
    # Set the rest of parameters
    TT.auto.set() 
    #
    # Basic checks
    if( dev.cur() == 1 ) 
    {   #
        stop("Cannot add points unless the TT.plot has been drawn")
    }   #
    # 
    if( css.transf )
    {   #
        text.transf.fun <- get( text.transf.fun ) 
        #
        tri.data <- text.transf.fun( 
            tri.data        = tri.data,  
            base.css.ps.lim = base.css.ps.lim,  
            dat.css.ps.lim  = dat.css.ps.lim,  
            css.names       = css.names,  
            blr.tx          = blr.tx,  
            text.sum        = text.sum,  
            text.tol        = text.tol,  
            tri.sum.tst     = tri.sum.tst,  
            tri.pos.tst     = tri.pos.tst,  
            trsf.add.opt1   = trsf.add.opt1, 
            trsf.add.opt2   = trsf.add.opt2  
        )   #
    }   #   
    #
    xy.coord <- TT.css2xy( 
        tri.data    = tri.data, 
        css.names   = css.names, 
        geo         = geo, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        set.par     = FALSE, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    nobs <- dim(xy.coord)[1]
    #
    text( 
        x       = xy.coord$"xpos", 
        y       = xy.coord$"ypos", 
        labels  = labels, 
        col     = col, 
        #type   = points.type, 
        cex     = cex, 
        font    = font, 
        family  = family.op, 
        adj     = adj, 
        pos     = pos, 
        offset  = offset  
    )   # 
    #
    invisible( xy.coord ) 
}   # 






TT.baseplot <- function(# Internal. Create an empty plot scene for a texture triangle.
### Create an empty plot where a texture triangle can be drawn with 
### other secondary functions (frame, axis, ...). Also return the 
### 'geo' parameters needed by these secondary functions.
##keywords<< internal

    geo             = NULL, 
    class.sys       = "none",  
    # 
    # "GEO" parameters
    blr.clock       = NULL,  
    tlr.an          = NULL,  
    blr.tx          = NULL,  
    text.sum        = NULL,  
    base.css.ps.lim = NULL, 
    # 
    # DATA TESTS:
    tri.sum.tst     = NULL,  
    tri.pos.tst     = NULL,  
    #
    # ADDITIONAL parameters:
    text.tol        = NULL,  
    unit.ps         = NULL,  
    unit.tx         = NULL,  
    #
    b.lim           = NULL,   # default c(0,1) 
    l.lim           = NULL,   # default c(0,1) 
    # 
    main            = NULL,  
    # 
    new.mar         = NULL,  
    # 
    bg              = NULL,  
    fg              = NULL,  
    col             = NULL,  
    cex.main        = NULL,  
    #
    lang            = NULL   
){  # 
    css.names   <- c("CLAY","SILT","SAND")
    # 
    if( is.null(geo) )
    {   #
        geo <- TT.geo.get( 
            class.sys       = class.sys, 
            blr.clock       = blr.clock, 
            tlr.an          = tlr.an, 
            blr.tx          = blr.tx, 
            text.sum        = text.sum, 
            base.css.ps.lim = base.css.ps.lim  
        )   #
    }   #
    # 
    # Set geographical parameters:
    TT.geo.set( 
        geo     = geo  
        #p.env  = environment()  
    )   #
    # 
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Retrieve classes-system (texture triangle) parameters:
    if( class.sys != "none" ) 
    {   #
        TT.data <- TT.get(class.sys) 
    }   #
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Fix the plot limits:
    if( any(is.null(b.lim)) ){ b.lim <- TT.get("b.lim") * text.sum } 
    if( any(is.null(l.lim)) ){ l.lim <- TT.get("l.lim") * text.sum } 
    r.lim <- text.sum - c( b.lim[1] + l.lim[2], b.lim[1] + l.lim[1] ) 
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Create a "base frame", with coordinates expressed in CLAY SILT SAND
    base.frame  <- data.frame( 
        #   #   S1                          S2                          S3 
        "b" = c(b.lim[1],                   b.lim[2],                   b.lim[1] ), 
        "l" = c(l.lim[2],                   l.lim[1],                   l.lim[1] ), 
        "r" = c(text.sum-b.lim[1]-l.lim[2], text.sum-b.lim[2]-l.lim[1], text.sum-b.lim[1]-l.lim[1] ) 
    )   #
    colnames(base.frame)    <- blr.tx 
    #
    if( is.null(main) )
    {   #
        lang.par    <- TT.get("lang.par") 
        #
        if( is.null(lang) ){ lang <- TT.get("lang") } 
        #
        main        <- lang.par[ lang.par$"lang" == lang , "TT" ] 
        main        <- parse(text=main)[[1]]    # Added 2009/06/27 
        #
        if( class.sys != "none" )
        {   #
            main <- paste( 
                sep = "", 
                main, 
                ": ", 
                TT.data$"main" 
            )   #
        }   #
    }   #   
    #
    # Auto-set parameters that are not in par() 
    TT.auto.set(set.par=FALSE) 
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Convert CLAY SILT SAND coordinates into xy coordinates
    ghost.TT    <- TT.css2xy( 
        tri.data    = base.frame, 
        geo         = geo, 
        css.names   = css.names, 
        text.tol    = text.tol, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        set.par     = FALSE, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Setup new graph margins:
    #
    # # default c(5, 2, 4, 2)
    # # c(bottom, left, top, right)
    if( is.null(new.mar) )
    {   #
        #              c(bottom, left, top, right)
        new.mar     <- c(5.0, 3.5, 3.0, 3.5)+0.1 
        #
        if( if(is.null(main)){FALSE}else{is.na(as.character(main))} )
        {   #
            new.mar[3]  <- 0.1 
        }   #
        if( tlr.an[2] > tlr.an[3] )
        {   #
            new.mar[4] <- 0.1 
        }else{ 
            if( tlr.an[2] < tlr.an[3] )
            {   #
                new.mar[2] <- 0.1 
            }else{  # Equality case
                new.mar[c(2,4)] <- 0.1 
            }   #
        }   #
    }   #
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Setup other graphical parameters:
    par.list    <- list( 
        "mar"   = new.mar, 
        "pty"   = "s", 
        "xpd"   = TRUE, 
        "bg"    = bg, 
        "fg"    = fg  
        #"col"  = col  
    )   #
    #
    par.list    <- par.list[ unlist(lapply(X=par.list,FUN=function(X){!is.null(X)})) ] 
    #
    # Sets new par() values
    do.call( what = "par", args = par.list ) 
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # | Ghost plot to set the limits of |
    # | the graph                       | 
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    plot( 
        x           = ghost.TT$"xpos", 
        y           = ghost.TT$"ypos",
        type        = "n", 
        axes        = FALSE, 
        xlim        = range(ghost.TT$"xpos"), 
        ylim        = min(ghost.TT$"ypos")+c(0,diff(range(ghost.TT$"xpos"))), 
        main        = main, 
        cex.main    = cex.main, 
        xlab        = "", 
        ylab        = ""  
    )   #
    #
    # Return the geo(metrical) specifications of the plot
    return( invisible( geo ) ) 
}   #

#   TT.baseplot() 






TT.edges <- function(# Internal. Plot the edges (bare axis) of a soil texture triangle. 
### Plot the edges (bare axis) of a soil texture triangle. This 
### is not a primary plot function, TT.baseplot() must have been 
### called before (usually inside TT.plot()).
##keywords<< internal

    geo, 
    #
    text.tol        = NULL, 
    text.sum        = NULL, 
    blr.clock       = NULL, 
    #
    col.axis        = NULL, 
    plot.axis       = TRUE,     # plot the axis (not only background)
    frame.bg.col    = NULL, 
    lwd.axis        = NULL, 
    #
    tri.sum.tst     = NULL, 
    tri.pos.tst     = NULL, 
    bg              = NULL  
){  # 
    css.names   <- c("CLAY","SILT","SAND")  
    # 
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic sets remaining NULL varaibles 
    TT.auto.set() 
    #
    blr.tx <- TT.blr.tx.check( 
        blr.tx      = blr.tx, 
        css.names   = css.names  
    )   # 
    #
    # Set the base frame
    base.frame <- TT.get("base.frame") * text.sum 
    colnames(base.frame)    <- blr.tx 
    # # base.frame is not in the options so has not been called at 
    # # this stage...
    #
    if( is.null(frame.bg.col) ) 
    {   #
        frame.bg.col    <- TT.DJ.col( 
            cl      = c(bg,col.axis), 
            w       = c(0.9,0.1), 
            gray.l  = FALSE 
        )   #
    }   # 
    #
    tri.TT <- TT.css2xy( 
        tri.data    = base.frame, 
        geo         = geo, 
        css.names   = css.names, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        set.par     = FALSE, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    xpos    <- tri.TT$"xpos"
    ypos    <- tri.TT$"ypos"
    #
    if( !plot.axis ){ col.axis <- NA }
    #
    polygon( 
        x       = tri.TT$"xpos", 
        y       = tri.TT$"ypos", 
        border  = col.axis, 
        col     = frame.bg.col, 
        lwd     = lwd.axis  
    )   #
}   #






TT.lines <- function(# Internal. Used to plot line elements of a texture plot axis, ticks, arrows, etc.
### Used to plot line elements of a texture plot axis, ticks, arrows, etc.

    geo         = geo, 
    #
    at.1.s      = TT.get("at"),         # Start values of lines on each side of the triangle
    at.2.s      = 1 - TT.get("at"),     # (Start values) For grid lines: reverse value, for ticks, at + ticks.shift
    at.3.s      = 0,                    # (Start values) For grid lines: 0 value, for ticks, 0 - ticks.shift
    at.1.e      = TT.get("at"),         # End values of lines on each side of the triangle
    at.2.e      = 0,                    # at.2.e: (End values) logically equal to at.3.s
    at.3.e      = 1 - TT.get("at"),     # at.3.e: (End values) logically equal to at.2.s
    #
    text.tol    = NULL, 
    text.sum    = NULL, 
    blr.clock   = NULL, 
    #
    tri.sum.tst = NULL, 
    tri.pos.tst = NULL  
){  #
    css.names   <- c("CLAY","SILT","SAND") 
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    TT.auto.set() 
    #
    blr.tx <- TT.blr.tx.check( 
        blr.tx      = blr.tx, 
        css.names   = css.names
    )   # 
    #
    # If one of the parameter has length one, expand it to at.1.s length:
    # + adapt the scale to text.sum
    for( j in c("at.1.s","at.2.s","at.3.s","at.1.e","at.2.e","at.3.e") )
    {   #
        if(  (length(get(j)) == 1) & (length(at.1.s) != 1)  )
        {   #
            assign( 
                x       = j, 
                value   = rep(get(j),length(at.1.s))
            )   #
        }   # 
        assign( 
            x       = j, 
            value   = text.sum * get(j)
        )   #
    }   #
    #
    # NEW NEW NEW NEW
    grid.lns <- list( 
        "B" = list( # Dataframe of BLR coordinates for grid lines starting from the bottom of the triangle:
            "start" = data.frame(   # BLR (or CSS) values of the segments starts
                #     #                     # TTT   # TXF   # FTX   # FFF   
                "B" = at.1.s, 
                "L" = TT.switch( blr.clock, at.3.s, at.2.s, at.3.s, at.2.s ), 
                "R" = TT.switch( blr.clock, at.2.s, at.3.s, at.2.s, at.3.s ) 
            ),  #
            "end"   = data.frame(   # BLR (or CSS) values of the segments ends
                "B" = at.1.e, 
                "L" = TT.switch( blr.clock, at.3.e, at.2.e, at.3.e, at.2.e ), 
                "R" = TT.switch( blr.clock, at.2.e, at.3.e, at.2.e, at.3.e ) 
            )   #
        ),  #
        "L" = list( # Dataframe of BLR coordinates for grid lines starting from the left of the triangle:
            "start" = data.frame(   # BLR (or CSS) values of the segments starts
                #     #                     # TTT   # TXF   # FTX   # FFF   
                "B" = TT.switch( blr.clock, at.2.s, at.2.s, at.3.s, at.3.s ), 
                "L" = at.1.s, 
                "R" = TT.switch( blr.clock, at.3.s, at.3.s, at.2.s, at.2.s )  
            ),  #
            "end"   = data.frame(   # BLR (or CSS) values of the segments ends
                "B" = TT.switch( blr.clock, at.2.e, at.2.e, at.3.e, at.3.e ), 
                "L" = at.1.e, 
                "R" = TT.switch( blr.clock, at.3.e, at.3.e, at.2.e, at.2.e )  
            )   #
        ),  #
        "R" = list( # Dataframe of BLR coordinates for grid lines starting from the right of the triangle:
            "start" = data.frame(   # BLR (or CSS) values of the segments starts
                #     #                     # TTT   # TXF   # FTX   # FFF   
                "B" = TT.switch( blr.clock, at.3.s, at.3.s, at.2.s, at.2.s ), 
                "L" = TT.switch( blr.clock, at.2.s, at.2.s, at.3.s, at.3.s ), 
                "R" = at.1.s  
            ),  #
            "end"   = data.frame(   # BLR (or CSS) values of the segments ends
                "B" = TT.switch( blr.clock, at.3.e, at.3.e, at.2.e, at.2.e ), 
                "L" = TT.switch( blr.clock, at.2.e, at.2.e, at.3.e, at.3.e ), 
                "R" = at.1.e  
            )   #
        )   #
    )   #
    #
    grid.lns <- lapply(  
        X   = grid.lns, 
        FUN = function(X){ 
            assign("x",X)
            lapply( 
                X   = x, 
                FUN = function(X){ 
                    colnames(X) <- blr.tx
                    TT.css2xy( 
                        tri.data    = X, 
                        geo         = geo, 
                        css.names   = css.names, 
                        text.tol    = text.tol, 
                        tri.sum.tst = tri.sum.tst, 
                        tri.pos.tst = tri.pos.tst, 
                        set.par     = FALSE, 
                        text.sum    = text.sum, 
                        blr.clock   = blr.clock  
                    )   #
                }   #
            )   #
        }   #
    )   #
    #
    return( grid.lns ) 
}   #






TT.grid <- function(# Plot a grid at regular texture intervals inside an existing soil texture triangle. 
### Plot a grid at regular texture intervals inside an existing soil texture triangle.

    geo             = geo, 
    at              = NULL, 
    #
    text.tol        = NULL, 
    text.sum        = NULL, 
    blr.clock       = NULL, 
    #
    grid.col        = NULL, 
    grid.lty        = NULL, 
    lwd.axis        = NULL, 
    #
    tri.sum.tst     = NULL, 
    tri.pos.tst     = NULL, 
    #
    # Parameters for auto adaptation of the grid color
    # to the class polygon background colors
    class.p.bg.col  = NULL,     # added 2009/05/18 
    class.p.bg.hue  = NULL,     # added 2009/05/18 
    frame.bg.col    = NULL,     # added 2009/05/19 
    bg              = NULL,     # added 2009/05/22 
    col.axis        = NULL      # added 2009/05/22 
){  #
    TT.auto.set() 
    #
    # 1. Setting the colors
    if( is.null(grid.col) ) 
    {   #
        class.p.bg.col.test <- is.logical( class.p.bg.col ) 
        if( class.p.bg.col.test )
        {   #
            class.p.bg.col.test <- class.p.bg.col 
        }else{ 
            class.p.bg.col.test <- TRUE 
        }   #
        #
        # There is a color gradient in the texture classes polygons?
        if( class.p.bg.col.test ) 
        {   # Check "darkness" of the background:
            night.cols  <-  TT.col2hsv(bg)[,"v"] < 0.5 
            #
            grid.col    <- hsv( 
                h   = class.p.bg.hue, 
                # Below, check with class polygon color: consitency:
                s   = ifelse(night.cols,0.45,0.45), # a little less than the min saturation (night or day)
                v   = ifelse(night.cols,0.20,0.80)  # equal the min value night or max value day
            )   # 
        # No color gradient
        }else{ 
            # Frame backgound color is NULL (so default = light gray)
            if( is.null(frame.bg.col) )
            {   #  # Added 20090617 
                # Step that will "remove" transparency:
                bg.hsv      <- col2rgb( bg, alpha = FALSE )[,1]/255 
                #
                grid.col    <- rgb( 
                    red   = bg.hsv["red"], 
                    green = bg.hsv["green"], 
                    blue  = bg.hsv["blue"]  
                )   #
            # Frame backgound color is not NULL 
            }else{ 
                # gid color as a mix of frame background and frame line colors
                grid.col    <- TT.DJ.col( 
                    cl      = as.character( c(frame.bg.col,col.axis) ), 
                    w       = c(0.9,0.1), 
                    gray.l  = FALSE 
                )   #
            }   #
        }   #
    }   #
    #
    at.r    <- 1 - at 
    at.0    <- 0
    #
    grid.lns <- TT.lines( 
        geo         = geo, 
        at.1.s      = at, 
        at.2.s      = at.r, 
        at.3.s      = at.0, 
        at.1.e      = at, 
        at.2.e      = at.0, 
        at.3.e      = at.r, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    invisible( lapply( 
            X   = grid.lns, 
            FUN = function(X){ 
                segments( 
                    x0      = X$"start"$"xpos", 
                    y0      = X$"start"$"ypos", 
                    x1      = X$"end"$"xpos", 
                    y1      = X$"end"$"ypos", 
                    col     = grid.col, 
                    lty     = grid.lty, 
                    lwd     = lwd.axis 
                )   #
            }   #
    )   )   #
    #
    return( invisible( grid.lns ) )
}   #





TT.ticks <- function(# Internal. Plot the axis' ticks of a texture triangle plot. 
### Plot the axis' ticks of a texture triangle plot. 
##keywords<< internal

    geo, 
    at          = NULL, 
    #
    text.tol    = NULL, 
    text.sum    = NULL, 
    blr.clock   = NULL, 
    #
    tk.s        = NULL, 
    #
    tri.sum.tst = NULL, 
    tri.pos.tst = FALSE,    # Ticks are outside the triangle 
    lwd.axis    = NULL, 
    col.axis    = NULL  
){  #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    TT.auto.set() 
    #
    at.2.s  <- 1 - at
    at.3.s  <- 0
    at.2.e  <- at.2.s + tk.s
    at.3.e  <- at.3.s - tk.s
    #
    grid.lns <- TT.lines( 
        geo         = geo, 
        at.1.s      = at, 
        at.2.s      = at.2.s, 
        at.3.s      = at.3.s, 
        at.1.e      = at, 
        at.2.e      = at.2.e, 
        at.3.e      = at.3.e, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    # NEW NEW NEW
    plot.TF <- c( 
        #     #                     # TTT  # TXF  # FTX  # FFF   
        "B" = TT.switch( blr.clock, TRUE,  TRUE,  TRUE,  TRUE ), 
        "L" = TT.switch( blr.clock, TRUE,  FALSE, TRUE,  TRUE ), 
        "R" = TT.switch( blr.clock, TRUE,  TRUE,  FALSE, TRUE )  
    )   #
    # #
    # if( is.null(col.axis) ){ col.axis <- col }
    #
    invisible( lapply( 
            X   = names(grid.lns), 
            FUN = function(X){ 
                if( plot.TF[X] )
                {   #
                    segments( 
                        x0      = grid.lns[[X]]$"start"$"xpos", 
                        y0      = grid.lns[[X]]$"start"$"ypos", 
                        x1      = grid.lns[[X]]$"end"$"xpos", 
                        y1      = grid.lns[[X]]$"end"$"ypos", 
                        col     = col.axis, 
                        lty     = 1, 
                        lwd     = lwd.axis  
                    )   #
                }   #
            }   #
    )   )   #
}   #






TT.ticks.lab <- function(# Internal. Plot the axis ticks' labels of a texture triangle plot. 
### Plot the axis ticks' labels of a texture triangle plot. 
##keywords<< internal

    geo, 
    at          = NULL, 
    #
    text.tol    = NULL, 
    text.sum    = NULL, 
    blr.clock   = NULL, 
    tlr.an      = NULL, 
    #
    tk.ls       = NULL, 
    #
    tri.sum.tst = NULL, 
    tri.pos.tst = FALSE,    # Ticks are outside the triangle
    #
    col.axis    = NULL, 
    font.axis   = NULL, 
    cex.axis    = NULL, 
    family.op   = NULL  
){  #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    TT.auto.set() 
    #
    at.2.s  <- 1 - at
    at.3.s  <- 0
    at.2.e  <- at.2.s + tk.ls
    at.3.e  <- at.3.s - tk.ls
    #
    grid.lns <- TT.lines( 
        geo         = geo, 
        at.1.s      = at, 
        at.2.s      = at.2.s, 
        at.3.s      = at.3.s, 
        at.1.e      = at, 
        at.2.e      = at.2.e, 
        at.3.e      = at.3.e, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    # NEW NEW NEW
    angle   <- c( 
        #     #                     # TTT       # TXF   # FTX   # FFF   
        "B" = TT.switch( blr.clock, -tlr.an[3], +90,    -90,    +tlr.an[2]  ), 
        "L" = TT.switch( blr.clock, +00,        NA,     +00,    -tlr.an[3]  ), 
        "R" = TT.switch( blr.clock, +tlr.an[2], +00,     NA,    +00         )  
    )   #
    #
    # NEW NEW NEW
    plot.TF <- c( 
        #     #                     # TTT  # TXF  # FTX  # FFF   
        "B" = TT.switch( blr.clock, TRUE,  TRUE,  TRUE,  TRUE ), 
        "L" = TT.switch( blr.clock, TRUE,  FALSE, TRUE,  TRUE ), 
        "R" = TT.switch( blr.clock, TRUE,  TRUE,  FALSE, TRUE )  
    )   #
    #
    invisible( lapply( 
            X   = names(grid.lns), 
            FUN = function(X){ 
                if( plot.TF[X] )
                {   #
                    text( 
                        x           = grid.lns[[X]]$"end"$"xpos", 
                        y           = grid.lns[[X]]$"end"$"ypos", 
                        labels      = at * text.sum, 
                        col         = col.axis, 
                        srt         = angle[X], 
                        font        = font.axis, 
                        cex         = cex.axis, 
                        family      = family.op  
                    )   #
                }   #
            }   #
    )   )   #
}   #






TT.axis.arrows <- function(# Internal. Plot the axis' arrows of a texture triangle plot. 
### Plot the axis' arrows of a texture triangle plot. 
##keywords<< internal

    geo, 
    #css.names       = NULL, 
    css.lab         = NULL, 
    a.l             = TT.get("arrows.lims"), 
    a.h.s           = TT.get("arrows.head.shift"), 
    a.t.s           = TT.get("arrows.text.shift"), 
    a.t.s2          = TT.get("arrows.text.shift2"), 
    a.b.s           = TT.get("arrows.base.shift"), 
    # 
    text.tol        = NULL, 
    text.sum        = NULL, 
    blr.clock       = NULL, 
    tlr.an          = NULL, 
    base.css.ps.lim = NULL, 
    # 
    tri.sum.tst     = FALSE,    # !!!! Set to FALSE because sums are not equal to 1 with right-triangles oblic-side
    tri.pos.tst     = FALSE,    # !!!! Set to FALSE because values are outside the triangle
    #
    lwd.lab         = NULL, 
    arrows.lty      = NULL, 
    col.lab         = NULL, 
    font.lab        = NULL, 
    cex.lab         = NULL, 
    family.op       = NULL, 
    unit.ps         = NULL, 
    unit.tx         = NULL, 
    lang            = NULL  
){  #
    css.names   <- TT.get("css.names") 
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    TT.auto.set() 
    #
    blr.tx <- TT.blr.tx.check( 
        blr.tx      = blr.tx, 
        css.names   = css.names
    )   # 
    #
    tri.blr.ps.lim  <- TT.blr.ps.lim( 
        blr.tx      = blr.tx, 
        css.ps.lim  = base.css.ps.lim  
    )   # 
    # 
    lang.par <- TT.get("lang.par") 
    lang.par <- unlist( lang.par[ lang.par$"lang" == lang ,  ] )
    #
    if( any(is.null(css.lab)) )
    {   #
        base.expr   <- expression( bold('[') * UNITa * bold(']') ~ bold(FRACTION) ~ bold(MIN) * bold('-') * bold(MAX) ~ UNITb ) 
        #
        blr.lab     <- c( 
            "B" = eval( substitute( 
                expression( bold('[') * UNITa * bold(']') ~ bold(FRACTION) ~ bold(MIN) * bold('-') * bold(MAX) ~ UNITb ), 
                list( 
                    "FRACTION"  = parse(text=lang.par[names(blr.tx)[1]])[[1]], 
                    "MIN"       = as.character( tri.blr.ps.lim[1,"B"] ), 
                    "MAX"       = as.character( tri.blr.ps.lim[2,"B"] ), 
                    "UNITa"     = unit.tx, 
                    "UNITb"     = unit.ps  
                ) ) #
            ),  #
            "L" = eval( substitute( 
                expression( bold('[') * UNITa * bold(']') ~ bold(FRACTION) ~ bold(MIN) * bold('-') * bold(MAX) ~ UNITb ), 
                list( 
                    "FRACTION"  = parse(text=lang.par[names(blr.tx)[2]])[[1]], 
                    "MIN"       = as.character( tri.blr.ps.lim[1,"L"] ), 
                    "MAX"       = as.character( tri.blr.ps.lim[2,"L"] ), 
                    "UNITa"     = unit.tx, 
                    "UNITb"     = unit.ps  
                ) ) #
            ),  #
            "R" = eval( substitute( 
                expression( bold('[') * UNITa * bold(']') ~ bold(FRACTION) ~ bold(MIN) * bold('-') * bold(MAX) ~ UNITb ), 
                list( 
                    "FRACTION"  = parse(text=lang.par[names(blr.tx)[3]])[[1]], 
                    "MIN"       = as.character( tri.blr.ps.lim[1,"R"] ), 
                    "MAX"       = as.character( tri.blr.ps.lim[2,"R"] ), 
                    "UNITa"     = unit.tx, 
                    "UNITb"     = unit.ps  
                ) ) #
            )   #
        )   #
    }else{ 
        names(css.lab)  <- css.names 
        blr.lab         <- css.lab[ blr.tx ] 
    }   #
    #
    names(blr.lab) <- c("B","L","R") 
    #
    # NEW NEW NEW
    plot.TF <- c( 
        #     #                     # TTT  # TXF  # FTX  # FFF   
        "B" = TT.switch( blr.clock, TRUE,  TRUE,  TRUE,  TRUE ), 
        "L" = TT.switch( blr.clock, TRUE,  FALSE, TRUE,  TRUE ), 
        "R" = TT.switch( blr.clock, TRUE,  TRUE,  FALSE, TRUE )  
    )   #
    #
    # 1. Drawing the arrows base:
    at.1.s  <- a.l[1] 
    at.2.s  <- (1-at.1.s) + a.b.s 
    at.3.s  <- 0 - a.b.s 
    at.1.e  <- a.l[2] 
    at.2.e  <- (1-at.1.e) + a.b.s 
    at.3.e  <- 0 - a.b.s 
    #
    grid.lns1 <- TT.lines( 
        geo         = geo, 
        at.1.s      = at.1.s, 
        at.2.s      = at.2.s, 
        at.3.s      = at.3.s, 
        at.1.e      = at.1.e, 
        at.2.e      = at.2.e, 
        at.3.e      = at.3.e, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    invisible( lapply( 
            X   = names(grid.lns1), 
            FUN = function(X){ 
                if( plot.TF[X] )
                {   #
                    x0  <- grid.lns1[[X]]$"start"$"xpos" 
                    y0  <- grid.lns1[[X]]$"start"$"ypos" 
                    x1  <- grid.lns1[[X]]$"end"$"xpos" 
                    y1  <- grid.lns1[[X]]$"end"$"ypos" 
                    #
                    segments( 
                        x0      = x0,   y0  = y0, 
                        x1      = x1,   y1  = y1, 
                        col     = col.lab, 
                        lty     = arrows.lty, 
                        lwd     = lwd.lab 
                    )   #
                    #
                    # Add a nice little point at the arrow start point:
                    points( 
                        x       = c(x0), 
                        y       = c(y0), 
                        col     = col.lab, 
                        pch     = 16, 
                        cex     = 0.5 
                    )   #
                }   #
            }   #
    )   )   #
    #
    # 2. Drawing the arrows 2nd part and head:
    at.1.s  <- a.l[2] 
    at.2.s  <- (1-at.1.s) + a.b.s 
    at.3.s  <- 0 - a.b.s 
    at.1.e  <- a.l[2] 
    at.2.e  <- (1-at.1.e) + a.h.s 
    at.3.e  <- 0 - a.h.s 
    #
    grid.lns2 <- TT.lines( 
        geo         = geo, 
        at.1.s      = at.1.s, 
        at.2.s      = at.2.s, 
        at.3.s      = at.3.s, 
        at.1.e      = at.1.e, 
        at.2.e      = at.2.e, 
        at.3.e      = at.3.e, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    if( any( !plot.TF ) ) 
    {   #
        at.1.s  <- - a.b.s 
        at.2.s  <- 0.5 + a.b.s 
        at.3.s  <- 0.5 + a.b.s 
        at.1.e  <- - a.h.s 
        at.2.e  <- 0.5 + a.h.s 
        at.3.e  <- 0.5 + a.h.s 
        #
        grid.lnsX <- TT.lines( 
            geo         = geo, 
            at.1.s      = at.1.s, 
            at.2.s      = at.2.s, 
            at.3.s      = at.3.s, 
            at.1.e      = at.1.e, 
            at.2.e      = at.2.e, 
            at.3.e      = at.3.e, 
            tri.sum.tst = tri.sum.tst, 
            tri.pos.tst = tri.pos.tst, 
            text.sum    = text.sum, 
            blr.clock   = blr.clock  
        )   #
        #
        grid.lns2[ !plot.TF ] <- grid.lnsX[ !plot.TF ]
    }   #
    #
    invisible( lapply( 
            X   = names(grid.lns2), 
            FUN = function(X){ 
                arrows( 
                    x0      = grid.lns2[[X]]$"start"$"xpos", 
                    y0      = grid.lns2[[X]]$"start"$"ypos", 
                    x1      = grid.lns2[[X]]$"end"$"xpos", 
                    y1      = grid.lns2[[X]]$"end"$"ypos", 
                    code    = 2, 
                    length  = 0.15, 
                    col     = col.lab, 
                    lty     = arrows.lty, 
                    lwd     = lwd.lab  
                )   #
            }   #
    )   )   #   
    #
    #
    # 3. Drawing the arrows text/label:
    at.1.s  <- a.l[2] + a.t.s2 
    at.2.s  <- (1-at.1.s) + a.t.s 
    at.3.s  <- 0 - a.t.s 
    at.1.e  <- a.l[2] + a.t.s2      # in fact usesell/for compatibility
    at.2.e  <- (1-at.1.e) + a.t.s   # in fact usesell/for compatibility
    at.3.e  <- 0 - a.t.s            # in fact usesell/for compatibility
    #
    grid.lns3 <- TT.lines( 
        geo         = geo, 
        at.1.s      = at.1.s, 
        at.2.s      = at.2.s, 
        at.3.s      = at.3.s, 
        at.1.e      = at.1.e, 
        at.2.e      = at.2.e, 
        at.3.e      = at.3.e, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    if( any( !plot.TF ) ) 
    {   #
        at.1.s  <- -a.b.s 
        at.2.s  <- 0.5 - at.1.s 
        at.3.s  <- 0.5 
        at.1.e  <- -a.b.s 
        at.2.e  <- 0.5 - at.1.e 
        at.3.e  <- 0.5 
        #
        grid.lnsX <- TT.lines( 
            geo         = geo, 
            at.1.s      = at.1.s, 
            at.2.s      = at.2.s, 
            at.3.s      = at.3.s, 
            at.1.e      = at.1.e, 
            at.2.e      = at.2.e, 
            at.3.e      = at.3.e, 
            tri.sum.tst = tri.sum.tst, 
            tri.pos.tst = tri.pos.tst, 
            text.sum    = text.sum, 
            blr.clock   = blr.clock  
        )   #
        #
        grid.lns3[ !plot.TF ] <- grid.lnsX[ !plot.TF ]
    }   #
    # 
    # NEW NEW NEW
    angle   <- c( 
        #     #                     # TTT       # TXF       # FTX       # FFF   
        "B" = TT.switch( blr.clock, +00,        +00,        +00,        +00         ), 
        "L" = TT.switch( blr.clock, +tlr.an[2], +tlr.an[2], +tlr.an[2], +tlr.an[2]  ), 
        "R" = TT.switch( blr.clock, -tlr.an[3], +tlr.an[3], -tlr.an[3], -tlr.an[3]  )  
    )   #
    # 
    # NEW NEW NEW
    adj2    <-  c( 
        #     #                     # TTT   # TXF   # FTX   # FFF   
        "B" = TT.switch( blr.clock, 1,      1,      0,      0   ), 
        "L" = TT.switch( blr.clock, 0,      1,      0,      1   ), 
        "R" = TT.switch( blr.clock, 0,      0,      0,      1   )  
    )   #
    #
    invisible( lapply( 
            X   = names(grid.lns3), 
            FUN = function(X){ 
                text( 
                    x       = grid.lns3[[X]]$"start"$"xpos", 
                    y       = grid.lns3[[X]]$"start"$"ypos", 
                    labels  = blr.lab[X], # Changed 2009/06/30 
                    col     = col.lab, 
                    font    = font.lab, 
                    cex     = cex.lab, 
                    family  = family.op, 
                    srt     = angle[X], 
                    adj     = c(adj2[X],.5), 
                )   #
            }   #
    )   )   #   
}   #



# TT.iwd <- function( 
#     x, 
#     y, 
#     z, 
#     grid.nx     = 10, 
#     grid.ny     = 10, 
#     pow         = 2, 
#     min.dist    = 100 
# ){  #
#     x.rng   <- range(x) 
#     y.rng   <- range(y) 
#     #
#     x.val       <- seq(from=x.rng[1],to=x.rng[2],length.out=grid.nx)
#     y.val       <- seq(from=y.rng[1],to=y.rng[2],length.out=grid.ny)
#     #
#     point.df    <- data.frame(  "x" = x,        "y" = y,        "z" = z     ) 
#     grid.df     <- expand.grid( "x" = x.val,    "y" = y.val,    "z" = NA    ) 
#     #
#     dist.mat    <- dist( 
#         x       = rbind(point.df[,c("x","y")],grid.df[,c("x","y")]), 
#         method  = "euclidean", 
#         diag    = FALSE, 
#         upper   = FALSE  
#     )   #
#     dist.mat    <- as.matrix( dist.mat )
#     dist.mat    <- dist.mat[ (dim(point.df)[1]+1):dim(dist.mat)[1] , 1:dim(point.df)[1] ]
#     dist.mat[ dist.mat > min.dist ] <- NA
#     #
#     grid.df$"z" <- unlist( apply( 
#         X       = dist.mat, 
#         MARGIN  = 1, 
#         FUN     = function(X,z.vec){ 
#             sel <- !is.na(X)
#             if( any(sel) )
#             {   #
#                 w.z <- 1/(X^pow)
#                 sum( z.vec * w.z, na.rm = T )/sum( w.z, na.rm = T ) 
#             }else{ NA }
#         },  #
#         z.vec   = z 
#     )   )   #
#     #
#     return( grid.df ) 
# }   #



# TT.image.shade <- function(
#     iwd.res, 
#     pol.x, 
#     pol.y 
# ){  #
#     require(sp) 
#     #
#     pip <- point.in.polygon(
#         point.x = iwd.res$"x", 
#         point.y = iwd.res$"y", 
#         pol.x   = pol.x, 
#         pol.y   = pol.y 
#     )   #
#     #
#     pip <- pip != 0
#     #
#     iwd.res$"z"[!pip] <- NA
#     #
#     return(iwd.res)
# }   #



# TT.as.image <- function(iwd.res)
# {   #
#     iwd.res <- iwd.res[ order(iwd.res$"x",iwd.res$"y") ,  ] 
#     x   <- unique(iwd.res$"x") 
#     y   <- unique(iwd.res$"y") 
#     z   <- matrix( 
#         data    = iwd.res$"z", 
#         nrow    = length(y), 
#         ncol    = length(x), 
#         byrow   = TRUE 
#     )   #
#     list( 
#         "x" = x, 
#         "y" = y, 
#         "z" = z 
#     )   #
# }   #







TT.dataset <- function(# Genetates a virtual cross correlated clay silt sand + Z dataset. 
### Genetates a virtual cross correlated clay silt sand + Z dataset, 
### where Z is a virtual 4th variable correlated to the texture.

    n, 
    seed.val    = NULL, 
    css.names   = NULL, 
    text.sum    = NULL 
){  #
    TT.auto.set(set.par=FALSE) 
    
    # require( "MASS" ) 
    
    nm      <- c(css.names,"Z") 
    
    CorMat  <- as.matrix(   data.frame( 
            #   #   "1",   "2",   "3",   "Z"
            "1" = c(+1.000,+0.300,-0.600,+0.500), 
            "2" = c(+0.300,+1.000,-0.300,+0.200), 
            "3" = c(-0.600,-0.300,+1.000,-0.500), 
            "Z"     = c(+0.500,+0.200,-0.500,+1.000) 
    )   )   #
    #
    colnames(CorMat)    <- nm 
    rownames(CorMat)    <- nm 
    #
    sd.var              <- rep(6,length(nm))
    names(sd.var)       <- nm
    CovMat              <- CorMat
    for( i in nm )
    {   #
        for( j in nm )
        {   #
            CovMat[i,j] <- CorMat[i,j] * sd.var[i] * sd.var[j] 
        }   #
    }   #
    #
    if(!is.null(seed.val)){ set.seed(seed.val) } 
    #
    rand.text <- mvrnorm(n=n,mu=rep(10,4),Sigma=CovMat)
    #
    rand.text[,1:3] <- t(   apply( 
            X       = rand.text[,1:3], 
            MARGIN  = 1, 
            FUN     = function(X){ 
                X   <- abs(X)
                (X / sum(X))*text.sum
            }   #
    )   )   #
    #
    return( as.data.frame( rand.text ) )
}   #






TT.classes.tbl <- function(# Returns the table of classes of a texture classification system.
### Returns the table of classes of a texture classification system. 
### Returns the classes abbreviations, names and the vertices numbers 
### that defines each class. Use TT.vertices.tbl() to retrieve the 
### clay silt sand coordinates of the triangle classes vertices. 
###  See also TT.vertices.plot().

    class.sys       = "HYPRES.TT", 
    collapse        = NULL 
){  #
    TT.data <- TT.get( class.sys ) 
    #
    if( is.null(collapse) ){ collapse = ", " }
    #
    tbl <- do.call( 
        what    = "rbind", 
        args    = lapply( 
            X   = names(TT.data[["tt.polygons"]]), 
            FUN = function(X){ 
                Y   <- TT.data[["tt.polygons"]][[X]] 
                c( 
                    "abbr"      = X, 
                    "name"      = Y[["name"]], 
                    "points"    = paste( Y[["points"]], collapse = collapse )  
                )   #
            }   #
        )   #
    )   # 
    #
    return( tbl ) 
}   #

#     TT.classes.tbl() 
#     TT.classes.tbl( class.sys = "USDA.TT" ) 






TT.vertices.tbl <- function(# Returns the table of vertices of a texture classification system.
### Returns the table of vertices of a texture classification system. 
### Returns the clay silt sand coordinates of each vertices. Use 
### TT.classes.tbl() to know the vertices that bounds each texture 
### class. See also TT.vertices.plot().

    class.sys       = "HYPRES.TT"  
){  #
    TT.data <- TT.get( class.sys ) 
    #
    TT.data <- TT.data[["tt.points"]] 
    #
    TT.data <- data.frame( 
        "points"    = 1:dim(TT.data)[1], 
        TT.data  
    )   #
    #
    return( TT.data ) 
}   #


#     TT.vertices.tbl() 
#     TT.vertices.tbl( class.sys = "USDA.TT" ) 






TT.vertices.plot <- function(# Internal. Plot the vertices of a texture classification system. 
### Plot the vertices of a texture classification system, on top 
### of an already drawn texture triangle plot. Also plot the 
### vertices numbers. See TT.vertices.tbl() and TT.classes.tbl() 
### for a non graphic, tabular equivalent of the plot.
##keywords<< internal

    geo, 
    class.sys   = "HYPRES.TT", 
    fg          = NULL, 
    col         = NULL, 
    cex         = NULL, 
    font        = NULL, 
    family.op   = NULL, 
    adj         = NULL, 
    pos         = NULL, 
    offset      = NULL, 
    blr.tx      = NULL, 
    text.sum    = NULL, 
    blr.clock   = NULL  
){  #
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    # Set the rest of parameters
    TT.auto.set() 
    # 
    TT.data <- TT.get( class.sys ) 
    #
    TT.data <- TT.data[["tt.points"]] 
    #
    TT.data <- data.frame( 
        "points"    = 1:dim(TT.data)[1], 
        TT.data * geo[["text.sum"]] 
    )   #
    #
    TT.text( 
        geo         = geo, 
        tri.data    = TT.data, 
        labels      = TT.data[,"points"], 
        fg          = fg, 
        col         = col, 
        cex         = cex, 
        font        = font, 
        family.op   = family.op, 
        adj         = adj, 
        pos         = pos, 
        offset      = offset, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock, 
        blr.tx      = blr.tx  
    )   #
    #
    return( invisible( TT.data ) ) 
}   #






TT.polygon.area <- function(# Internal. Determines the area of 1 polygon (in x-y coordinates). 
### Determines the area of 1 non-intersecting polygon (in x-y 
### coordinates). Used by TT.polygon.centroids(). pol.x[1]:pol.y[1] 
### is supposed different from pol.x[n]:pol.y[n] (i.e. the polygon 
### is NOT closed). 
### After "http://local.wasp.uwa.edu.au/~pbourke/geometry/polyarea/
### Calculating The Area And Centroid Of A Polygon. Written by 
### Paul Bourke, July 1988".
##keywords<< internal

 pol.x,
### Vector of numericals. X coordinates of each vertices of the 
### polygon.

 pol.y
### Vector of numericals. Y coordinates of each vertices of the 
### polygon.

){  # Index range:
    i <- 1:length(pol.x)
    #
    # Close the polygon
    pol.x <- c(pol.x,pol.x[1]) 
    pol.y <- c(pol.y,pol.y[1]) 
    #
    # Calculate the area
    A <- 0.5 * sum( pol.x[i] * pol.y[i+1] - pol.x[i+1] * pol.y[i] )
    #
    return( A ) 
    ### Returns a single numerical: area of the polygon.
}   #






TT.polygon.centroids <- function(# Internal. Determines the centroid of 1 polygon (in x-y coordinates). 
### Determines the centroid of 1 non-intersecting polygon (in x-y 
### coordinates). Used to determine the centroid of each texture 
### class in the texture triangle onces its clay silt sand 
### coordinates have been converted to x-y coordinates. pol.x[1]:pol.y[1] 
### is supposed different from pol.x[n]:pol.y[n] (i.e. the polygon 
### is NOT closed). 
### After "http://local.wasp.uwa.edu.au/~pbourke/geometry/polyarea/ 
### Calculating The Area And Centroid Of A Polygon. Written by 
### Paul Bourke, July 1988".
##keywords<< internal

 pol.x,
### Vector of numericals. X coordinates of each vertices of the 
### polygon.

 pol.y
### Vector of numericals. Y coordinates of each vertices of the 
### polygon.

){  #
    i <- 1:length(pol.x)
    #
    # Calculate the area:
    A <- TT.polygon.area( 
        pol.x = pol.x, 
        pol.y = pol.y  
    )   #
    #
    # Close the polygon
    pol.x <- c(pol.x,pol.x[1]) 
    pol.y <- c(pol.y,pol.y[1]) 
    #
    # Calculate the area
    Cx <- (1/(6*A)) * sum( (pol.x[i]+pol.x[i+1]) * (pol.x[i]*pol.y[i+1] - pol.x[i+1]*pol.y[i]) )  
    Cy <- (1/(6*A)) * sum( (pol.y[i]+pol.y[i+1]) * (pol.x[i]*pol.y[i+1] - pol.x[i+1]*pol.y[i]) )  
    #
    return( c("x"=Cx,"y"=Cy) ) 
    ### Returns a vector of 2 numericals: x and y coordinates of 
    ### the polygon's centroid. Vector items are names "x" and "y". 
}   #





TT.classes <- function(# Plot the texture classes polygons in a texture triangle plot.
### Plot the texture classes ploygons in an existing texture 
### triangle plot. Draw the polygons and the labels inside each 
### polygons.

    geo, 
    class.sys, 
    tri.css.ps.lim  = NULL, 
    css.transf      = NULL, 
    text.transf.fun = NULL, 
    trsf.add.opt1   = NULL,   # Additionnal option 1 
    trsf.add.opt2   = NULL,   # Additionnal option 2 
    text.tol        = NULL, 
    text.sum        = NULL, 
    base.css.ps.lim = NULL, 
    blr.tx          = NULL, 
    blr.clock       = NULL, 
    tri.sum.tst     = NULL, 
    tri.pos.tst     = NULL, 
    bg              = NULL, 
    class.lab.col   = NULL, 
    class.p.bg.col  = NULL, 
    class.p.bg.hue  = NULL, 
    class.line.col  = NULL, 
    class.lty       = NULL, 
    class.lab.show  = NULL, 
    cex.lab         = NULL, 
    font.lab        = NULL, 
    family.op       = NULL, 
    lwd.axis        = NULL, 
    col.axis        = NULL, # NB: not the color of the polygon line. par("col.axis")
 new.centroid=TRUE  
### Single logical. If TRUE (default) the new method (Paul Bourke) 
### is used to calculate the centroid. If FALSE the centroid is 
### taken as the mean x and y coordinates of the vertices.

){  # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Retrieve classes-system (texture triangle) parameters:
    TT.data     <- TT.get(class.sys) 
    # 
    css.names   <- c("CLAY","SILT","SAND") 
    #
    if( any( is.null(tri.css.ps.lim) ) )
    {   #
        tri.css.ps.lim  <- TT.data[["tri.css.ps.lim"]] 
    }   #
    #
    # Set geographical parameters:
    TT.geo.set( 
        geo     = geo  
        #p.env  = environment()  
    )   #
    # 
    # Set remaining parameters:
    TT.auto.set() 
    # 
    tri.data    <- TT.data$"tt.points" * text.sum 
    # 
    if( css.transf )
    {   #
        text.transf.fun <- get( text.transf.fun ) 
        #
        tri.data <- text.transf.fun( 
            tri.data        = tri.data,  
            base.css.ps.lim = base.css.ps.lim,  
            dat.css.ps.lim  = tri.css.ps.lim,  
            css.names       = css.names,  
            blr.tx          = blr.tx,  
            text.sum        = text.sum,  
            text.tol        = text.tol,  
            tri.sum.tst     = tri.sum.tst,  
            tri.pos.tst     = tri.pos.tst,  
            trsf.add.opt1   = trsf.add.opt1,  
            trsf.add.opt2   = trsf.add.opt2   
        )   #
    }   #
    # 
    # CLASES-POLYGONS: background plot
    xy.new <- TT.css2xy( 
        tri.data    = tri.data, 
        geo         = geo, 
        css.names   = css.names, 
        text.tol    = text.tol, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        set.par     = FALSE, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    poly.nm <- names( TT.data$"tt.polygons" )
    #
    # Class centroid computation: used for polygon labels position and background colors
    if( new.centroid ) 
    {   #
        cent.xy <- do.call(# Changed the 2010/06/10
            "cbind", 
            lapply( 
                X   = poly.nm, 
                FUN = function(X){ 
                    sel.vec <- (TT.data$"tt.polygons"[[ X ]])$"points"
                    #
                    TT.polygon.centroids(
                        pol.x = xy.new[ sel.vec , "xpos" ], 
                        pol.y = xy.new[ sel.vec , "ypos" ] 
                    )   #
                }   #
            )   #
        )   #
    }else{ 
        cent.xy <- do.call( 
            "cbind", 
            lapply( 
                X   = poly.nm, 
                FUN = function(X){ 
                    sel.vec <- (TT.data$"tt.polygons"[[ X ]])$"points"
                    cent.x  <- mean( xy.new[ sel.vec , "xpos" ] ) 
                    cent.y  <- mean( xy.new[ sel.vec , "ypos" ] ) 
                    return(c("x"=cent.x,"y"=cent.y))
                }   #
            )   #
        )   #
    }   #
    colnames(cent.xy) <- poly.nm 
    #
    # Set the "night colors" parameter
    night.cols  <-  TT.col2hsv(bg)[,"v"] < 0.5 
    #
    # - Test the type of class.p.bg.col parameter:
    class.p.bg.col.test <- is.logical( class.p.bg.col ) 
    if( class.p.bg.col.test )
    {   #
        class.p.bg.col.test <- class.p.bg.col 
    }else{ 
        class.p.bg.col.test <- TRUE 
    }   #
    #
    for( i.pol in 1:length(poly.nm) )
    {   #
        pol <- poly.nm[i.pol]
        #
        sel.vec <- (TT.data$"tt.polygons"[[ pol ]])$"points"
        #
        # Compute a classes-polygon background HSV color range:
        if( class.p.bg.col.test )
        {   #
            if( !is.character( class.p.bg.col ) )
            {   #
                x.range <- range(cent.xy["x",]) 
                y.range <- range(cent.xy["y",]) 
                a.x     <- (1-0)/(diff(x.range)) 
                b.x     <- 1-a.x*x.range[2]
                a.y     <- (1-0)/(diff(y.range)) 
                b.y     <- 1-a.y*y.range[2]
                #
                # Below: check consitency with TT.grid, if changes
                if( night.cols )
                {   #
                    class.sat <- 0.95 - (cent.xy["x",pol]*a.x+b.x)*0.45 # range 0.50;0.95 
                    class.val <- 0.35 - (cent.xy["y",pol]*a.y+b.y)*0.15 # range 0.20;0.35 
                }else{ 
                    class.sat <- 0.50 + (cent.xy["x",pol]*a.x+b.x)*0.45 # range 0.50;0.95 
                    class.val <- 0.85 + (cent.xy["y",pol]*a.y+b.y)*0.15 # range 0.85;1.00 
                }   #
                #
                class.p.bg.col2 <- hsv( 
                    h = class.p.bg.hue,  
                    s = class.sat,  # max = 0.9 
                    v = class.val   # max = 0.9 
                )   #
            }else{ 
                if( length(class.p.bg.col) != 1 ) 
                {   #
                    class.p.bg.col2 <- class.p.bg.col[i.pol] 
                }else{ 
                    class.p.bg.col2 <- class.p.bg.col 
                }   #
            }   #
        }else{ 
            class.p.bg.col2 <- NA 
        }   #
        #
        # Compute the polygon border colors
        if( is.null(class.line.col) ) 
        {   #
            if( night.cols )
            {   #
                wg  <- c(0.60,0.40) 
            }else{ 
                wg  <- c(0.65,0.35) 
            }   #
            #
            class.line.col <- TT.DJ.col( 
                cl      = c(col.axis,bg), 
                w       = wg, 
                gray.l  = FALSE 
            )   #
        }   #
        #
        if( is.null( class.lty ) )
        {   #
            class.lty   <- par("lty") 
        }   #
        #
        # Plot the classes-polygons: phase 1, the background
        #   the lines are drawn after the grid
        polygon( 
            x       = xy.new[ sel.vec , "xpos" ], 
            y       = xy.new[ sel.vec , "ypos" ], 
            border  = class.line.col, 
            lwd     = lwd.axis, 
            col     = class.p.bg.col2, 
            lty     = class.lty  
        )   #
        #
        #
        # Plot the classes-polygons: phase 3, the polygon labels
        if( class.lab.show != "none" )
        {   #
            if( class.lab.show == "full" )
            {   #
                class.lab <- (TT.data$"tt.polygons"[[ pol ]])$"name"
            }else{ 
                if( class.lab.show == "abr" )
                {   #
                    class.lab <- pol
                }else{ stop("class.lab.show must be 'none', 'abr' or 'full'") }
            }   #
            #
            if( is.null(class.lab.col) )
            {   #
                class.lab.col <- class.line.col 
            }   #
            #
            if( "class.lab.cex" %in% names(TT.data) ){ # New 2015-03-06
                class.lab.cex <- TT.data[[ "class.lab.cex" ]]
            }else{ 
                class.lab.cex <- 1
            }   
            #
            text( 
                x       = cent.xy["x",pol], 
                y       = cent.xy["y",pol], 
                labels  = class.lab, 
                adj     = c(0.5,0.5), 
                offset  = 0, 
                col     = class.lab.col, 
                cex     = cex.lab * class.lab.cex, 
                font    = font.lab, 
                family  = family.op  
            )   #
        }   #
    }   #
}   #

#   geo <- TT.baseplot() 
#   TT.classes(geo=geo,class.sys="USDA.TT")

#     TT.classes(
#         geo             = geo, 
#         class.sys       = "FR.AISNE.TT", 
#         tri.css.ps.lim  = c(0,2,63,2000), 
#         css.transf      = TRUE, 
#         # Additional "graphical" options
#         class.line.col  = "blue", 
#         class.lab.col   = "blue", 
#         lwd.axis        = 2, 
#         class.lab.show  = FALSE, 
#     )   #






TT.plot <- function(# Plot soil texture triangles / diagrams.
### Plot a soil texture triangle (also called soil texture 
### diagrams, or soil texture ternary plots), with or without 
### background soil texture classes boundaries, and with or without
### soil texture data points. The triangle geometry depends on the 
### soil texture classification system chosen ('class.sys' argument) 
### or on 'forcing' parameters (see below). 
### Both the boundaries of the background texture classification 
### system and the texture data points can be transformed from 
### one particle size limits system to another (the particle size 
### limits system of the plot). Default behaviour is no transformation 
### (set 'css.transf' argument to TRUE to allow transformation). 
### There are 3 different way to set the triangle geometry and 
### characteristics (1) setting the 'class.sys' argument [lowest 
### priority], (2) changing one or several values of the 'geo' 
### list of arguments or (3) setting the corresponding arguments 
### of TT.plot() [highest priority]. These arguments are "blr.clock", 
### "tlr.an", "blr.tx", "text.sum", and "base.css.ps.lim". Different 
### geometry arguments can be set at different levels (1, 2 or 3). 
### Case (1) should be used when one wants to use the 'default' triangle 
### geometry associated with a given texture classification system 
### (chosen with the 'class.sys' argument). Case (2) should be used 
### when TT.plot() has been called previously, with a call like 
### geo <- TT.plot(), so the 'geo' object returned can be used 
### for setting the geometry of a new texture triangle TT.plot( 
### geo = geo ) identical to the previous one. Case (3) should be 
### used whenever the user wants to set the geometry of a texture 
### triangle plot different from default values of the texture 
### classification system chosen, and without re-using the geometry 
### from a previous plot. 
### ON DEFAULT VALUES OF TT.plot() ARGUMENTS? As TT.plot() shares 
### its arguments with many other functions, their default value 
### is not defined in TT.plot() source code, but rather in 
### a dedicated list object called 'TT.par' and stored in the 
### environment TT.env. The function TT.get() is used to retrieve 
### the default value of the arguments defined in TT.par (see 
### ?TT.get). For instance, to know the default value of 'class.sys', 
### you can type TT.get("class.sys"). To set a different default 
### value for a given argument in R, use TT.set() (see ?TT.set). 
### For instance to change the default value of 'class.sys', type 
### TT.set( "class.sys" = "USDA.TT" ).

# GENERAL Parameters:

 geo=NULL,
### List. 'geo' is one of the 3 way to set the texture triangle 
### geometry. See there description and hierarchy in the function 
### description. If geo != NULL, then geo must be 
### a list containing 1 or several of the following items: 
### "blr.clock", "tlr.an", "blr.tx", "text.sum", and "base.css.ps.lim". 
### See the options with the same name for a description of the 
### expected values and effects. The list can be created manually 
### (like list( "text.sum" = 1000 ) ), or taken from the output of 
### a previous call to TT.plot(), TT.baseplot() or TT.geo.get() 
### (that return a 'geo' list).

 tri.data=NULL,
### Data frame. Data frame containing the CLAY, SILT and SAND 
### 'coordinates' of the texture data points to be plotted on top 
### of the texture triangle and texture class boundaries. The data 
### frame can contain more column than needed (ignored). The data 
### frame must have column named CLAY, SILT and SAND (uppercase, 
### the order has no importance) or named after the 'css.names' 
### argument (alternative names). If 'z.name' argument is not NULL, 
### the data frame must also contain a column named after 'z.name' 
### value. The sum of CLAY, SILT and SAND must be equal to 'text.sum' 
### ('text.tol' determines the error tolerance).

 add=FALSE,
### Single logical. If FALSE, a new plot is created. If FALSE, 
### the plot is added to the existing one.

 css.names=NULL,
### Vector of 3 character strings. Name of the columns in 'tri.data' 
### that contains the CLAY SILT and SAND values, respectively. 
### If NULL, default c("CLAY","SILT","SAND") value is assumed. Not 
### to be confused with 'css.lab' that defines the labels of the 
### CLAY SILT and SAND axes in the plot.

 z.name=NULL,
### Single character string. Name of the column in 'tri.data' 
### that contains the '4th quantitative variable' whose value 
### must be used to define the points expansion factor and 
### color (bubble plot). If NULL, a simple plot is drawn (no 
### 'bubbles')

 main=NULL,
### Single character string or expression. Main title of the plot.

# INPUT DATA description

 blr.tx=NULL,
### Vector of 3 character strings. The 1st, 2nd and 3rd values must 
### be either CLAY, SILT or SAND, and determines the particle size classes 
### associated with the BOTTOM, LEFT and RIGHT axis, respectively. 
### CLAY, SILT and SAND order in the vector is free, but they should 
### all be used one time. The CLAY, SILT and SAND names must appear 
### whatever the corresponding columns names in 'tri.data' (eventually 
### set by 'css.names') and whatever the labels of the axis in the 
### plot (eventually set by 'css.lab') 

 css.lab=NULL,
### Vector of 3 character strings or 3 expressions. The 1st, 2nd 
### and 3rd values must be text strings or expressions, and determines 
### the axes plot labels for the CLAY, SILT and SAND particle size classes, 
### respectively. 'css.lab' values are independent from columns 
### names in 'tri.data' (eventually set by 'css.names') and from 
### whatever the placement of particle size classes on each axis 
### (eventually set by 'blr.tx') 

 text.sum=NULL,
### Single numerical. Sum of the 3 particle size classes for each texture 
### value (fixed). The real sum of the 3 particle size classes in 'tri.data' 
### should be >= text.sum * (1-text.tol) OR  <= text.sum * (1+text.tol), 
### where 'text.tol' is an argument that can be changed. If some 
### of the texture values don't match this requirement, an error 
### occur (function fails) and TT.plot returns a of bad values with 
### their actual particle size classes sum. You can 'normalise' you data 
### table () prior to the use of TT.plot, by using the function 
### TT.normalise.sum(), so all values match the 'text.sum' criteria. 
### See also 'tri.sum.tst' that can be set to FALSE to avoid 
### sum of particle size classes tests.

# blr.psize.lim  = NULL,

 base.css.ps.lim=NULL,
### Vector of 4 numericals. Particle size boundaries (upper and lower) 
### of the 3 particle size classes (CLAY, SILT and SAND, starting from 
### the lower size of CLAY particles, 0, to the upper size of the 
### SAND particles, 2000), in micrometers, FOR THE BASE PLOT. These 
### particles size class limits are the references and all other 
### texture values with different limits will be converted into 
### that reference if (and only if) css.transf == TRUE (not default). 
### If NULL, 'base.css.ps.lim' will be set to the default value of the 
### texture classification system chosen ('class.sys'). The 
### transformation function is set by 'text.transf.fun' and is 
### a log-linear interpolation by default.

 tri.css.ps.lim=NULL,
### Vector of 4 numericals. Particle size boundaries (upper and lower) 
### of the 3 particle size classes (CLAY, SILT and SAND, starting from 
### the lower size of CLAY particles, 0, to the upper size of the 
### SAND particles, 2000), in micrometers, FOR THE TEXTURE TRIANGLE. 
### If not NULL, different from 'base.css.ps.lim', and 
### css.transf == TRUE (not default), then the CLAY SILT and SAND 
### coordinates of the texture triangle will be converted into 
### the 'base.css.ps.lim' reference. If NULL, 'tri.css.ps.lim' will 
### be set to the default value of the texture classification system 
### chosen ('class.sys'). The transformation function is set by 
### 'text.transf.fun' and is a log-linear interpolation by default.

 dat.css.ps.lim=NULL,
### Vector of 4 numericals. Particle size boundaries (upper and lower) 
### of the 3 particle size classes (CLAY, SILT and SAND, starting from 
### the lower size of CLAY particles, 0, to the upper size of the 
### SAND particles, 2000), in micrometers, FOR THE TEXTURE DATA TABLE
### ('tri.data'). If not NULL, different from 'base.css.ps.lim', and 
### css.transf == TRUE (not default), then the CLAY SILT and SAND 
### coordinates of the texture data in tri.data will be converted into 
### the 'base.css.ps.lim' reference. If NULL, 'tri.css.ps.lim' will 
### be set to the default value of the texture classification system 
### chosen ('class.sys'). The transformation function is set by 
### 'text.transf.fun' and is a log-linear interpolation by default.

 css.transf=NULL,
### Single logical. Set to TRUE to transform the texture coordinates 
### of the texture triangle ('class.sys') or the texture data 
### ('tri.data') into the base particle size class limits. 
### See 'base.css.ps.lim' for the base plot particle size class limits, 
### 'tri.css.ps.lim' for the triangle particle size class limits 
### and 'dat.css.ps.lim' for the data table particle size class limits. 
### The transformation function is set by 'text.transf.fun' and 
### is a log-linear interpolation by default. The default value is 
### FALSE, so no transformation is made.

 text.transf.fun=NULL,
### R function with the same argument names and same output as 
### the function TT.text.transf(). 'text.transf.fun' is the function 
### that transform the texture values from one system of particle 
### class size limits to another. Only used if css.transf == TRUE. 
### Default value is text.transf.fun=TT.text.transf. See also 
### 'base.css.ps.lim', 'tri.css.ps.lim' and 'dat.css.ps.lim'.

 trsf.add.opt1=NULL,
### Non pre-defined format. If the user specifies its own texture 
### transformation function in 'text.transf.fun' (not TT.text.transf()), 
### then he can use 'trsf.add.opt1' and 'trsf.add.opt1' as 
### new, additional, argument for his function. So the format of 
### 'trsf.add.opt1' depends on the function defined by the user 
### in 'text.transf.fun'.

 trsf.add.opt2=NULL,
### Non pre-defined format. If the user specifies its own texture 
### transformation function in 'text.transf.fun' (not TT.text.transf()), 
### then he can use 'trsf.add.opt1' and 'trsf.add.opt1' as 
### new, additional, argument for his function. So the format of 
### 'trsf.add.opt1' depends on the function defined by the user 
### in 'text.transf.fun'.

 unit.ps=NULL,
### Single text string or expression. Unit of particle size class 
### limits displayed on the plot (= part of the axis labels). Does 
### not affect the other calculations. Default micrometers.

 unit.tx=NULL,
### Single text string or expression. Unit of particle texture values
### displayed on the plot (= part of the axis labels). Does 
### not affect the other calculations. Default is percentage.

# TRIANGLE FRAME geometry: Bottom, Left and Right side properties

 blr.clock=NULL,
### Vector of logicals, eventually with NA values. Direction of 
### increasing texture values on the BOTTOM, LEFT and RIGHT axis, 
### respectively. A value of TRUE means that the axis direction is 
### clockwise. A value of FALSE means that the axis direction is 
### counterclockwise. A value of NA means that the axis direction 
### is centripetal. Possible combinations are c(T,T,T); c(F,F,F); 
### c(F,T,NA) and c(T,NA,F), for fully clockwise, fully counterclockwise, 
### right centripetal and left centripetal orientations, respectively.

 tlr.an=NULL,
### Vector of numericals. Value - in degrees - of the TOP, LEFT and 
### RIGHT angles of the triangle. Any value between 0 and 90 is possible, 
### but values belonging to 0 or 45 or 60 or 90 give a better graphical 
### rendering.

# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# | General graphical parameters

# FONTS:
 font=NULL,
### Single integer. Not used yet.

 font.axis=NULL,
### Single integer. Same definition as par("font.axis"). Font of 
### the triangle axis's numbering.

 font.lab=NULL,
### Single integer. Same definition as par("font.lab"). Font of 
### the triangle axis's labels.


 font.main=NULL,
### Single integer. Same definition as par("font.main"). Font of 
### the triangle main title.

# COLORS:

 bg=NULL,
### Text string containing an R color code. Background color of the 
### plot (= outside the triangle). See 'frame.bg.col' for the background 
### color inside the triangle frame.

 fg=NULL,
### Text string containing an R color code. DEPRECATED. foreground 
### color of the plot (= point fill color).

 col=NULL,
### Text string containing an R color code. Same definition as par("col"). Color 
### of the points plotted in the triangle.

 col.axis=NULL,
### Text string containing an R color code. Color of the triangle's 
### axis (line and tics) The color of the texture classes boundaries 
### is set by 'class.line.col'.

 col.lab=NULL,
### Text string containing an R color code. Color of the triangle's 
### labels (text) and arrows. The color of the texture classes labels 
### is set by 'class.lab.col'.

 col.main=NULL,
### Text string containing an R color code. Color of the main title.

# CEX:
 cex=NULL,
### Vector of numerical or single numerical. Same definition as par("cex"). 
### Expansion factor for the points plotted on the triangle.

 cex.axis=NULL,
### Single numerical. Same definition as par("cex.axis"). Expansion factor for 
### the axis tics labels (numbering).

 cex.lab=NULL,
### Single numerical. Same definition as par("cex.lab"). Expansion factor for 
### the axis labels AND the texture classes labels.

 cex.main=NULL,
### Single numerical. Same definition as par("cex.main"). Expansion factor for the 
### main title. 

# LWD: 
 lwd=NULL,
### Single numerical. Same definition as par("lwd"). Line width for 
### the graphical elements inside the triangle (points plotted).

 lwd.axis=NULL,
### Single numerical. Same definition as par("lwd.axis"). Line width 
### for the axis lines, tics and the grid lines inside the triangle.

 lwd.lab=NULL,
### Single numerical. Same definition as par("lwd"). Line width for 
### the direction arrows.

# FAMILY: 
 family.op=NULL,
### Single text string. Same definition as par("family"). Font type 
### to be used in the plot text elements (title, labels)

# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# | Specific graphical parameters

# TRIANGLE FRAME parameters
 frame.bg.col=NULL,
### Text string containing an R color code. Background color of the 
### triangle plot (= inside the triangle). See 'bg' for the background 
### color outside the triangle frame.

# GRID LINES parameters:
 at=NULL,
### Vector of numericals. Location of the grid line start points 
### on all 3 triangles axis. At the moment values are identical for 
### all 3 axis, and changes to that parameter have not been tested. 

 grid.show=NULL,
### Single logical. If set to TRUE (the default) a grid is drawn 
### on the background of the texture triangle. Set to FALSE to 
### remove the grid.

 grid.col=NULL,
### Text string containing an R color code. Color of the grid lines. 
### If equal to NULL, then an appropriate color is used. Appropriate 
### means (i) if 'class.p.bg.col' is FALSE (no color gradient in 
### texture class polygons), then grid.col is equal to 'bg' (without
### transparency) unless a color is specified for the triangle 
### frame background ('frame.bg.col'), in which case grid.col is a 
### mix of 'frame.bg.col' and 'col.axis'. (ii) if 'class.p.bg.col' 
### is TRUE, then grid.col is a light or dark color based on 
### 'class.p.bg.hue' (light if 'bg' is dark and dark if 'bg' is light).

 grid.lty=NULL,
### Single numerical. Line type of the grid lines (same types as 
### par("lty")).

# CLASSES POLYGON (texture) parameters:

 class.sys=NULL,
### Single text string. Text code of the texture classification 
### system to be plotted on the background of the texture triangle. 
### That texture classification system will determines the triangle 
### geometry and particle class size system of the plot, unless 
### higher level options are chosen (see the function definition). 
### Possible values are "none" (no classification plotted), "USDA.TT" 
### (USDA texture triangle), "HYPRES.TT" (texture triangle of the 
### European Sil Map), "FR.AISNE.TT" (French texture triangle 
### of the Aisne region soil survey), "FR.GEPPA.TT" (French GEPPA 
### texture triangle), "DE.BK94.TT" (German texture triangle), 
### "UK.SSEW.TT" (Soil Survey of England and Wales), "AU.TT" 
### (Australian texture triangle), "BE.TT" (Belgium texture triangle), 
### "CA.EN.TT" (Canadian texture triangle, with English class abbreviations) and 
### "CA.FR.TT" (Canadian texture triangle, with French class abbreviations).

 class.lab.show=NULL,
### Single text string. If equal to "abr" (default) or "full", labels 
### are drawn inside texture class polygons with their full name 
### ("full") or abbreviated name ("abr"). If equal to "none", no label 
### is drawn.

 class.lab.col=NULL,
### Text string containing an R color code. Color of the text label 
### drawn inside texture class polygons.

 class.line.col=NULL,
### Text string containing an R color code. Color of the texture 
### class polygon boundary lines. 

 class.p.bg.col=NULL,
### Single logical OR vector of R colors (character strings). 
### If FALSE (the default), no color gradient 
### is used inside the texture class polygons. If TRUE, a color 
### gradient is drawn, with the color hue specified in 'class.p.bg.hue' 
### and with saturation and values that vary with texture. If 
### 'class.p.bg.col' is a vector of R colors of the same length 
### as the number of classes in the triangle, these colors 
### will be used as background color for each texture classe plygons.

 class.p.bg.hue=NULL,
### Single numerical. Only used if class.p.bg.col == TRUE (no default). 
### Color hue (between 0 and 1) used to create a color gradient 
### between the different texture class polygons.

# ARROWS (axis) and arrows LABELS parameters:

 arrows.show=NULL,
### Single logical. If TRUE (default), 3 arrows are drawn outside 
### the triangle, along each axis, that show the direction of 
### increasing values (arrow base) and of isovalue (arrow tip) 
### of the texture class. If FALSE no arrows are drawn.

 arrows.lty=NULL,
### Single numerical. Line type of the arrows drawn outside 
### the triangle, along each axis. Same possible types as par("lty").

# POINTS DATA parameters:

 points.type=NULL,
### Single text letter. Point type. Either "p" (points only), "l" 
### (lines only) or "b" (both points and lines), as for plot() or 
### points(). Point refer here to soil texture values plotted on 
### the triangle.

 pch=NULL,
### Single numerical or vector of numericals, or single text string 
### or vector of text string. Point shape number(s) or point character(s) 
### to be plotted. Point refer here to soil texture values plotted on 
### the triangle.

# Added 20090617 
# ADDITIONAL 4TH VARIABLE: z

 z.type=NULL,
### Single character string. Type of plot to be used for displaying 
### a 4th variable on the texture triangle (in addition to Clay, 
### Silt and Sand). Only used if 'z.name' is not NULL. Currently 
### only one value is supported, "bubble", for displaying a bubble 
### plot with bubble sizes and color saturation and values proportional 
### to the value of tri.data[,z.name]. 
### The value 'map' is deprecated and replaced by TT.iwd(), TT.image() or 
### TT.contour().

 z.col.hue=NULL,
### Single numerical. Hue of the bubble color ([0-1]) to be used 
### if 'z.name' is not NULL. A gradient of saturation and value is 
### automatically created for the bubbles (with this hue). 

 z.cex.range=NULL,
### Vector of 2 numericals. Minimum and maximum 'cex' of the bubbles 
### plotted on the triangle if 'z.name' is not NULL. 

 z.pch=NULL,
### Single numerical or vector of numericals. Point symbol number(s) 
### to be used for the bubbles if 'z.name' is not NULL. 

# TERNARY VARIABLES control tolerance

 text.tol=NULL,
### Single numerical. Tolerance on the sum of the 3 particle size classes. 
### The real sum of the 3 particle size classes in 
### 'tri.data' should be >= text.sum * (1-text.tol) OR 
### <= text.sum * (1+text.tol). See 'text.sum' for more details, as 
### well as 'tri.sum.tst' (to prevent texture sum tests).

# DATA TESTS:

 tri.sum.tst=NULL,
### Single logical. If TRUE (the default), the sum of the 3 texture 
### classes of each texture value in 'tri.data' will be checked 
### in regard to 'text.sum' and 'text.tol'. If FALSE, no test 
### is done.

 tri.pos.tst=NULL,
### Single logical. If TRUE (the default), the position of texture 
### values in 'tri.data' are tested to check that they are not 
### OUTSIDE the texture triangle (i.e. that some texture values may 
### be negative).

# ADDITIONAL parameters:

 b.lim=NULL,
### Vector of 2 numerical values. This is an equivalent to plot() 
### xlim argument. Minimum and maximum x / bottom value of the 
### texture triangle area, in FRACTION OF THE MAXIMAL EXTENSION. 
### Default is c(0,1). The real span is then b.lim * text.sum. 
### This is a minimal 'zoom' implementation (results are not 
### always perfect). 'b.lim' and 'l.lim' should be equal for 
### better rendering.

 l.lim=NULL,
### Vector of 2 numerical values. This is an equivalent to plot() 
### ylim argument. Minimum and maximum y / left value of the 
### texture triangle area, in FRACTION OF THE MAXIMAL EXTENSION. 
### Default is c(0,1). The real span is then l.lim * text.sum. 
### This is a minimal 'zoom' implementation (results are not 
### always perfect). 'b.lim' and 'l.lim' should be equal for 
### better rendering.

 lang=NULL,
### Single text string. Determines the language used for the plot 
### main title and axis labels. Possible values are 'en' (English, 
### the default), "fr" (French), "it" (Italian), "es" (Spanish), 
### "de" (German), "nl" (Dutch), "se" (Swedish) and "fl" (Flemish).

# Graph margins:

 new.mar=NULL,
### Vector of 4 numericals. Margin sizes of the plot. Default is 
### the same as par("mar"). See par("mar") for more details. Use 
### this at your own risks!

 new.centroid=TRUE  
### Single logical. If TRUE (default) the new method (Paul Bourke) 
### is used to calculate the centroid. If FALSE the centroid is 
### taken as the mean x and y coordinates of the vertices.

){  #
    if( is.null( class.sys ) ){ class.sys <- TT.get("class.sys") } 
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Create the "base" empty plot
    if( is.null( geo ) ) 
    {   #
        geo <- TT.baseplot( 
            geo             = geo, 
            class.sys       = class.sys, 
            # GEO PARAMETERS
            blr.clock       = blr.clock, 
            tlr.an          = tlr.an, 
            blr.tx          = blr.tx, 
            text.sum        = text.sum, 
            base.css.ps.lim = base.css.ps.lim, 
            # END GEO
            tri.sum.tst     = tri.sum.tst, 
            tri.pos.tst     = tri.pos.tst, 
            text.tol        = text.tol, 
            unit.ps         = unit.ps, 
            unit.tx         = unit.tx, 
            b.lim           = b.lim, 
            l.lim           = l.lim, 
            main            = main, 
            new.mar         = new.mar, 
            bg              = bg, 
            fg              = fg, 
            col             = col, 
            cex.main        = cex.main, 
            lang            = lang  
        )   #
    }   #
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    # --> Maybe useless in the future
    TT.geo.set( 
        geo     = geo  
        #p.env  = environment()  
    )   # 
    # 
    if( any(is.null(tri.css.ps.lim)) & (class.sys != "none") )
    {   #
        tri.css.ps.lim  <- TT.get(class.sys)[["tri.css.ps.lim"]] 
    }   #
    # 
    if( any(is.null(dat.css.ps.lim)) )
    {   #
        dat.css.ps.lim  <- base.css.ps.lim 
    }   #
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatically sets remaining parameters 
    TT.auto.set() 
    # 
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Test the data table provided:
    if( !is.null(tri.data) ) 
    {   #
        TT.data.test( 
            tri.data    = tri.data, 
            css.names   = css.names, 
            text.sum    = text.sum, 
            text.tol    = text.tol, 
            tri.sum.tst = tri.sum.tst, 
            tri.pos.tst = tri.pos.tst  
            #set.par     = FALSE  # Useless?? # Added 2009/06/27 
        )   #
    }   #   
    # 
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Setup might.cols option (dark or light background) 
    night.cols  <-  TT.col2hsv(bg)[,"v"] < 0.5 
    #
    if( z.type == "map" )
    {   #
        message( 
            'z.type == "map" is a deprecated parameter. Consider using TT.iwd(), TT.image() or TT.contour()'  
        )   #
    }   #
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # (optional) create a heatmap with z data:
    # if( !is.null(z.name) & z.type == "map" ) 
    # {   #
    #     blr.tx <- TT.blr.tx.check( 
    #         blr.tx      = blr.tx, 
    #         css.names   = css.names
    #     )   # 
    #     #
    #     base.frame2             <- TT.get("base.frame") * text.sum
    #     colnames(base.frame2)   <- blr.tx 
    #     #
    #     TT.iwd.image( 
    #         tri.data        = tri.data, 
    #         z.name          = z.name, 
    #         geo             = geo, 
    #         css.names       = css.names, 
    #         tri.pol.data    = base.frame2, 
    #         grid.nx         = 20, 
    #         grid.ny         = 20, 
    #         pow             = 2, 
    #         min.dist        = 25, 
    #         add             = TRUE, 
    #         z.col.hue       = z.col.hue, 
    #         image.n.col     = 20, 
    #         text.tol        = text.tol, 
    #         tri.sum.tst     = tri.sum.tst, 
    #         tri.pos.tst     = tri.pos.tst, 
    #         bg              = bg  
    #     )   #
    # }   #
    # 
    # TRIANGLE-FRAME BACKGROUND PLOT
    if( !add )
    {   #
        TT.edges( 
            geo             = geo, 
            text.tol        = text.tol, 
            col.axis        = col.axis, 
            plot.axis       = FALSE, 
            frame.bg.col    = frame.bg.col, 
            lwd.axis        = lwd.axis, 
            tri.sum.tst     = tri.sum.tst, 
            tri.pos.tst     = tri.pos.tst, 
            bg              = bg, 
            text.sum        = text.sum, 
            blr.clock       = blr.clock  
        )   #
    }   #
    #
    # +----------------------------------------------------+
    # Plot the classes polygon fill/background, but not border
    # before/below the grid lines
    #
    # - Define if class.p.bg.col is TRUE (or not a logical)
    class.p.bg.col.test <- is.logical( class.p.bg.col ) 
    if( class.p.bg.col.test ) 
    {   #
        class.p.bg.col.test <- class.p.bg.col 
    }else{ 
        class.p.bg.col.test <- TRUE 
    }   #
    #
    if( (class.sys != "none") & class.p.bg.col.test ) 
    {   #
        TT.classes(
            geo             = geo, 
            class.sys       = class.sys, 
            tri.css.ps.lim  = tri.css.ps.lim, 
            css.transf      = css.transf, 
            text.transf.fun = text.transf.fun, 
            trsf.add.opt1   = trsf.add.opt1, 
            trsf.add.opt2   = trsf.add.opt2, 
            text.tol        = text.tol, 
            tri.sum.tst     = tri.sum.tst, 
            tri.pos.tst     = tri.pos.tst, 
            class.lab.col   = class.lab.col, 
            class.p.bg.col  = class.p.bg.col, 
            class.p.bg.hue  = class.p.bg.hue, 
            class.line.col  = NA, 
            class.lab.show  = class.lab.show, 
            cex.lab         = cex.lab, 
            font.lab        = font.lab, 
            family.op       = family.op, 
            lwd.axis        = lwd.axis, 
            col.axis        = col.axis, 
            bg              = bg, 
            text.sum        = text.sum, 
            base.css.ps.lim = base.css.ps.lim, 
            blr.tx          = blr.tx, 
            blr.clock       = blr.clock, 
            new.centroid    = new.centroid  
        )   #
    }   #
    #
    if( grid.show )
    {   #
        TT.grid( 
            geo             = geo, 
            at              = at, 
            text.tol        = text.tol, 
            grid.col        = grid.col, 
            grid.lty        = grid.lty, 
            lwd.axis        = lwd.axis, 
            tri.sum.tst     = tri.sum.tst, 
            tri.pos.tst     = tri.pos.tst, 
            class.p.bg.col  = class.p.bg.col,   # added 2009/05/18 
            class.p.bg.hue  = class.p.bg.hue,   # added 2009/05/18 
            frame.bg.col    = frame.bg.col,     # added 2009/05/19 
            bg              = bg,               # added 2009/05/22 
            col.axis        = col.axis,         # added 2009/05/22 
            text.sum        = text.sum, 
            blr.clock       = blr.clock  
        )   #
    }   #
    #
    # +----------------------------------------------------+
    # Plot the classes polygon lines, but no background
    # after/over the grid lines
    if( class.sys != "none" )
    {   #
        TT.classes(
            geo             = geo, 
            class.sys       = class.sys, 
            tri.css.ps.lim  = tri.css.ps.lim, 
            css.transf      = css.transf, 
            text.transf.fun = text.transf.fun, 
            trsf.add.opt1   = trsf.add.opt1, 
            trsf.add.opt2   = trsf.add.opt2, 
            text.tol        = text.tol, 
            tri.sum.tst     = tri.sum.tst, 
            tri.pos.tst     = tri.pos.tst, 
            class.lab.col   = class.lab.col, 
            class.p.bg.col  = FALSE, 
            class.p.bg.hue  = class.p.bg.hue, 
            class.line.col  = class.line.col, 
            class.lab.show  = class.lab.show, 
            cex.lab         = cex.lab, 
            font.lab        = font.lab, 
            family.op       = family.op, 
            lwd.axis        = lwd.axis, 
            col.axis        = col.axis, 
            bg              = bg, 
            text.sum        = text.sum, 
            base.css.ps.lim = base.css.ps.lim, 
            blr.tx          = blr.tx, 
            blr.clock       = blr.clock, 
            new.centroid    = new.centroid  
        )   #
    }   #
    #
    # +----------------------------------------------------+
    #
    TT.ticks( 
        geo         = geo, 
        at          = at, 
        text.tol    = text.tol, 
        tk.s        = TT.get("ticks.shift"), 
        tri.sum.tst = tri.sum.tst, 
        lwd.axis    = lwd.axis,         # Added   2009/05/15 
        col.axis    = col.axis,         # Changed 2009/05/18 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    TT.ticks.lab( 
        geo         = geo, 
        at          = at, 
        text.tol    = text.tol, 
        tk.ls       = TT.get("ticks.lab.shift"), 
        tri.sum.tst = tri.sum.tst, 
        col.axis    = col.axis, 
        font.axis   = font.axis, 
        cex.axis    = cex.axis, 
        family.op   = family.op, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock, 
        tlr.an      = tlr.an  
    )   #
    #
    if( arrows.show )
    {   #
        TT.axis.arrows( 
            geo             = geo, 
            css.lab         = css.lab, 
            a.l             = TT.get("arrows.lims"), 
            a.h.s           = TT.get("arrows.head.shift"), 
            a.t.s           = TT.get("arrows.text.shift"), 
            a.b.s           = TT.get("arrows.base.shift"), 
            text.tol        = text.tol,     # useless
            base.css.ps.lim = base.css.ps.lim, 
            tlr.an          = tlr.an, 
            lwd.lab         = lwd.lab, 
            arrows.lty      = arrows.lty, 
            col.lab         = col.lab, 
            font.lab        = font.lab, 
            cex.lab         = cex.lab, 
            family.op       = family.op, 
            unit.ps         = unit.ps, 
            unit.tx         = unit.tx, 
            lang            = lang, 
            text.sum        = text.sum, 
            blr.clock       = blr.clock  
        )   #
    }   #
    #
    # TRIANGLE-FRAME FOREGROUND-BORDER PLOT
    TT.edges( 
        geo             = geo, 
        text.tol        = text.tol, 
        col.axis        = col.axis, 
        plot.axis       = TRUE, 
        frame.bg.col    = NA, 
        lwd.axis        = lwd.axis, 
        tri.sum.tst     = tri.sum.tst, 
        tri.pos.tst     = tri.pos.tst, 
        bg              = bg, 
        text.sum        = text.sum, 
        blr.clock       = blr.clock  
    )   #
    #
    if( !is.null(tri.data) )
    {   #
        points.out <- TT.points( 
            tri.data        = tri.data, 
            geo             = geo, 
            base.css.ps.lim = base.css.ps.lim, 
            dat.css.ps.lim  = dat.css.ps.lim, 
            css.transf      = css.transf, 
            text.transf.fun = text.transf.fun, 
            trsf.add.opt1   = trsf.add.opt1, 
            trsf.add.opt2   = trsf.add.opt2, 
            css.names       = css.names, 
            z.name          = z.name, 
            text.tol        = text.tol, 
            pch             = pch, 
            fg              = fg, 
            col             = col, 
            bg              = bg, 
            cex             = cex, 
            lwd             = lwd, 
            points.type     = points.type, 
            tri.sum.tst     = tri.sum.tst, 
            tri.pos.tst     = tri.pos.tst, 
            z.type          = z.type, 
            z.col.hue       = z.col.hue, 
            z.cex.range     = z.cex.range,  
            z.pch           = z.pch, 
            text.sum        = text.sum, 
            blr.clock       = blr.clock, 
            blr.tx          = blr.tx  
        )   #
    }else{ 
        points.out <- NULL
    }   #
    #
    return( invisible( geo ) ) 
}   #






TT.points.in.classes <- function(# Classify a table of soil texture data according to a soil texture triangle.
### The function calculate in which classe(s) of a texture triangle 
### (classification system defined by 'class.sys') lies each soil 
### sample (with texture data) in the table 'tri.data'. As a sample 
### may lie inside a texture class, but also at the edge of 2 or 
### more texture classes, the function does not only output 
### one single texture class per sample. If 'PiC.type' is 'n' or 
### 'l', it rather output a table where each column is a texture 
### class and each row a texture sample, and yes / no information 
### about the belonging of the sample to each texture class. 
### Alternatively, If 'PiC.type' is 't'it will output a text 
### string (per sample) containing all the texture classes 
### to which that point belong.
### The texture data in 'tri.data' can be transformed into 
### another particle size system prior to their classification 
### if needed. See the options  base.css.ps.lim, tri.css.ps.lim, 
### dat.css.ps.lim, css.transf and text.transf.fun.
### ON DEFAULT VALUES OF TT.points.in.classes() ARGUMENTS? As 
### TT.points.in.classes() shares 
### its arguments with many other functions, their default value 
### is not defined in TT.points.in.classes() source code, but rather in 
### a dedicated list object called 'TT.par' and stored in the 
### environment TT.env. The function TT.get() is used to retrieve 
### the default value of the arguments defined in TT.par (see 
### ?TT.get). For instance, to know the default value of 'class.sys', 
### you can type TT.get("class.sys"). To set a different default 
### value for a given argument in R, use TT.set() (see ?TT.set). 
### For instance to change the default value of 'class.sys', type 
### TT.set( "class.sys" = "USDA.TT" ).

 tri.data,
### Data frame. Data frame containing the CLAY, SILT and SAND 
### 'coordinates' of the texture data points to be classified The data 
### frame can contain more column than needed (ignored). The data 
### frame must have column named CLAY, SILT and SAND (uppercase, 
### the order has no importance) or named after the 'css.names' 
### argument (alternative names). The sum of CLAY, SILT and SAND 
### must be equal to 'text.sum' 
### ('text.tol' determines the error tolerance).

 class.sys=NULL,
### Single text string. Text code of the texture classification 
### system to be used for the classification of 'tri.data'. 
### Possible values are "none" (no classification plotted), "USDA.TT" 
### (USDA texture triangle), "HYPRES.TT" (texture triangle of the 
### European Soil Map), "FR.AISNE.TT" (French texture triangle of 
### the Aisne region soil survey), "FR.GEPPA.TT" (French GEPPA 
### texture triangle), "DE.BK94.TT" (German texture triangle), 
### "UK.SSEW.TT" (Soil Survey of England and Wales), "AU.TT" 
### (Australian texture triangle), "BE.TT" (Belgium texture triangle), 
### "CA.EN.TT" (Canadian texture triangle, with English class abbreviations) and 
### "CA.FR.TT" (Canadian texture triangle, with French class abbreviations)
### (see the package vignette for a complete list).

 PiC.type=NULL,
### Single character string. If equal to 'n', then a table of 0, 
### 1, 2 or 3 is outputed (0 if the sample does not belong to a class, 
### 1 if it does, 2 if it lies on an edge and 3 if it lies on a 
### vertex). Notice that the accuracy of the classification is 
### not garanteed for samples lying very close to an edge, or right 
### on it. See <http://www.mail-archive.com/r-help@r-project.org/msg96180.html>

 css.names=NULL,
### Vector of 3 character strings. Name of the columns in 'tri.data' 
### that contains the CLAY SILT and SAND values, respectively. 
### If NULL, default c("CLAY","SILT","SAND") value is assumed. Not 
### to be confused with 'css.lab' that defines the labels of the 
### CLAY SILT and SAND axes in the plot.

 text.sum=NULL,
### Single numerical. Sum of the 3 particle size classes for each texture 
### value (fixed). The real sum of the 3 particle size classes in 'tri.data' 
### should be >= text.sum * (1-text.tol) OR  <= text.sum * (1+text.tol), 
### where 'text.tol' is an argument that can be changed. If some 
### of the texture values don't match this requirement, an error 
### occur (function fails) and TT.points.in.classes returns a of bad values with 
### their actual particle size classes sum. You can 'normalise' you data 
### table () prior to the use of TT.points.in.classes, by using the function 
### TT.normalise.sum(), so all values match the 'text.sum' criteria. 
### See also 'tri.sum.tst' that can be set to FALSE to avoid 
### sum of particle size classes tests.

 base.css.ps.lim=NULL,
### Vector of 4 numericals. Particle size boundaries (upper and lower) 
### of the 3 particle size classes (CLAY, SILT and SAND, starting from 
### the lower size of CLAY particles, 0, to the upper size of the 
### SAND particles, 2000), in micrometers, FOR THE BASE SYSTEM. These 
### particles size class limits are the references and all other 
### texture values with different limits will be converted into 
### that reference if (and only if) css.transf == TRUE (not default). 
### If NULL, 'base.css.ps.lim' will be set to the default value of the 
### texture classification system chosen ('class.sys'). The 
### transformation function is set by 'text.transf.fun' and is 
### a log-linear interpolation by default.

 tri.css.ps.lim=NULL,
### Vector of 4 numericals. Particle size boundaries (upper and lower) 
### of the 3 particle size classes (CLAY, SILT and SAND, starting from 
### the lower size of CLAY particles, 0, to the upper size of the 
### SAND particles, 2000), in micrometers, FOR THE TEXTURE TRIANGLE. 
### If not NULL, different from 'base.css.ps.lim', and 
### css.transf == TRUE (not default), then the CLAY SILT and SAND 
### coordinates of the texture triangle will be converted into 
### the 'base.css.ps.lim' reference. If NULL, 'tri.css.ps.lim' will 
### be set to the default value of the texture classification system 
### chosen ('class.sys'). The transformation function is set by 
### 'text.transf.fun' and is a log-linear interpolation by default.

 dat.css.ps.lim=NULL,
### Vector of 4 numericals. Particle size boundaries (upper and lower) 
### of the 3 particle size classes (CLAY, SILT and SAND, starting from 
### the lower size of CLAY particles, 0, to the upper size of the 
### SAND particles, 2000), in micrometers, FOR THE TEXTURE DATA TABLE
### ('tri.data'). If not NULL, different from 'base.css.ps.lim', and 
### css.transf == TRUE (not default), then the CLAY SILT and SAND 
### coordinates of the texture data in tri.data will be converted into 
### the 'base.css.ps.lim' reference. If NULL, 'tri.css.ps.lim' will 
### be set to the default value of the texture classification system 
### chosen ('class.sys'). The transformation function is set by 
### 'text.transf.fun' and is a log-linear interpolation by default.

 css.transf=NULL,
### Single logical. Set to TRUE to transform the texture coordinates 
### of the texture triangle ('class.sys') or the texture data 
### ('tri.data') into the base particle size class limits. 
### See 'base.css.ps.lim' for the base plot particle size class limits, 
### 'tri.css.ps.lim' for the triangle particle size class limits 
### and 'dat.css.ps.lim' for the data table particle size class limits. 
### The transformation function is set by 'text.transf.fun' and 
### is a log-linear interpolation by default. The default value is 
### FALSE, so no transformation is made.

 text.transf.fun=NULL,
### R function with the same argument names and same output as 
### the function TT.text.transf(). 'text.transf.fun' is the function 
### that transform the texture values from one system of particle 
### class size limits to another. Only used if css.transf == TRUE. 
### Default value is text.transf.fun=TT.text.transf. See also 
### 'base.css.ps.lim', 'tri.css.ps.lim' and 'dat.css.ps.lim'.

 trsf.add.opt1=NULL,
### Non pre-defined format. If the user specifies its own texture 
### transformation function in 'text.transf.fun' (not TT.text.transf()), 
### then he can use 'trsf.add.opt1' and 'trsf.add.opt1' as 
### new, additional, argument for his function. So the format of 
### 'trsf.add.opt1' depends on the function defined by the user 
### in 'text.transf.fun'.

 trsf.add.opt2=NULL,
### Non pre-defined format. If the user specifies its own texture 
### transformation function in 'text.transf.fun' (not TT.text.transf()), 
### then he can use 'trsf.add.opt1' and 'trsf.add.opt1' as 
### new, additional, argument for his function. So the format of 
### 'trsf.add.opt1' depends on the function defined by the user 
### in 'text.transf.fun'.

 text.tol=NULL,
### Single numerical. Tolerance on the sum of the 3 particle size classes. 
### The real sum of the 3 particle size classes in 
### 'tri.data' should be >= text.sum * (1-text.tol) OR 
### <= text.sum * (1+text.tol). See 'text.sum' for more details, as 
### well as 'tri.sum.tst' (to prevent texture sum tests).

 tri.sum.tst=NULL,
### Single logical. If TRUE (the default), the sum of the 3 texture 
### classes of each texture value in 'tri.data' will be checked 
### in regard to 'text.sum' and 'text.tol'. If FALSE, no test 
### is done.

 tri.pos.tst=NULL,
### Single logical. If TRUE (the default), the position of texture 
### values in 'tri.data' are tested to check that they are not 
### OUTSIDE the texture triangle (i.e. that some texture values may 
### be negative).

 collapse=NULL,
### Single character string. If PiC.type = "t" and a sample lie 
### on the edge of 2 texture classes, then both will be outputed 
### in a single character string, separated by 'collapse'. Example of 
### output: [1] "C" "VF, F" "C" "C" "M"

 texture2xy=FALSE,
### Single logical. Set to FALSE to avoid any transformation of the 
### texture data (trigonometric) prior to testure data classification. 
### Setting to FALSE avoid some numerical accuracy problems when 
### a point is on the border of a texture class.

 blr.tx=NULL,
### Vector of 3 character strings. The 1st, 2nd and 3rd values must 
### be either CLAY, SILT or SAND, and determines the particle size classes 
### associated with the BOTTOM, LEFT and RIGHT axis, respectively. 
### CLAY, SILT and SAND order in the vector is free, but they should 
### all be used one time. The CLAY, SILT and SAND names must appear 
### whatever the corresponding columns names in 'tri.data' (eventually 
### set by 'css.names') and whatever the labels of the axis in the 
### plot (eventually set by 'css.lab') 

 blr.clock=NULL
### Vector of logicals, eventually with NA values. Direction of 
### increasing texture values on the BOTTOM, LEFT and RIGHT axis, 
### respectively. A value of TRUE means that the axis direction is 
### clockwise. A value of FALSE means that the axis direction is 
### counterclockwise. A value of NA means that the axis direction 
### is centripetal. Possible combinations are c(T,T,T); c(F,F,F); 
### c(F,T,NA) and c(T,NA,F), for fully clockwise, fully counterclockwise, 
### right centripetal and left centripetal orientations, respectively.

){  #
    if( is.null( class.sys ) ){ class.sys <- TT.get("class.sys") } 
    
    # require( "sp" ) 
    
    TT.data <- TT.get(class.sys) 
    #
    if( is.null(collapse) ){ collapse = ", " } 
    #
    geo <- TT.geo.get( 
        class.sys       = class.sys, 
        blr.clock       = NULL, 
        tlr.an          = NULL, 
        blr.tx          = NULL, 
        text.sum        = text.sum, 
        base.css.ps.lim = base.css.ps.lim  
    )   #
    #
    if( any( is.null(tri.css.ps.lim) ) )
    {   #
        tri.css.ps.lim  <- TT.data[["tri.css.ps.lim"]] 
    }   #
    #
    if( any( is.null(base.css.ps.lim) ) )
    {   #
        base.css.ps.lim  <- tri.css.ps.lim 
    }   #
    #
    if( any( is.null(dat.css.ps.lim) ) )
    {   #
        dat.css.ps.lim  <- base.css.ps.lim  # Equal to tri.css.ps.lim, 
        #                                   # except if base.css.ps.lim
        #                                   # is specified
    }   #
    #
    # Set geographical parameters:
    TT.geo.set( 
        geo     = geo  
        #p.env  = environment()  
    )   #
    #
    TT.auto.set( set.par = FALSE ) 
    #
    TT.data.test( 
        tri.data    = tri.data, 
        css.names   = css.names, 
        text.sum    = text.sum, 
        text.tol    = text.tol, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst  
    )   #
    # 
    if( css.transf )
    {   #
        text.transf.fun <- get( text.transf.fun ) 
        #
        tri.data <- text.transf.fun( 
            tri.data        = tri.data,  
            base.css.ps.lim = base.css.ps.lim,  
            dat.css.ps.lim  = dat.css.ps.lim,  
            css.names       = css.names, 
            blr.tx          = blr.tx,  
            text.sum        = text.sum,  
            text.tol        = text.tol,  
            tri.sum.tst     = tri.sum.tst,  
            tri.pos.tst     = tri.pos.tst,  
            trsf.add.opt1   = trsf.add.opt1, 
            trsf.add.opt2   = trsf.add.opt2  
        )   #
        #
        TT.data$"tt.points" <- text.transf.fun( 
            tri.data        = TT.data$"tt.points",  
            base.css.ps.lim = base.css.ps.lim,  
            dat.css.ps.lim  = tri.css.ps.lim,  
            css.names       = NULL, # TT.data$"tt.points" names are not necessarily identical to those in css.name
            blr.tx          = blr.tx,  
            text.sum        = 1, 
            text.tol        = text.tol,  
            tri.sum.tst     = tri.sum.tst,  
            tri.pos.tst     = tri.pos.tst,  
            trsf.add.opt1   = trsf.add.opt1, 
            trsf.add.opt2   = trsf.add.opt2  
        )   #
    }   #
    #
    # Operate trigonometric transformation of the values into 
    # the x-y pland, or just used 2 texture values?
    if( texture2xy )
    {   #
        classes.points.xy <- TT.css2xy( 
            tri.data    = TT.data$"tt.points" * text.sum, 
            geo         = geo, 
            # 
            css.names   = css.names, 
            text.tol    = text.tol, 
            #
            tri.sum.tst = tri.sum.tst, 
            tri.pos.tst = tri.pos.tst, 
            set.par     = FALSE, 
            text.sum    = text.sum, 
            blr.clock   = blr.clock  
        )   #
        #
        data.points.xy <- TT.css2xy( 
            tri.data    = tri.data, 
            geo         = geo, 
            # 
            css.names   = css.names, 
            text.tol    = text.tol, 
            #
            tri.sum.tst = tri.sum.tst, 
            tri.pos.tst = tri.pos.tst, 
            set.par     = FALSE, 
            text.sum    = text.sum, 
            blr.clock   = blr.clock  
        )   #
    }else{ 
        classes.points.xy <- data.frame( 
            "xpos" = TT.data$"tt.points"[,"SILT"] * text.sum, 
            "ypos" = TT.data$"tt.points"[,"CLAY"] * text.sum  
        )   #
        #
        data.points.xy <- data.frame( 
            "xpos" = tri.data[,css.names[2]], 
            "ypos" = tri.data[,css.names[1]]  
        )   #
    }   #
    
    # require( "sp" ) 
    
    # Vectorisable and custom wrapper for point.in.polygon():
    points.in.class <- function( 
        X,                  # X = class name
        polygons.list       = TT.data$"tt.polygons", 
        classes.points.xy   = classes.points.xy, 
        data.points.xy      = data.points.xy  
    ){  #
        sel.vec <- (polygons.list[[ X ]])$"points"
        #
        xpol    <- classes.points.xy[ sel.vec , "xpos" ]
        ypol    <- classes.points.xy[ sel.vec , "ypos" ]
        #
        PiP.res <- point.in.polygon( 
            point.x = data.points.xy$"xpos",  
            point.y = data.points.xy$"ypos",  
            pol.x   = xpol,  
            pol.y   = ypol  
            # mode.checked = TRUE 
        )   # 
        #
        return( PiP.res ) 
    }   #
    #
    classes.names <- names( TT.data$"tt.polygons" )
    #
    PiP.res2 <- do.call( 
        "cbind", 
        lapply( 
            X                   = classes.names, 
            FUN                 = points.in.class, 
            polygons.list       = TT.data$"tt.polygons",  
            classes.points.xy   = classes.points.xy, 
            data.points.xy      = data.points.xy  
        )   #
    )   #
    #
    if( PiC.type %in% c("l","t") )
    {   #
        PiP.res2 <- t( apply( 
            X       = PiP.res2, 
            MARGIN  = 1, 
            FUN     = function(X){ 
                as.logical( X ) 
            }   #
        ) )   #
    }   #
    #
    colnames(PiP.res2) <- classes.names
    #
    if( PiC.type == "t" )
    {   #
        PiP.res2 <- unlist( apply( 
            X       = PiP.res2, 
            MARGIN  = 1, 
            FUN     = function(X){ 
                paste( names(X)[X], collapse = collapse ) 
            }   #
        ) )   #
    }   #
    #
    return( PiP.res2 ) 
}   #

#     my.text2 <- my.text

#     my.text2 <- do.call( 
#         what    = "rbind", 
#         args    = lapply( 
#             X   = 1:8500, 
#             FUN = function(X){ 
#                 my.text
#             }   #
#         )   #
#     )   #

#     dim(my.text2)

#     system.time( 
#         TT.points.in.classes( 
#             tri.data    = my.text2, 
#             class.sys   = "USDA.TT"  
#         )   #
#     )   #
#     # 14000 texture classified per second






TT.xy2css <- function(# Internal. Convert point-data duplets (2 variables, x-y coordinaes) in Clay silta and sand coordinates. 
### Internal. Convert point-data duplets (2 variables, x-y 
### coordinaes) in Clay silta and sand coordinates. 

 xy.data,
### a data.frame with xpos and ypos columns

 geo, 
 css.names       = NULL, 
 
 text.tol        = NULL, 
 
 tri.sum.tst     = NULL, 
 tri.pos.tst     = NULL, 
 set.par         = FALSE, 
 blr.clock       = NULL, 
 text.sum        = NULL  

){  #
    if( !(class(xy.data) %in% c("data.frame","matrix")) )
    {   #
        stop( 
            paste( 
                sep = "", 
                "xy.data must be a data.frame or a matrix (now ", 
                class(xy.data), ")"   
            )   #
        )   #
    }   #
    #
    if( any(!(colnames(xy.data) %in% c("xpos","ypos"))) )
    {   #
        stop( 
            paste( 
                sep = "", 
                "xy.data must be a have xpos and ypos in column names (now ", 
                paste(colnames(xy.data),collapse=", "), ")"   
            )   #
        )   #
    }   #
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    # Set rest of variables:
    TT.auto.set(set.par=set.par) 
    #
    blr.tx <- TT.blr.tx.check( 
        blr.tx      = blr.tx, 
        css.names   = css.names
    )   # 
    #
    # Check for tlr.an: data type
    if( any( is.null(tlr.an) | is.na(tlr.an) ) | !is.numeric(tlr.an) | (length(tlr.an) != 3) )
    {   #
        stop( paste( 
            sep = "", 
            "tlr.an (=", 
            paste(tlr.an,collapse=";"), 
            ") must be a numeric, non-null, non-na vector of length 3" 
        )   )  
    }   #
    #
    # Check for tlr.an: angle sum must be 180 degrees
    if( sum(tlr.an) != 180 ) 
    {   #
        stop( paste( 
            sep = "", 
            "sum(tlr.an) (=", 
            paste(tlr.an,collapse=";"), 
            ") must be 180 (degrees)"  
        )   )  
    }   #
    # 
    # Test (anti)clock settings
    ok.clock <- list( 
        #       #    Bottom Left    Right  
        "TTT"   = c( TRUE,  TRUE,   TRUE    ), 
        "TXF"   = c( TRUE,  NA,     FALSE   ), 
        "FTX"   = c( FALSE, TRUE,   NA      ), 
        "FFF"   = c( FALSE, FALSE,  FALSE   )  
    )   #
    #
    ok.clock <- unlist( lapply( 
            X           = ok.clock, 
            FUN         = function(X,blr.clock){ 
                identical(blr.clock,X) 
            },  #
            blr.clock   = blr.clock
    )   )   #
    #
    if( !any(ok.clock) )
    {   #
        stop( paste( 
                sep = "", 
                "blr.clock (=", 
                paste(as.character(blr.clock),collapse=";"), 
                ") MUST be one of: ", 
                paste(names(ok.clock),collapse=";"), 
                "; [with X == NA]. consider revising" 
        )   )   #
    }   #
    #
    # Angle transformation: degree to radian
    tlr.an <- TT.deg2rad(tlr.an)
    #
    # # "reverse" the bottom and right orientation to fit x and y orientation
    rev.dt <- function( 
        i, 
        blr.c   = blr.clock, 
        tri.d   = tri.data, 
        blr.t   = blr.tx, 
        text.s  = text.sum  
    ){  #
        val <- tri.d[  , blr.t[i] ]
        if( !is.na(blr.c[i]) )  # Do not reverse NA sides
        {   #
            if( (blr.c[i] & (i != 2)) | (!blr.c[i] & (i == 2)) ) 
            {   #
                val <- ( text.s - val )
            }   #
        }   #
        val 
    }   #
    #
    # The x,y coordnates calculation is 1st separated depending on blr.clock[2]
    if( !is.na(blr.clock[2]) ){ cond2 <- blr.clock[2] }else{ cond2 <- FALSE } 
    #
    tri.data <- matrix( 
        data    = NA, 
        ncol    = 3, 
        nrow    = dim(xy.data)[1]  
    )   #
    tri.data <- as.data.frame( tri.data ) 
    colnames( tri.data ) <- css.names 
    #
    if( cond2 )         # if cond2 == TRUE: Either the TTT or the FTX case:
    {   #               #   * ypos depends on the Left (2) angle (clock)
        # ypos    <- tri.data[  , blr.tx[2] ] * sin(tlr.an[2])
        tri.data[  , blr.tx[2] ] <- xy.data[,"ypos"] / sin(tlr.an[2]) 
        #               #
        if( blr.clock[1] == TRUE )  
        {   #           #   * This is the TTT case:
            #xpos    <- tri.data[  , blr.tx[1] ] - ypos/tan(tlr.an[3])
            tri.data[  , blr.tx[1] ] <- xy.data[,"xpos"] + xy.data[,"ypos"]/tan(tlr.an[3])
        }else{          #   * This is the FTX case:
            #xpos    <- tri.data[  , blr.tx[1] ] + ypos/tan(tlr.an[2])
            tri.data[  , blr.tx[1] ] <- xy.data[,"xpos"] - xy.data[,"ypos"]/tan(tlr.an[2])
        }   #
        #
        for( j in 1:3 ){ tri.data[,blr.tx[j]] <- rev.dt("i"=j) } 
        #
        tri.data[  , blr.tx[3] ] <- text.sum - tri.data[  , blr.tx[2] ] - tri.data[  , blr.tx[1] ] 
    }else{              # if cond2 == FALSE: Either the FFF or the TXF case:
        #               #   * ypos depends on the Right (3) angle (clock)
        #ypos    <- tri.data[  , blr.tx[3] ] * sin(tlr.an[3])
        tri.data[  , blr.tx[3] ] <- xy.data[,"ypos"] / sin(tlr.an[3]) 
        #  
        if( blr.clock[1] == FALSE )  
        {   #           #   * This is the FFF case:
            #xpos    <- tri.data[  , blr.tx[1] ] + ypos/tan(tlr.an[2]) 
            tri.data[  , blr.tx[1] ] <- xy.data[,"xpos"] - xy.data[,"ypos"]/tan(tlr.an[2]) 
        }else{          #   * This is the TXF case:
            #xpos    <- tri.data[  , blr.tx[1] ] - ypos/tan(tlr.an[3]) 
            tri.data[  , blr.tx[1] ] <- xy.data[,"xpos"] + xy.data[,"ypos"]/tan(tlr.an[3]) 
        }   #
        #
        for( j in 1:3 ){ tri.data[,blr.tx[j]] <- rev.dt("i"=j) } 
        #
        tri.data[  , blr.tx[2] ] <- text.sum - tri.data[  , blr.tx[3] ] - tri.data[  , blr.tx[1] ] 
    }   #
    #
    # Test the data provided:
    TT.data.test( 
        tri.data    = tri.data, 
        css.names   = css.names, 
        text.sum    = text.sum, 
        text.tol    = text.tol, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst  
        #set.par     = set.par  
    )   #
    #
    return( tri.data )
}   #

#     rand.text <- TT.dataset(n=100,seed.val=1980042401) 

#     xy.data <- TT.css2xy( 
#         geo         = TT.geo.get( class.sys = "USDA.TT" ), 
#         tri.data    = rand.text
#     )   #

#     tri.data <- TT.xy2css( 
#         geo         = TT.geo.get( class.sys = "USDA.TT" ), 
#         xy.data    = xy.data
#     )   #

#     all( (rand.text - tri.data) < 1e-6 ) 
#     
#     
#     xy.data <- TT.css2xy( 
#         geo         = TT.geo.get( class.sys = "FR.GEPPA.TT" ), 
#         tri.data    = rand.text
#     )   #

#     tri.data <- TT.xy2css( 
#         geo         = TT.geo.get( class.sys = "FR.GEPPA.TT" ), 
#         xy.data    = xy.data
#     )   #

#     all( (rand.text - tri.data) < 1e-6 ) 






TT.locator <- function(# Interactive (mouse clic) retrieval the CLAY SILT SAND coordinate of points on a texture triangle.
### Interactive (mouse clic) retrieval the CLAY SILT SAND coordinate of points on a texture triangle. 

 geo, 
 css.names       = NULL, 
 
 text.tol        = NULL, 
 
 tri.sum.tst     = NULL, 
 tri.pos.tst     = FALSE, 
 set.par         = FALSE, 
 n               = 512, 
 type            = "n", 
 ... 
### Further argumets passed to locator() 

){  #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    # Set rest of variables:
    TT.auto.set(set.par=set.par) 
    # 
    xy.data <- locator(n=n,type=type,...) 
    # 
    xy.data <- data.frame( 
        "xpos"  = xy.data[["x"]], 
        "ypos"  = xy.data[["y"]]  
    )   #
    #
    tri.data <- TT.xy2css( 
        xy.data         = xy.data, 
        geo             = geo, 
        css.names       = css.names, 
        text.tol        = text.tol, 
        tri.sum.tst     = tri.sum.tst, 
        tri.pos.tst     = tri.pos.tst, 
        set.par         = set.par  
    )   #
    #
    return( tri.data ) 
}   #

#     geo <- TT.plot( class.sys = "USDA.TT" ) 
#     
#     TT.locator(geo=geo) 
#     
#     geo <- TT.plot( class.sys = "FR.GEPPA.TT" ) 
#     
#     TT.locator(geo=geo) 






TT.xy.grid <- function(# Internal. Create a grid in the x-y coordinate system. 
### Create a grid in the x-y coordinate system. Most of the function 
### is a reshaped extract from kde2d() from the MASS package, by 
### Venables & Ripley (+ modifications)
##keywords<< internal

    x,  
    y,  
    n   = 25  
){  #
    nx <- length(x)
    #
    if(length(y) != nx)
    {   #
        stop("data vectors must be the same length")
    }   #
    #
    if(any(!is.finite(x)) || any(!is.finite(y)))
    {   #
        stop("missing or infinite values in the data are not allowed") 
    }   #
    #
    lims <- c(range(x), range(y))
    #
    # gx <- seq.int( lims[1], lims[2], length.out = n ) 
    gx <- seq( from = lims[1], to = lims[2], length.out = n ) 
    # gy <- seq.int( lims[3], lims[4], length.out = n ) 
    gy <- seq( from = lims[3], to = lims[4], length.out = n ) 
    #
    return( 
        list( 
            "xpos"  = gx, 
            "ypos"  = gy, 
            "xypos" = expand.grid( 
                "xpos"  = gx, 
                "ypos"  = gy  
            )   #
        )   #
    )   #
}   #






TT.chemometrics.alr <- function(# Compute the additive log-ratio transformation of compositional data.
### Function that compute the additive 
### log-ratio transformation of compositional data (here texture 
### data). This a a copy-paste-and-rename of the alr function provided 
### by the package chemometrics: P. Filzmoser and K. Varmuza (2008). 
### chemometrics: Multivariate Statistical Analysis in Chemometrics. 
### R package version 0.4.
### The function has been modified so it returns NA when a value 
### is below or equal to zero (this happens when using a regular 
### grid of texture data, for practical reasons).
### The function has also been modified so it uses column name
### rather than column index.

    X, 
    divisorvar, 
    css.names  
){  #
    # additive logratio transformation
    # INPUT: X...data matrix, divisorvar...number of ratioing variable
    #
    X.alr <- X[,css.names] 
    #
    sel.exc <- apply( 
        X       = X[,css.names], 
        MARGIN  = 1, 
        FUN     = function(X){ 
            any( X <= 0 ) 
        }   #
    )   #
    #
    X.alr[ sel.exc,  ] <- NA 
    #
    X.alr[ !sel.exc, css.names ] <-   
        log10( 
            X[ !sel.exc, css.names ] /  
            X[ !sel.exc, css.names[divisorvar] ] 
        )   #
    #
    return( X.alr[, css.names[-divisorvar] ] ) 
}   #






TT.mahalanobis <- function(# Calculates the Mahalanobis distance between clay silt and sand.
### Function that calculated the Mahalanobis 
### distance between clay silt and sand, on a regular x-y grid 
### (back-transformed to Clay silt and sand for Mahalanobis 
### calculation). The underlying function is mahalanobis() by 
### R Development Core Team (2009)

    # Parameter for TT.css2xy() and TT.xy2css() 
    geo,  
    tri.data,  
    css.names   = NULL,  
    text.tol    = NULL,  
    text.sum    = NULL,  
    blr.clock   = NULL,  
    tri.sum.tst = NULL,  
    tri.pos.tst = NULL,  
    set.par     = FALSE,  
    # x-y Grid parameter:
    n           = 25,  
    # Mahalanobis parameter:
    center      = NULL,  #  default colMeans(tri.data[,css.names]) 
    cov.mat     = NULL,  #  default cov(tri.data[,css.names]) 
    inverted    = FALSE,  
    ..., 
    # Additional parameter related to Mahalanobis
    alr         = FALSE, #  If TRUE an additive log-ratio transformation
    #                    #  of the data is performed, and the Mahalanobis 
    #                    #  distance is computed on all classes but css.names[divisorvar] 
    divisorvar  = 2      #  The Mahalanobis distance will be computed 
    #                    #  on all the texture class but css.names[divisorvar]  
){  #
    if( length(divisorvar) != 1 |   
        !is.numeric(divisorvar) |   
        !(divisorvar %in% 1:3)   
    ){  #
        stop( 
            "divisorvar should be equal to 1, 2 or 3"
        )   #
    }   #
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    # Set rest of variables:
    TT.auto.set(set.par=set.par) 
    #
    xy.bound <- TT.get("base.frame") * text.sum 
    colnames( xy.bound ) <- css.names 
    #
    xy.bound <- TT.css2xy( 
        tri.data        = xy.bound,   
        geo             = geo, 
        text.tol        = text.tol, 
        css.names       = css.names, 
        tri.sum.tst     = tri.sum.tst, 
        tri.pos.tst     = tri.pos.tst, 
        set.par         = set.par, 
        text.sum        = text.sum, 
        blr.clock       = blr.clock  
    )   #
    #
    xy.grid <- TT.xy.grid( 
        x   = xy.bound[,"xpos"], 
        y   = xy.bound[,"ypos"], 
        n   = n  
    )   #
    #
    tri.grid <- TT.xy2css( 
        xy.data     = xy.grid[["xypos"]], 
        geo         = geo, 
        css.names   = css.names, 
        text.tol    = text.tol, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = FALSE, 
        set.par     = set.par  
    )   #
    #
    if( alr )
    {   #
        tri.grid <- TT.chemometrics.alr( 
            X           = tri.grid, 
            divisorvar  = divisorvar, 
            css.names   = css.names  
        )   #
        #
        tri.data <- TT.chemometrics.alr( 
            X           = tri.data, 
            divisorvar  = divisorvar, 
            css.names   = css.names  
        )   #
        #
        if( is.null(center)  ){ center  <- colMeans(tri.data) } 
        if( is.null(cov.mat) ){ cov.mat <- cov(tri.data) } 
    }else{ 
        tri.data <- tri.data[,css.names[-divisorvar]] 
        tri.grid <- tri.grid[,css.names[-divisorvar]] 
        #
        if( is.null(center)  ){ center  <- colMeans(tri.data) } 
        if( is.null(cov.mat) ){ cov.mat <- cov(tri.data) } 
    }   #
    #
    sel.exc <- apply( 
        X       = tri.grid, 
        MARGIN  = 1, 
        FUN     = function(X){ 
            any( is.na(X) ) 
        }   #
    )   #
    #
    maha <- rep(NA,dim(tri.grid)[1])
    #
    maha[!sel.exc] <- mahalanobis( 
        x           = tri.grid[!sel.exc,], 
        center      = center, 
        cov         = cov.mat, 
        inverted    = FALSE, 
        ...  
    )   #
    
    # require( "sp" ) 
    
    PiP <- as.logical(  point.in.polygon( 
        point.x = xy.grid[["xypos"]][,"xpos"],
        point.y = xy.grid[["xypos"]][,"ypos"], 
        pol.x   = xy.bound[,"xpos"], 
        pol.y   = xy.bound[,"ypos"]  
        #mode.checked    = FALSE  #  Deprecated
    )   )   #
    #
    maha[ !PiP ] <- NA 
    #
    return( 
        list( 
            "x" = xy.grid[["xpos"]], 
            "y" = xy.grid[["ypos"]], 
            "z" = matrix( 
                data    = maha, 
                nrow    = length(xy.grid[["ypos"]]), 
                ncol    = length(xy.grid[["xpos"]]),
                byrow   = FALSE 
            )   #
        )   #
    )   #
}   #

#     rand.text <- TT.dataset(n=100,seed.val=1980042401) 

#     geo <- TT.geo.get( class.sys = "USDA.TT" )

#     maha <- TT.mahalanobis( 
#         geo         = geo, 
#         tri.data    = rand.text
#     )   #






TT.kde2d <- function(# Calculated the 2D probabilty density on an x-y grid.
### Function that calculated the 2D probabilty 
### density on an x-y grid (and NOT on the clay silt sand 
### reference system). Wrapper around the kde2d function from the 
### MASS package.
    # Parameter for TT.css2xy() and TT.xy2css() 
    geo, 
    tri.data, 
    css.names   = NULL, 
    text.tol    = NULL, 
    text.sum    = NULL, 
    blr.clock   = NULL, 
    tri.sum.tst = NULL, 
    tri.pos.tst = NULL, 
    set.par     = FALSE, 
    # x-y Grid parameter:
    n           = 25, 
    #
    lims        = c("points","triangle")[2]  
){  #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    # Set rest of variables:
    TT.auto.set(set.par=set.par) 
    #
    # xy.bound <- TT.css2xy( 
    #     tri.data        = data.frame(
    #         "CLAY"  = c(1, 0, 0), 
    #         "SILT"  = c(0, 1, 0), 
    #         "SAND"  = c(0, 0, 1)  
    #     ) * text.sum,   
    #     geo             = geo  
    # )   #
    #
    xy.bound <- TT.get("base.frame") * text.sum 
    colnames(xy.bound) <- css.names 
    #
    xy.bound <- TT.css2xy( 
        tri.data    = xy.bound,   
        geo         = geo, 
        text.tol    = text.tol, 
        css.names   = css.names, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        set.par     = set.par, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    xy.data <- TT.css2xy( 
        tri.data    = tri.data, 
        geo         = geo, 
        text.tol    = text.tol, 
        css.names   = css.names, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        set.par     = set.par, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    if( lims == "points" )
    {   #
        lims <- c( 
            range( xy.data[,"xpos"] ), 
            range( xy.data[,"ypos"] )
        )   #
    }else{ 
        lims <- c( 
            range( xy.bound[,"xpos"] ), 
            range( xy.bound[,"ypos"] )
        )   #
    }   #
    #
    dens.xy <- kde2d( 
        x       = xy.data[,"xpos"], 
        y       = xy.data[,"ypos"], 
        n       = n, 
        lims    = lims  
    )   #
    #
    ex.dens.xy <- expand.grid( 
        "x" = dens.xy[["x"]], 
        "y" = dens.xy[["y"]]  
    )   #
    
    # require( "sp" ) 
    
    PiP <- as.logical(  point.in.polygon( 
        point.x = ex.dens.xy[,"x"],
        point.y = ex.dens.xy[,"y"], 
        pol.x   = xy.bound[,"xpos"], 
        pol.y   = xy.bound[,"ypos"]  
        #mode.checked    = FALSE  #  Deprecated
    )   )   #
    #
    dens.xy[["z"]][ !PiP ] <- NA 
    #
    return( dens.xy )  
}   #

#     rand.text <- TT.dataset(n=100,seed.val=1980042401) 

#     geo <- TT.geo.get( class.sys = "USDA.TT" )

#     kde.res <- TT.kde2d( 
#         geo         = geo, 
#         tri.data    = rand.text
#     )   #






TT.iwd <- function(# Inverse weighted distance interpolation on a grid. 
### Inverse weighted distance interpolation on a grid. 

 tri.data, 
 z.name, 
 geo, 
 
 css.names       = NULL, 
 tri.pol.data    = NULL,     # edges data, same format as tri.data
 
 text.tol        = NULL, 
 text.sum        = NULL, 
 blr.clock       = NULL, 
 
 tri.sum.tst     = NULL, 
 tri.pos.tst     = NULL, 
 #bg             = NULL, 
 
 set.par         = FALSE, 
 n               = 25, 
 lims            = c("points","triangle")[1], 
 max.dist        = NULL, 
 q.max.dist      = 0.5, 
 pow             = 0.5  

){  # 
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic takes the values in geo and 
    # attributes them to same name variables
    TT.geo.set(
        geo     = geo  
        #p.env  = environment()  
    )   # 
    #
    # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # Automatic sets remaining NULL varaibles 
    TT.auto.set( set.par = set.par ) 
    #
    xy.bound <- TT.get("base.frame") * text.sum 
    colnames( xy.bound ) <- css.names 
    #
    xy.bound <- TT.css2xy( 
        tri.data    = xy.bound,   
        geo         = geo, 
        text.tol    = text.tol, 
        css.names   = css.names, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        set.par     = set.par, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    xy.data <- TT.css2xy( 
        tri.data    = tri.data, 
        geo         = geo, 
        text.tol    = text.tol, 
        css.names   = css.names, 
        tri.sum.tst = tri.sum.tst, 
        tri.pos.tst = tri.pos.tst, 
        set.par     = set.par, 
        text.sum    = text.sum, 
        blr.clock   = blr.clock  
    )   #
    #
    if( lims == "points" )
    {   #
        xy.grid <- TT.xy.grid( 
            x   = xy.data[,"xpos"], 
            y   = xy.data[,"ypos"], 
            n   = n  
        )   #
    }else{ 
        xy.grid <- TT.xy.grid( 
            x   = xy.bound[,"xpos"], 
            y   = xy.bound[,"ypos"], 
            n   = n  
        )   #
    }   #
    #
    n.text <- dim(tri.data)[1] 
    #
    xy.data <- data.frame( 
        xy.data, 
        "z" = tri.data[,z.name]  
    )   #
    #
    xy.grid[["xypos"]] <- data.frame( 
        xy.grid[["xypos"]], 
        "z" = NA  
    )   #
    #
    dist.mat    <- dist( 
        x       = rbind( 
            xy.data[,c("xpos","ypos")], 
            xy.grid[["xypos"]][,c("xpos","ypos")] 
        ),  # 
        method  = "euclidean", 
        diag    = FALSE, 
        upper   = FALSE  
    )   #
    #
    if( is.null( max.dist ) )
    {   #
        max.dist    <- as.numeric( quantile( x = dist.mat, probs = q.max.dist ) ) 
    }   #
    #
    dist.mat    <- as.matrix( dist.mat )
    dist.mat    <- dist.mat[ (n.text+1):dim(dist.mat)[1] , 1:n.text ]
    dist.mat[ dist.mat > max.dist ] <- NA
    #
    xy.grid[["xypos"]][,"z"] <- unlist( apply( 
            X       = dist.mat, 
            MARGIN  = 1, 
            FUN     = function(X,z.vec){ 
                sel <- !is.na(X)
                if( any(sel) )
                {   #
                    w.z <- 1/(X^pow)
                    sum( z.vec * w.z, na.rm = TRUE )/sum( w.z, na.rm = TRUE ) 
                }else{ NA } 
            },  #
            z.vec   = xy.data[,"z"] 
    )   )   #
    
    # require( "sp" ) 
    
    PiP <- as.logical(  point.in.polygon( 
        point.x = xy.grid[["xypos"]][,"xpos"],
        point.y = xy.grid[["xypos"]][,"ypos"], 
        pol.x   = xy.bound[,"xpos"], 
        pol.y   = xy.bound[,"ypos"]  
        #mode.checked    = FALSE  #  Deprecated
    )   )   #
    #
    xy.grid[["xypos"]][ !PiP, "z" ] <- NA 
    #
    return( 
        list( 
            "x" = xy.grid[["xpos"]], 
            "y" = xy.grid[["ypos"]], 
            "z" = matrix( 
                data    = xy.grid[["xypos"]][, "z" ], 
                nrow    = length(xy.grid[["ypos"]]), 
                ncol    = length(xy.grid[["xpos"]]),
                byrow   = FALSE 
            )   #
        )   #
    )   #
}   #

#     rand.text <- TT.dataset(n=100,seed.val=1980042401) 

#     geo <- TT.geo.get( class.sys = "USDA.TT" )

#     izd.res <- TT.iwd( 
#         geo         = geo, 
#         tri.data    = rand.text, 
#         z.name      = "Z"  
#     )   #






TT.contour <- function(# Wrapper for the contour() function adapted to texture triangles.
### A wrapper for the contour() function 
### adapted to texture triangles (plot preparation).
### designed to plot the results of TT.mahalanobis() or 
### TT.kde2d(), before or after plot.

    geo, 
    x,  #  passed to contour
    add             = FALSE,  #  also passed to contour
    tri.sum.tst     = NULL, 
    tri.pos.tst     = NULL, 
    text.tol        = NULL, 
    unit.ps         = NULL, 
    unit.tx         = NULL, 
    b.lim           = NULL, 
    l.lim           = NULL, 
    main            = NULL, 
    new.mar         = NULL, 
    bg              = NULL, 
    fg              = NULL, 
    col             = NULL, 
    cex.main        = NULL, 
    lang            = NULL, 
    #
    # contour() parameters
    nlevels         = 10, 
    levels          = NA, 
    labels          = NULL, 
    xlim            = NA, 
    ylim            = NA, 
    zlim            = NA, 
    labcex          = 1, 
    drawlabels      = TRUE, 
    method          = "flattest", 
    #vfont, 
    axes            = TRUE, 
    frame.plot      = NA, 
    #col            = par("fg"), 
    lty             = NA, 
    lwd             = NA, 
    blr.clock       = NULL, 
    tlr.an          = NULL, 
    blr.tx          = NULL, 
    text.sum        = NULL, 
    base.css.ps.lim = NULL, 
    ...
){  #
    if( !is.list( x ) )
    {   #
        stop( 
            paste( 
                sep = "", 
                "x must be a list (now ", 
                class(x), ")"   
            )   #
        )   #
    }   #
    #
    if( !all(names(x) %in% c("x","y","z")) )
    {   #
        stop( 
            paste( 
                sep = "", 
                "x names must be x, y and z (now ", 
                paste(names(x),collapse=", "), ")"   
            )   #
        )   #
    }   #
    #
    # # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # # Automatic takes the values in geo and 
    # # attributes them to same name variables
    # TT.geo.set(
    #     geo     = geo  
    #     #p.env  = environment()  
    # )   # 
    # #
    # # +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ 
    # # Automatic sets remaining NULL variables 
    # TT.auto.set( set.par = TRUE ) 
    #
    class.sys <- "none" 
    #
    if( !add )
    {   #
        TT.baseplot( 
            geo             = geo, 
            class.sys       = class.sys, 
            # GEO PARAMETERS
            blr.clock       = blr.clock, 
            tlr.an          = tlr.an, 
            blr.tx          = blr.tx, 
            text.sum        = text.sum, 
            base.css.ps.lim = base.css.ps.lim, 
            # END GEO
            tri.sum.tst     = tri.sum.tst, 
            tri.pos.tst     = tri.pos.tst, 
            text.tol        = text.tol, 
            unit.ps         = unit.ps, 
            unit.tx         = unit.tx, 
            b.lim           = b.lim, 
            l.lim           = l.lim, 
            main            = main, 
            new.mar         = new.mar, 
            bg              = bg, 
            fg              = fg, 
            col             = col, 
            cex.main        = cex.main, 
            lang            = lang  
        )   #
    }   #
    #
    if( any( is.na( xlim       )) ){ xlim       <- range( x[["x"]], finite = TRUE) }  
    if( any( is.na( ylim       )) ){ ylim       <- range( x[["y"]], finite = TRUE) }  
    if( any( is.na( zlim       )) ){ zlim       <- range( x[["z"]], finite = TRUE) }  
    if( any( is.na( levels     )) ){ levels     <- pretty(zlim, nlevels) }  
    if( any( is.na( frame.plot )) ){ frame.plot <- axes }  
    if( any( is.na( lty        )) ){ lty        <- par("lty") }  
    if( any( is.na( lwd        )) ){ lwd        <- par("lwd") }  
    #
    contour( 
        x           = x, 
        add         = TRUE, 
        nlevels     = nlevels, 
        levels      = levels, 
        labels      = labels, 
        xlim        = xlim, 
        ylim        = ylim, 
        zlim        = zlim, 
        labcex      = labcex, 
        drawlabels  = drawlabels, 
        method      = method, 
        axes        = axes, 
        frame.plot  = frame.plot, 
        lty         = lty, 
        lwd         = lwd, 
        col         = col, 
        ...
    )   #
}   #

#     rand.text <- TT.dataset(n=100,seed.val=1980042401) 

#     geo <- TT.geo.get( class.sys = "USDA.TT" )

#     maha <- TT.mahalanobis( 
#         geo         = geo, 
#         tri.data    = rand.text
#     )   #

#     TT.contour( x = maha, geo = geo, add = FALSE, lwd = 2 ) 

#     TT.plot( 
#         class.sys   = "FR.GEPPA.TT",
#         geo         = geo, 
#         grid.show   = FALSE, 
#         add         = TRUE  
#    )    #



#     rand.text <- TT.dataset(n=100,seed.val=1980042401) 

#     geo <- TT.geo.get( class.sys = "USDA.TT" )

#     kde.res <- TT.kde2d( 
#         geo         = geo, 
#         tri.data    = rand.text
#     )   #

#     TT.contour( x = kde.res, geo = geo, add = FALSE, lwd = 2 ) 

#     TT.plot( 
#         class.sys   = "HYPRES.TT",
#         geo         = geo, 
#         grid.show   = FALSE, 
#         add         = TRUE  
#    )    #






TT.image <- function(# Wrapper for the contour() function adapted to texture triangles.
### A wrapper for the contour() function 
### adapted to texture triangles (plot preparation).
### designed to plot the results of TT.mahalanobis() or 
### TT.kde2d() [to be written], before or after plot.

    geo, 
    x,  #  passed to contour
    add             = FALSE,  #  also passed to contour
    tri.sum.tst     = NULL, 
    tri.pos.tst     = NULL, 
    text.tol        = NULL, 
    unit.ps         = NULL, 
    unit.tx         = NULL, 
    b.lim           = NULL, 
    l.lim           = NULL, 
    main            = NULL, 
    new.mar         = NULL, 
    bg              = NULL, 
    fg              = NULL, 
    cex.main        = NULL, 
    lang            = NULL, 
    #
    # contour() parameters
    xlim            = NA, 
    ylim            = NA, 
    zlim            = NA, 
    col             = rev( heat.colors(12) ),
    oldstyle        = FALSE, 
    blr.clock       = NULL, 
    tlr.an          = NULL, 
    blr.tx          = NULL, 
    text.sum        = NULL, 
    base.css.ps.lim = NULL, 
 ...
### Additional parameters passed to image().
){  #
    if( !is.list( x ) )
    {   #
        stop( 
            paste( 
                sep = "", 
                "x must be a list (now ", 
                class(x), ")"   
            )   #
        )   #
    }   #
    #
    if( !all(names(x) %in% c("x","y","z")) )
    {   #
        stop( 
            paste( 
                sep = "", 
                "x names must be x, y and z (now ", 
                paste(names(x),collapse=", "), ")"   
            )   #
        )   #
    }   #
    #
    class.sys <- "none"  
    #
    if( !add )
    {   #
        TT.baseplot( 
            geo             = geo, 
            class.sys       = class.sys, 
            # GEO PARAMETERS
            blr.clock       = blr.clock, 
            tlr.an          = tlr.an, 
            blr.tx          = blr.tx, 
            text.sum        = text.sum, 
            base.css.ps.lim = base.css.ps.lim, 
            # END GEO
            tri.sum.tst     = tri.sum.tst, 
            tri.pos.tst     = tri.pos.tst, 
            text.tol        = text.tol, 
            unit.ps         = unit.ps, 
            unit.tx         = unit.tx, 
            b.lim           = b.lim, 
            l.lim           = l.lim, 
            main            = main, 
            new.mar         = new.mar, 
            bg              = bg, 
            fg              = fg, 
            col             = col, 
            cex.main        = cex.main, 
            lang            = lang  
        )   #
    }   #
    #
    if( any( is.na( xlim       )) ){ xlim       <- range( x[["x"]], finite = TRUE) }  
    if( any( is.na( ylim       )) ){ ylim       <- range( x[["y"]], finite = TRUE) }  
    if( any( is.na( zlim       )) ){ zlim       <- range( x[["z"]], finite = TRUE) }  
    #
    image( 
        x           = x, 
        add         = TRUE, 
        xlim        = xlim, 
        ylim        = ylim, 
        zlim        = zlim, 
        col         = col,
        oldstyle    = oldstyle, 
        ...  
    )   #
}   #

#     rand.text <- TT.dataset(n=100,seed.val=1980042401) 

#     geo <- TT.geo.get( class.sys = "USDA.TT" )

#     maha <- TT.mahalanobis( 
#         geo         = geo, 
#         tri.data    = rand.text, 
#         n           = 50, 
#         divisorvar  = 2  
#     )   #

#     TT.contour( x = maha, geo = geo, add = FALSE, lwd = 6, col = "red" ) 

#     maha <- TT.mahalanobis( 
#         geo         = geo, 
#         tri.data    = rand.text, 
#         n           = 50, 
#         divisorvar  = 1  
#     )   #

#     TT.contour( x = maha, geo = geo, add = TRUE, lwd = 3, col = "green" ) 

#     maha <- TT.mahalanobis( 
#         geo         = geo, 
#         tri.data    = rand.text, 
#         n           = 50, 
#         divisorvar  = 3  
#     )   #

#     TT.contour( x = maha, geo = geo, add = TRUE, lwd = 1, col = "blue" ) 

#     TT.plot( 
#         class.sys   = "FR.GEPPA.TT",
#         geo         = geo, 
#         grid.show   = FALSE, 
#         add         = TRUE  
#    )    #



#     rand.text <- TT.dataset(n=100,seed.val=1980042401) 

#     geo <- TT.geo.get( class.sys = "USDA.TT" )

#     maha <- TT.mahalanobis( 
#         geo         = geo, 
#         tri.data    = rand.text, 
#         n           = 50, 
#         divisorvar  = 2, 
#         alr         = TRUE 
#     )   #

#     TT.contour( 
#         x = maha, 
#         geo = geo, 
#         add = FALSE, 
#         lwd = 6, 
#         col = "red", 
#         levels  = c(0.5,1,2,4,8) 
#     )   #

#     maha <- TT.mahalanobis( 
#         geo         = geo, 
#         tri.data    = rand.text, 
#         n           = 50, 
#         divisorvar  = 1, 
#         alr         = TRUE 
#     )   #

#     TT.contour( 
#         x = maha, 
#         geo = geo, 
#         add = TRUE, 
#         lwd = 3, 
#         col = "green", 
#         levels  = c(0.5,1,2,4,8) 
#     )   #


#     maha <- TT.mahalanobis( 
#         geo         = geo, 
#         tri.data    = rand.text, 
#         n           = 50, 
#         divisorvar  = 3, 
#         alr         = TRUE 
#     )   #

#     TT.contour( 
#         x = maha, 
#         geo = geo, 
#         add = TRUE, 
#         lwd = 1, 
#         col = "blue", 
#         levels  = c(0.5,1,2,4,8) 
#     )   #

#     TT.plot( 
#         class.sys   = "FR.GEPPA.TT",
#         geo         = geo, 
#         grid.show   = FALSE, 
#         add         = TRUE  
#    )    #



#     rand.text <- TT.dataset(n=100,seed.val=1980042401) 

#     geo <- TT.geo.get( class.sys = "USDA.TT" )

#     kde.res <- TT.kde2d( 
#         geo         = geo, 
#         tri.data    = rand.text, 
#         n           = 50  
#     )   #

#     TT.image( x = kde.res, geo = geo, add = FALSE, col = hsv(h=0.21,s=0.9,v=seq(0.2,1,length.out=12)))  
#     TT.contour( x = kde.res, geo = geo, add = TRUE, lwd = 2 ) 

#     TT.plot( 
#         class.sys   = "HYPRES.TT",
#         geo         = geo, 
#         grid.show   = FALSE, 
#         add         = TRUE  
#    )    #



#     rand.text <- TT.dataset(n=100,seed.val=1980042401) 

#     geo <- TT.geo.get( class.sys = "USDA.TT" )

#     iwd.res <- TT.iwd( 
#         geo         = geo, 
#         tri.data    = rand.text, 
#         z.name      = "Z", 
#         pow         = 0.5, 
#         q.max.dist  = 0.50  
#     )   #

#     TT.image( x = iwd.res, geo = geo, add = FALSE ) 
#     TT.contour( x = iwd.res, geo = geo, add = TRUE, lwd = 2 ) 

#     TT.plot( 
#         class.sys   = "HYPRES.TT",
#         geo         = geo, 
#         grid.show   = FALSE, 
#         add         = TRUE  
#    )    #






TT.normalise.sum <- function(# Normalises the sum of the 3 particle size classes.
### Normalises the sum of the 3 particle size classes in tri.data 
### to text.sum (100%).
    tri.data, 
    css.names   = NULL, 
    text.sum    = NULL, 
    text.tol    = NULL, 
    #
    tri.pos.tst = NULL, 
    residuals   = FALSE  
){  #
    # Set rest of variables:
    TT.auto.set( set.par = FALSE ) 
    #
    TT.data.test( 
        tri.data    = tri.data, 
        css.names   = css.names, 
        text.sum    = text.sum, 
        text.tol    = text.tol, 
        tri.sum.tst = FALSE, 
        tri.pos.tst = tri.pos.tst  
    )   #
    #
    class.ini <- class( tri.data ) 
    #
    sum.ini <- apply( 
        X       = tri.data[,css.names], 
        MARGIN  = 1, 
        FUN     = sum  
    )   #
    #
    tri.data <- t( apply( 
        X       = tri.data[,css.names], 
        MARGIN  = 1, 
        FUN     = function(X){ 
            X * (text.sum/sum(X)) 
        }   #
    )   )   #
    #
    if( class.ini == "data.frame" ) 
    {   #
        tri.data <- as.data.frame( tri.data )
    }   #
    #
    #class( tri.data ) <- class.ini 
    #
    sum.end <- apply( 
        X       = tri.data[,css.names], 
        MARGIN  = 1, 
        FUN     = sum  
    )   #
    #
    if( residuals )
    {   #
        tri.data <- cbind( 
            tri.data, 
            "residuals" = sum.ini - sum.end  
        )   #
    }   #
    #
    return( tri.data )
}   #

#     rand.text <- TT.dataset(n=100,seed.val=1980042401) 

#     TT.normalise.sum("tri.data"=rand.text,"residuals"=TRUE)[1:20,] 






TT.normalise.sum.X <- function(# Normalises the sum of the X particle size classes.
### Normalises the sum of the X particle size classes
### in tri.data to text.sum (100%).

    tri.data, 
    text.sum    = NULL, 
    text.tol    = NULL, 
    #
    tri.pos.tst = NULL, 
    residuals   = FALSE  
){  #
    # Set rest of variables:
    TT.auto.set( set.par = FALSE ) 
    #
    TT.data.test.X( 
        tri.data    = tri.data, 
        text.sum    = text.sum, 
        text.tol    = text.tol, 
        tri.sum.tst = FALSE, 
        tri.pos.tst = tri.pos.tst  
    )   #
    sum.ini <- apply( 
        X       = tri.data[,], 
        MARGIN  = 1, 
        FUN     = sum  
    )   #
    #
    print(sum.ini)
    tri.data <- t( apply( 
        X       = tri.data[,], 
        MARGIN  = 1, 
        FUN     = function(X){ 
            X * (text.sum/sum(X)) 
        }   #
    )   )   #
    sum.end <- apply( 
        X       = tri.data[,], 
        MARGIN  = 1, 
        FUN     = sum  
    )   #
    #
    if( residuals )
    {   #
        tri.data <- cbind( 
            tri.data, 
            "residuals" = sum.ini - sum.end  
        )   #
    }   #
    #
    return( tri.data )
}   #

#    my.text4 <- data.frame(
#        "CLAY"  = c(05,60,15,04.9,25,05,25,45,65,75,13,47),
#        "FSILT" = c(02,04.3,10,15,25,40,35,20,10,05,10,20),
#        "CSILT" = c(03,04,05,10,30,45,30,25,05,10,07.2,23.3),
#        "SAND"  = c(90.5,32,70,70,20.3,10.9,9.3,9.4,20,10,70,10)
#    )   #
#    res <- TT.normalise.sum.X(
#        tri.data    = my.text4,
#        residuals   = TRUE
#    )   #
#    res


# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# | END OF THE R CODE                   |
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+



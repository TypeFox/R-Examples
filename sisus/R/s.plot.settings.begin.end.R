s.plot.settings.begin.end <-
function# Plot settings for all plots
### internal function for sisus
(filename.prefix = "SISUS"
### internal variable
, plot.filename = "default_filename"
### internal variable
, plot.mode = "settings"
### internal variable
, plot.format = 1
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

#
#   if modifying plot file types, need to update:
#     Excel workbook, for input settings
#     assign.variables.R, for reading workbook input and setting plot.format.list
#     s.plot.settings.begin.end.R, this function for settings
#
  ##library("Cairo"); # for raster formats under unix

  # determine if running in windows or unix environment
  #OS = .Platform$OS.type;

  ################################################################################
  if (plot.mode == "settings")
  {
    # filenames
    #file.temp     = paste(filename.prefix, "_", plot.filename, "_%02d", sep="");
    file.temp     = paste(filename.prefix, "_", plot.filename, sep="");  # no "_%02d" for convert to work
    file.png      = paste(file.temp, ".png" , sep="");
    file.eps      = paste(file.temp, ".eps" , sep="");
    file.pdf      = paste(file.temp, ".pdf" , sep="");
    file.bmp      = paste(file.temp, ".bmp" , sep="");
    file.jpeg     = paste(file.temp, ".jpeg", sep="");
    file.tiff     = paste(file.temp, ".tiff", sep="");    #(added at R 2.7)
    # file.gif      = paste(file.temp, ".gif" , sep="");  #(removed at R 2.7)

    convert.path  = "/usr/bin/convert";  # for converting from png to other raster formats

    # devices settings
    onefile       = TRUE;  # FALSE
    family        = "ComputerModern";
    encoding      = "TeXtext.enc";
    width.in      = 8.0;
    height.in     = 8.0;
    resolution    = 300;
    width.px      = 800;
    height.px     = 800;
    #width.px      = width.in *resolution; #800;
    #height.px     = height.in*resolution; #800;
    horizontal    = FALSE;
    paper         = "letter";
    pagecentre    = TRUE;
    print.it      = FALSE;
    quality.jpeg  = 100;
    compress.tiff = "none";

    # par settings
    par.oma       = c(5,4,4,4);
    par.mar       = c(4,4,2,2);
    par.mfrow     = c(1,1);
    par.bg        = "transparent";
    par.xpd       = TRUE;


    # list
    PLOT.SETTINGS = new.env();  # create a list to return with all data

      # filenames
    PLOT.SETTINGS$file.png      = file.png     ;
    PLOT.SETTINGS$file.eps      = file.eps     ;
    PLOT.SETTINGS$file.pdf      = file.pdf     ;
    PLOT.SETTINGS$file.bmp      = file.bmp     ;
    PLOT.SETTINGS$file.jpeg     = file.jpeg    ;
    # PLOT.SETTINGS$file.gif      = file.gif     ;
    PLOT.SETTINGS$file.tiff     = file.tiff    ;

    PLOT.SETTINGS$convert.path  = convert.path ;

      # pdf and ps parameters
    PLOT.SETTINGS$onefile       = onefile      ;
    PLOT.SETTINGS$family        = family       ;
   #PLOT.SETTINGS$title         = title        ;
   #PLOT.SETTINGS$fonts         = fonts        ;
    PLOT.SETTINGS$encoding      = encoding     ;
   #PLOT.SETTINGS$bg            = bg           ;
   #PLOT.SETTINGS$fg            = fg           ;
    PLOT.SETTINGS$width.in      = width.in     ;
    PLOT.SETTINGS$height.in     = height.in    ;
    PLOT.SETTINGS$horizontal    = horizontal   ;
   #PLOT.SETTINGS$pointsize     = pointsize    ;
    PLOT.SETTINGS$paper         = paper        ;
    PLOT.SETTINGS$pagecentre    = pagecentre   ;
    PLOT.SETTINGS$print.it      = print.it     ;
   #PLOT.SETTINGS$command       = command      ;
   #PLOT.SETTINGS$colormodel    = colormodel   ;

      # rasters
    PLOT.SETTINGS$width.px      = width.px     ;
    PLOT.SETTINGS$height.px     = height.px    ;
    PLOT.SETTINGS$resolution    = resolution;

      # jpeg
    PLOT.SETTINGS$quality.jpeg  = quality.jpeg ;
      # tiff


    PLOT.SETTINGS$par.oma       = par.oma      ;
    PLOT.SETTINGS$par.mar       = par.mar      ;
    PLOT.SETTINGS$par.mfrow     = par.mfrow    ;
    PLOT.SETTINGS$par.bg        = par.bg       ;
    PLOT.SETTINGS$par.xpd       = par.xpd      ;

    return( as.list(PLOT.SETTINGS) );
  } # settings

  ################################################################################
  if (plot.mode == "begin")
  {
    PLOT.SETTINGS = s.plot.settings.begin.end(filename.prefix, plot.filename);

    # to interactive display
    if (plot.format == 0) {
      #if (OS == "unix") { X11(); }
      #if (OS == "windows") { windows(); }
      dev.new();
    }

    # png is default output since quick display on web and used to convert to other raster formats
    if (plot.format == 1) {
      png       (filename   = PLOT.SETTINGS$file.png
               , width      = PLOT.SETTINGS$width.px
               , height     = PLOT.SETTINGS$height.px
      #         , res        = PLOT.SETTINGS$resolution
                 );
    }

    # ps
    if (plot.format == 2) {
      postscript(file       = PLOT.SETTINGS$file.eps   ,
                 onefile    = PLOT.SETTINGS$onefile    ,
      #           family     = PLOT.SETTINGS$family     ,
      #           encoding   = PLOT.SETTINGS$encoding   ,
                 width      = PLOT.SETTINGS$width.in   ,
                 height     = PLOT.SETTINGS$height.in  ,
                 horizontal = PLOT.SETTINGS$horizontal ,
                 paper      = PLOT.SETTINGS$paper      ,
                 pagecentre = PLOT.SETTINGS$pagecentre ,
                 print.it   = PLOT.SETTINGS$print.it
                 );
    }

    # pdf
    if (plot.format == 3) {
      pdf       (file       = PLOT.SETTINGS$file.pdf   ,
                 onefile    = PLOT.SETTINGS$onefile    ,
                 #family     = PLOT.SETTINGS$family     ,
                 #encoding   = PLOT.SETTINGS$encoding   ,
                 width      = PLOT.SETTINGS$width.in   ,
                 height     = PLOT.SETTINGS$height.in  ,
                # horizontal = PLOT.SETTINGS$horizontal ,
                 paper      = PLOT.SETTINGS$paper      ,
                 pagecentre = PLOT.SETTINGS$pagecentre
                # print.it   = PLOT.SETTINGS$print.it
                 );
    }

    # bmp
    if (plot.format == 4) {
      bmp       (filename   = PLOT.SETTINGS$file.bmp
               , width      = PLOT.SETTINGS$width.px
               , height     = PLOT.SETTINGS$height.px
      #         , res        = PLOT.SETTINGS$resolution
                 );

      ## Convert from png        #system("/usr/bin/convert test.png test.bmp")
      #convert.command = paste(PLOT.SETTINGS$convert.path, PLOT.SETTINGS$file.png, PLOT.SETTINGS$file.bmp, sep=" ");
      #system(convert.command);
    }

    # jpeg
    if (plot.format == 5) {
      jpeg      (filename   = PLOT.SETTINGS$file.jpeg
               , width      = PLOT.SETTINGS$width.px
               , height     = PLOT.SETTINGS$height.px
               , quality    = PLOT.SETTINGS$quality.jpeg
      #         , res        = PLOT.SETTINGS$resolution
                 );
    }

    ## gif
    #if (plot.format == 6) {
    #  # Convert from png        #system("/usr/bin/convert test.png test.gif")
    #  convert.command = paste(PLOT.SETTINGS$convert.path, PLOT.SETTINGS$file.png, PLOT.SETTINGS$file.gif, sep=" ");
    #  system(convert.command);
    #}

    # tiff
    if (plot.format == 6) {
      tiff      (filename    = PLOT.SETTINGS$file.tiff
               , width       = PLOT.SETTINGS$width.px
               , height      = PLOT.SETTINGS$height.px
               , compression = PLOT.SETTINGS$compress.tiff
      #         , res        = PLOT.SETTINGS$resolution
                 );
    }


    #dev.control(displaylist = "enable");
    #dev.set(which = dev.next());

    # par settings
    par(oma = PLOT.SETTINGS$par.oma,
        mar = PLOT.SETTINGS$par.mar,
      mfrow = PLOT.SETTINGS$par.mfrow,
         bg = PLOT.SETTINGS$par.bg,
        xpd = PLOT.SETTINGS$par.xpd);

  } # begin

  ################################################################################
  if (plot.mode == "print")
  {
    ## automatically             # ps
    ##dev.copy(which = dev.cur());  #  device = postscript  );   # ps
    #dev.copy(which = dev.next());  #  device = pdf );          # pdf
    #dev.copy(which = dev.next());  #  device = bmp );          # bmp
    #dev.copy(which = dev.next());  #  device = jpeg);          # jpeg
    #dev.copy(which = dev.next());  #  device = png );          # png
    #dev.copy(which = dev.next());  #  device = tiff);          # tiff
    #dev.set(which = dev.next()); # loop to ps
  } # print

  ################################################################################
  if (plot.mode == "end")
  {
    if (plot.format != 0) { dev.off(); };
    #graphics.off(); # don't use
  } # end

  ### internal variable
}

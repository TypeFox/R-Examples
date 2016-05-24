# Function to extract outlines from PostScript paths

PScaptureText <- c()

PScaptureChars <- c()
    
PScaptureHead <- function(file, charpath, charpos, setflat, encoding) {
    c("%!PS-Adobe-2.0 EPSF-1.2",
      "%%BeginProcSet:convertToR 0 0",
      
      # XML file header info
      sprintf("(<?xml version='1.0' encoding='%s' ?>\n\n) print",
              encoding),
      paste("(<picture version='3' xmlns:rgml='http://r-project.org/RGML'",
            " source='", file, "'",
            " date='", as.character(Sys.time()), "'",
            " creator='R (", R.Version()$major, ".", R.Version()$minor, ")'",
            " >\n\n) print", sep=""),
      # useful definitions
      # Define my counters GLOBALLY so the postscript
      # file I run cannot reset my counters with a restore!
      "true setglobal",
      "/convertToR dup 100 dict def load begin",
      "/str 50 string def",
      "/id 1 def",
      # Do these need to be even "larger"?
      "/xmax -99999 def",
      "/xmin  99999 def",
      "/ymax -99999 def",
      "/ymin  99999 def",

      # Start with x, y and end with dw
      "/transwidth {",
      "  transform",
      "  0 0 transform",
      "  2 index exch sub", # dy
      "  3 index 2 index sub", # dx
      "  2 exp exch 2 exp add sqrt",
      "  3 1 roll pop pop pop",
      "  } def",
      
      # bounding box processing
      "/boxinit {",
      "  convertToR /bxmax -99999 put",
      "  convertToR /bxmin  99999 put",
      "  convertToR /bymax -99999 put",
      "  convertToR /bymin  99999 put",
      "  } def",
      "/boxmove {",
      "  transform",
      "  dup",
      "  convertToR exch /cury exch put",
      "  dup",
      "  convertToR /bymin get lt {convertToR /bymin cury put} if",
      "  convertToR /bymax get gt {convertToR /bymax cury put} if",
      "  dup",
      "  convertToR exch /curx exch put",
      "  dup",
      "  convertToR /bxmin get lt {convertToR /bxmin curx put} if",
      "  convertToR /bxmax get gt {convertToR /bxmax curx put} if",
      "  } def",
      "/boxline {",
      "  transform",
      "  dup",
      "  convertToR exch /cury exch put",
      "  dup",
      "  convertToR /bymin get lt {convertToR /bymin cury put} if",
      "  convertToR /bymax get gt {convertToR /bymax cury put} if",
      "  dup",
      "  convertToR exch /curx exch put",
      "  dup",
      "  convertToR /bxmin get lt {convertToR /bxmin curx put} if",
      "  convertToR /bxmax get gt {convertToR /bxmax curx put} if",
      "  } def",
      "/boxcurve {",
      " } def",
      "/boxclose {",
      " } def",

      # Update global bbox from bbox
      "/boxupdate {",
      "  convertToR /bymin get convertToR /ymin get lt {convertToR /ymin convertToR /bymin get put} if",
      "  convertToR /bymax get convertToR /ymax get gt {convertToR /ymax convertToR /bymax get put} if",
      "  convertToR /bxmin get convertToR /xmin get lt {convertToR /xmin convertToR /bxmin get put} if",
      "  convertToR /bxmax get convertToR /xmax get gt {convertToR /xmax convertToR /bxmax get put} if",
      " } def",

      # Font height calculation (for *current* font)
      "/fontsize {",
      # DO NOT just 'get' the third element of FontMatrix
      # because the font may, e.g., rotate the font outlines
      "  currentfont /FontMatrix get",
      # For Type 1 fonts, scale the matrix by 1000
      # (see PLRM.pdf version 3 page 324; comment on FontMatrix)
      # For Type 42 (TrueType), no multiplier necessary
      # (http://www.adobe.com/devnet/font/pdfs/5012.Type42_Spec.pdf)
      # For other font types, not sure what the initial font matrix is
      # COPY the FontMatrix to avoid 'invalidaccess'
      "  currentfont /FontType get 1 eq { 1000 1000 matrix scale matrix concatmatrix } if",
      # Transform (0, 1) using this matrix, then transform to user coords
      "  dup 0 exch 1 exch transform transform",
      # Do the same to (0, 0)
      "  2 index 0 exch 0 exch transform transform",
      # dy
      "  2 index exch sub",
      # dx
      "  3 index 2 index sub",
      "  2 exp exch 2 exp add sqrt",
      # clean up
      "  5 1 roll pop pop pop pop",
      " } def",
      
      # path processing
      "/mymove {",
      "  (\t<move) print",
      "  transform",
      "  dup",
      "  convertToR exch /ystart exch put",
      "  dup",
      "  convertToR exch /cury exch put",
      "  dup",
      "  convertToR /ymin get lt {convertToR /ymin cury put} if",
      "  dup",
      "  convertToR /ymax get gt {convertToR /ymax cury put} if",
      "  ( y=') print str cvs print (') print",
      "  dup",
      "  convertToR exch /xstart exch put",
      "  dup",
      "  convertToR exch /curx exch put",
      "  dup",
      "  convertToR /xmin get lt {convertToR /xmin curx put} if",
      "  dup",
      "  convertToR /xmax get gt {convertToR /xmax curx put} if",
      "  ( x=') print str cvs print (') print",
      "  (/>\n) print",
      "} def",
      "/myline {",
      "  (\t<line) print",
      "  transform",
      "  dup",
      "  convertToR exch /cury exch put",
      "  dup",
      "  convertToR /ymin get lt {convertToR /ymin cury put} if",
      "  dup",
      "  convertToR /ymax get gt {convertToR /ymax cury put} if",
      "  ( y=') print str cvs print (') print",
      "  dup",
      "  convertToR exch /curx exch put",
      "  dup",
      "  convertToR /xmin get lt {convertToR /xmin curx put} if",
      "  dup",
      "  convertToR /xmax get gt {convertToR /xmax curx put} if",
      "  ( x=') print str cvs print (') print",
      "  (/>\n) print",
      "  } def",
      "/mycurve {",
      "  } def",
      "/myclose {",
      # Convert 'closepath' to 'lineto'
      # "  (\t<close/>\n) print",
      "  (\t<line) print",
      "  ( y=') print convertToR /ystart get str cvs print (') print",
      "  ( x=') print convertToR /xstart get str cvs print (') print",
      "  (/>\n) print",
      "  } def",
      
      # echoing graphics state
      "/printcol {",
      # What colorspace are we in?
      # Return colorspace array and extract
      # first element which is colorspace name
      "  currentcolorspace 0 get",
      # If it's DeviceRGB or DeviceGray or DeviceCMYK we're ok
      # Separation and Pattern are special cases that are currently handled
      # by just specifying "grey" as the colour
      # Otherwise we "hail mary" and hope that currencolor throws back 3 values
      # that can be interpreted as RGB (e.g., R's sRGB!)
      "  dup dup (Separation) eq exch (Pattern) eq or {pop 0.5 0.5 0.5} {dup (DeviceGray) eq exch dup (DeviceRGB) eq exch (DeviceCMYK) eq or or {currentrgbcolor} {currentcolor} ifelse} ifelse",
      "  (\t\t<rgb) print",
      # make sure the colour is RGB not BGR
      "  ( r=') print 2 index str cvs print (') print",
      "  ( g=') print 1 index str cvs print (') print",
      "  ( b=') print str cvs print (') print",
      "  (/>\n) print",
      "  pop pop",
      "  } def",
      "/printlwd {",
      # lwd is in user coords so transform
      # This will need transforming to a grDevices "lwd" in R
      "  currentlinewidth 0 transform",
      "  0 0 transform",
      "  3 index 2 index sub", # dx
      "  3 index 2 index sub", # dy
      "  2 exp exch 2 exp add sqrt",      
      # clean up
      "  5 1 roll pop pop pop pop",
      "  ( lwd=') print str cvs print (') print",
      "} def",
      "/printlineend {",
      "  ( lineend=') print currentlinecap str cvs print (') print",
      "} def",
      "/printlinejoin {",
      "  ( linejoin=') print currentlinejoin str cvs print (') print",
      "} def",
      "/printlinemiter {",
      "  ( linemitre=') print currentmiterlimit str cvs print (') print",
      "} def",
      "/printdash {",
      "  currentdash",
      # Like lwd, transform user coords
      # Ignore offset
      "  ( lty=') print pop {0 transform pop 0 0 transform pop sub str cvs print ( ) print} forall (') print",
      "} def",      
      "/printstyle {",
      "  (\t\t<style) print",
      "  printlwd",
      "  printdash",
      "  printlineend",
      "  printlinemiter",
      "  printlinejoin",
      "  (/>\n) print",
      "} def",

      # print out "closestroke" marker plus graphics state info
      "/mystroke {",
      "  (<path type='stroke') print",
      "  ( id=') print convertToR /id get str cvs print ('>\n) print",
      "  (\t<context>\n) print",
      "  printcol  ",
      "  printstyle",
      "  (\t</context>\n\n) print",
      "  pathforall",
      "  convertToR /id get 1 add convertToR exch /id exch put",
      "  (</path>\n\n) print",
      "} def",
      "/myfill {",
      "  (<path type='fill') print",
      "  ( id=') print convertToR /id get str cvs print ('>\n) print",
      "  (\t<context>\n) print",
      "  printcol",
      "  printstyle",
      "  (\t</context>\n\n) print",
      "  pathforall",
      "  convertToR /id get 1 add convertToR exch /id exch put",
      "  (</path>\n\n) print",
      "} def",
      "/myeofill {",
      "  (<path type='eofill') print",
      "  ( id=') print convertToR /id get str cvs print ('>\n) print",
      "  (\t<context>\n) print",
      "  printcol",
      "  printstyle",
      "  (\t</context>\n\n) print",
      "  pathforall",
      "  convertToR /id get 1 add convertToR exch /id exch put",
      "  (</path>\n\n) print",
      "} def",
      "/mytext {",
      "  (<text ) print",
      "  ( id=') print convertToR /id get str cvs print (') print",
      "  convertToR /id get 1 add convertToR exch /id exch put",
      # Insert MARKS here so that postProcess() can easily locate strings
      "  ( string='###TEXT) print dup print (###TEXT') print",
      # Type of text element
      if (charpath) {
          "  ( type='charpath') print"
      } else if (charpos) {
          "  ( type='char') print"
      } else {
          "  ( type='text') print"
      },
      # Calculate any width adjustments
      "  convertToR /widthadjfun get cvx exec",
      # (x, y) location of text
      "  currentpoint",
      "  transform",
      "    dup",
      "    convertToR exch /cury exch put",
      "    dup",
      "    convertToR /ymin get lt {convertToR /ymin cury put} if",
      "    dup",
      "    convertToR /ymax get gt {convertToR /ymax cury put} if",
      "    ( y=') print str cvs print (') print",
      "    dup",
      "    convertToR exch /curx exch put",
      "    dup",
      "    convertToR /xmin get lt {convertToR /xmin curx put} if",
      "    dup",
      "    convertToR /xmax get gt {convertToR /xmax curx put} if",
      "    ( x=') print str cvs print (') print",
      # width of text 
      "  dup stringwidth",
      # stringwidth has put wx and wy on stack
      # If wy is non-zero, text is at an angle
      "  dup 0 ne { 1 index 1 index exch atan } { 0 } ifelse",
      # Save this angle
      "  convertToR exch /curangle exch put",
      # Calculate user-space width from wx and wy
      "  2 exp exch 2 exp add sqrt",      
      # Transform user-space width to device
      "  0 transform",
      # Subtract (0, 0)[x]
      "  exch 0 0 transform pop sub",
      # Subtract (0, 0)[y]
      "  exch 0 0 transform exch pop sub",
      # If transformed y is non-zero then text is at an angle 
      "  dup 0 ne { 1 index 1 index exch atan } { 0 } ifelse",
      # Save user-space angle plus angle
      "  convertToR /curangle get add convertToR exch /curangle exch put",
      # Record user-space angle plus angle
      "  ( angle=') print",
      "    convertToR /curangle get str cvs print (') print",
      # Calculate device width 
      "  2 exp exch 2 exp add sqrt",
      # Make any adjustments
      "  convertToR /curadj get add",
      # Print width
      "  ( width=') print",
      "    str cvs print (') print",
      # Height of text (size of current font)
      # Print height
      "  ( height=') print",
      "    fontsize str cvs print (') print",      
      # Calculate text bb
      "  convertToR /bboxfun get cvx exec",
      # Record bbox
      "  ( bbox=') print",
      "  convertToR /bxmin get str cvs print ( ) print",
      "  convertToR /bymin get str cvs print ( ) print",
      "  convertToR /bxmax get str cvs print ( ) print",
      "  convertToR /bymax get str cvs print (') print",
      # Record font (if available)
      "  currentfont /FontName known",
      "  { ( fontName=') print currentfont /FontName get str cvs print (') print} if",
      "  currentfont /FontInfo known",
      "  { currentfont /FontInfo get /FamilyName known",
      "    { ( fontFamilyName=') print currentfont /FontInfo get /FamilyName get print (') print } if } if",
      "  currentfont /FontInfo known",
      "  { currentfont /FontInfo get /FullName known",
      "    { ( fontFullName=') print currentfont /FontInfo get /FullName get print (') print } if } if",
      "  (>\n) print",
      # Graphics context
      "  (\t<context>\n) print",
      "  printcol",
      "  printstyle",
      "  (\t</context>\n\n) print",

      # If charpath or charloc, need to nest additional elements
      if (charpath) {
          "convertToR /charpathfun get cvx exec"
      } else if (charpos) {
          "convertToR /charfun get cvx exec"
      },
      
      "  (</text>\n\n) print",
      "} def",
      "/mychar {",
      "  (<path type='char') print",
      "  ( id=') print convertToR /id get str cvs print (') print",
      # Insert MARKS here so that postProcess() can easily locate strings
      "  ( char='###TEXT) print 4 index print (###TEXT' >\n) print",
      "  pathforall",
      "  (\t<context>\n) print",
      "  printcol",
      "  printstyle",
      "  (\t</context>\n\n) print",
      "  convertToR /id get 1 add convertToR exch /id exch put",
      "  (</path>\n\n) print",
      "} def",

      # override paint operators
      "/stroke {",
      "  flattenpath {mymove} {myline} {mycurve} {myclose}",
      "  mystroke",
      "  newpath",
      "} def",
      "/fill {",
      "  flattenpath {mymove} {myline} {mycurve} {myclose}",
      "  myfill",
      "  newpath",
      "} def",
      "/eofill {",
      "  flattenpath {mymove} {myline} {mycurve} {myclose}",
      "  myeofill",
      "  newpath",
      "} def",

      # Turn rectangles into paths
      "/rectfill {",
      "  4 dict begin",
      "  /rheight exch def",
      "  /rwidth exch def",
      "  /rx exch def",
      "  /ry exch def",
      "  rx ry moveto",
      "  rwidth 0 rlineto",
      "  0 rheight rlineto",
      "  rwidth neg 0 rlineto",
      "  end",
      "  closepath",
      "  fill",
      "} def",
      "/rectstroke {",
      "  4 dict begin",
      "  /rheight exch def",
      "  /rwidth exch def",
      "  /rx exch def",
      "  /ry exch def",
      "  rx ry moveto",
      "  rwidth 0 rlineto",
      "  0 rheight rlineto",
      "  rwidth neg 0 rlineto",
      "  end",
      "  closepath",
      "  stroke",
      "} def",

      # text is split into individual characters
      "/nullchar {} def",
      "/showchar {",
      "  exch dup 3 -1 roll",
      "  1 getinterval",
      "  convertToR /charfun (nullchar) put",
      "  convertToR /bboxfun (showbbox) put",
      "  mytext",
      # Generate charpath so that currentpoint returns correct location
      # (this also consumes the original char)
      "  true charpath", 
      "  currentpoint newpath moveto",
      "} def",
      "/xshowchar {",
      # Save copy of 'i'
      "  dup",
      # Get copy of 'str' out front
      "  2 index",
      "  exch",
      "  1 getinterval",
      # Save current location (starting position for THIS char)
      "  currentpoint 3 -1 roll",
      "  convertToR /charfun (nullchar) put",
      "  convertToR /bboxfun (showbbox) put",
      "  mytext",
      # Pop char
      "  pop",
      # Return to start point for THIS char and use 'numarray'
      # to find location for NEXT char
      "  moveto 2 index exch get 0 rmoveto",
      # Save current location (starting position for next char),
      # start new path for next char,
      # move to save location
      "  currentpoint newpath moveto",
      "} def",
      "/yshowchar {",
      # Save copy of 'i'
      "  dup",
      # Get copy of 'str' out front
      "  2 index",
      "  exch",
      "  1 getinterval",
      # Save current location (starting position for THIS char)
      "  currentpoint 3 -1 roll",
      "  convertToR /charfun (nullchar) put",
      "  convertToR /bboxfun (showbbox) put",
      "  mytext",
      # Pop char
      "  pop",
      # Return to start point for THIS char and use 'numarray'
      # to find location for NEXT char
      "  moveto 2 index exch get 0 exch rmoveto",
      # Save current location (starting position for next char),
      # start new path for next char,
      # move to save location
      "  currentpoint newpath moveto",
      "} def",
      "/xyshowchar {",
      # Save copy of 'i'
      "  dup",
      # Get copy of 'str' out front
      "  2 index",
      "  exch",
      "  1 getinterval",
      # Save current location (starting position for THIS char)
      "  currentpoint 3 -1 roll",
      "  convertToR /charfun (nullchar) put",
      "  convertToR /bboxfun (showbbox) put",
      "  mytext",
      # Pop char
      "  pop",
      # Return to start point for THIS char and use 'numarray'
      # to find location for NEXT char
      "  moveto 2 index exch 2 mul 2 getinterval dup 0 get exch 1 get rmoveto",
      # Save current location (starting position for next char),
      # start new path for next char,
      # move to save location
      "  currentpoint newpath moveto",
      "} def",
      "/widthshowchar {",
      "  exch dup 3 -1 roll",
      "  1 getinterval",
      "  convertToR /charfun (nullchar) put",
      # Just get bbox of char 
      "  convertToR /bboxfun (showbbox) put",
      "  mytext",
      "  dup 0 get 3 index eq {4 index 4 index rmoveto} if",
      # Generate charpath so that currentpoint returns correct location
      # (this also consumes the original char)
      "  true charpath", 
      "  currentpoint newpath moveto",
      "} def",
      "/ashowchar {",
      "  exch dup 3 -1 roll",
      "  1 getinterval",
      "  convertToR /charfun (nullchar) put",
      # Just get bbox of char 
      "  convertToR /bboxfun (showbbox) put",
      "  mytext",
      # Generate charpath so that currentpoint returns correct location
      # (this also consumes the original char)
      "  true charpath", 
      # Do adjustment on every char
      "  2 index 2 index rmoveto",
      "  currentpoint newpath moveto",
      "} def",
      "/awidthshowchar {",
      "  exch dup 3 -1 roll",
      "  1 getinterval",
      "  convertToR /charfun (nullchar) put",
      # Just get bbox of char 
      "  convertToR /bboxfun (showbbox) put",
      "  mytext",
      # Generate charpath so that currentpoint returns correct location
      # (this also consumes the original char)
      # Check for special char adjustment
      "  dup 0 get 5 index eq {6 index 6 index rmoveto} if",
      "  true charpath", 
      "  2 index 2 index rmoveto",
      "  currentpoint newpath moveto",
      "} def",
      
      # text is split into individual characters,
      # each character is converted to a path, flattened
      # and then stroked
      "/strokechar {",
      "  exch dup 3 -1 roll",
      "  1 getinterval",
      "  dup", # copy of char for "char="
      "  true charpath flattenpath",
      "  {mymove} {myline} {mycurve} {myclose}",
      "  mychar",
      "  pop", # copy of char for "char="
      # Save current location (starting position for next char),
      # start new path for next char,
      # move to save location
      "  currentpoint newpath moveto",
      "} def",
      "/xstrokechar {",
      # Save copy of 'i'
      "  dup",
      # Get copy of 'str' out front
      "  2 index",
      "  exch",
      "  1 getinterval",
      # Save current location (starting position for THIS char)
      "  currentpoint 3 -1 roll",
      "  dup", # copy of char for "char="
      "  true charpath flattenpath",
      "  {mymove} {myline} {mycurve} {myclose}",
      "  mychar",
      "  pop", # copy of char for "char="
      # Return to start point for THIS char and use 'numarray'
      # to find location for NEXT char
      "  moveto 2 index exch get 0 rmoveto",
      # Save current location (starting position for next char),
      # start new path for next char,
      # move to save location
      "  currentpoint newpath moveto",
      "} def",
      "/ystrokechar {",
      # Save copy of 'i'
      "  dup",
      # Get copy of 'str' out front
      "  2 index",
      "  exch",
      "  1 getinterval",
      # Save current location (starting position for THIS char)
      "  currentpoint 3 -1 roll",
      "  dup", # copy of char for "char="
      "  true charpath flattenpath",
      "  {mymove} {myline} {mycurve} {myclose}",
      "  mychar",
      "  pop", # copy of char for "char="
      # Return to start point for THIS char and use 'numarray'
      # to find location for NEXT char
      "  moveto 2 index exch get 0 exch rmoveto",
      # Save current location (starting position for next char),
      # start new path for next char,
      # move to save location
      "  currentpoint newpath moveto",
      "} def",
      "/xystrokechar {",
      # Save copy of 'i'
      "  dup",
      # Get copy of 'str' out front
      "  2 index",
      "  exch",
      "  1 getinterval",
      # Save current location (starting position for THIS char)
      "  currentpoint 3 -1 roll",
      "  dup", # copy of char for "char="
      "  true charpath flattenpath",
      "  {mymove} {myline} {mycurve} {myclose}",
      "  mychar",
      "  pop", # copy of char for "char="
      # Return to start point for THIS char and use 'numarray'
      # to find location for NEXT char
      "  moveto 2 index exch 2 mul 2 getinterval dup 0 get exch 1 get rmoveto",
      # Save current location (starting position for next char),
      # start new path for next char,
      # move to save location
      "  currentpoint newpath moveto",
      "} def",
      "/widthstrokechar {",
      "  exch dup 3 -1 roll",
      "  1 getinterval",
      "  dup",
      "  true charpath flattenpath",
      "  {mymove} {myline} {mycurve} {myclose}",
      "  mychar",
      # Check for special char adjustment
      # (the get converts char to int)
      "  0 get 2 index eq {3 index 3 index rmoveto} if",
      "  currentpoint newpath moveto",
      "} def",
      "/astrokechar {",
      "  exch dup 3 -1 roll",
      "  1 getinterval",
      "  dup", # copy of char for "char="
      "  true charpath flattenpath",
      "  {mymove} {myline} {mycurve} {myclose}",
      "  mychar",
      "  pop", # copy of char for "char="
      "  2 index 2 index rmoveto",
      "  currentpoint newpath moveto",
      "} def",
      "/awidthstrokechar {",
      "  exch dup 3 -1 roll",
      "  1 getinterval",
      "  dup",
      "  true charpath flattenpath",
      "  {mymove} {myline} {mycurve} {myclose}",
      "  mychar",
      # Do adjustment on every char
      "  3 index 3 index rmoveto",
      # Check for special char adjustment
      "  0 get 4 index eq {5 index 5 index rmoveto} if",
      "  currentpoint newpath moveto",
      "} def",

      "/showchars {",
      "  dup length -1 add 0 exch 1 exch {showchar} for",
      "} def",
      "/xshowchars {",
      "  dup length -1 add 0 exch 1 exch {xshowchar} for",
      "} def",
      "/yshowchars {",
      "  dup length -1 add 0 exch 1 exch {yshowchar} for",
      "} def",
      "/xyshowchars {",
      "  dup length -1 add 0 exch 1 exch {xyshowchar} for",
      "} def",
      "/widthshowchars {",
      "  dup length -1 add 0 exch 1 exch {widthshowchar} for",
      "} def",
      "/ashowchars {",
      "  dup length -1 add 0 exch 1 exch {ashowchar} for",
      "} def",
      "/awidthshowchars {",
      "  dup length -1 add 0 exch 1 exch {awidthshowchar} for",
      "} def",

      "/showpaths {",
      "  dup length -1 add 0 exch 1 exch {strokechar} for",
      "} def",
      "/xshowpaths {",
      "  dup length -1 add 0 exch 1 exch {xstrokechar} for",
      "} def",
      "/yshowpaths {",
      "  dup length -1 add 0 exch 1 exch {ystrokechar} for",
      "} def",
      "/xyshowpaths {",
      "  dup length -1 add 0 exch 1 exch {xystrokechar} for",
      "} def",
      "/widthshowpaths {",
      "  dup length -1 add 0 exch 1 exch {widthstrokechar} for",
      "} def",
      "/ashowpaths {",
      "  dup length -1 add 0 exch 1 exch {astrokechar} for",
      "} def",
      "/awidthshowpaths {",
      "  dup length -1 add 0 exch 1 exch {awidthstrokechar} for",
      "} def",

      "/showbbox {",
      "  gsave",
      "  currentpoint newpath moveto dup true charpath flattenpath",
      "  boxinit",
      "  {boxmove} {boxline} {boxcurve} {boxclose} pathforall",
      # Update global picture xmin/xmax/ymin/ymax
      "  boxupdate",
      "  grestore",
      "} def",
      "/xshowbbox {",
      "  gsave",
      "  currentpoint newpath moveto",
      # For each char, save currentpoint, flatten charpath,
      # move back to start point and shift across using 'numarray'
      "  dup length -1 add 0 exch 1 exch { dup 2 index exch 1 getinterval currentpoint 3 -1 roll true charpath flattenpath moveto 2 index exch get 0 rmoveto } for",
      "  boxinit",
      "  {boxmove} {boxline} {boxcurve} {boxclose} pathforall",
      # Update global picture xmin/xmax/ymin/ymax
      "  boxupdate",
      "  grestore",
      "} def",
      "/yshowbbox {",
      "  gsave",
      "  currentpoint newpath moveto",
      # For each char, save currentpoint, flatten charpath,
      # move back to start point and shift across using 'numarray'
      "  dup length -1 add 0 exch 1 exch { dup 2 index exch 1 getinterval currentpoint 3 -1 roll true charpath flattenpath moveto 2 index exch get 0 exch rmoveto } for",
      "  boxinit",
      "  {boxmove} {boxline} {boxcurve} {boxclose} pathforall",
      # Update global picture xmin/xmax/ymin/ymax
      "  boxupdate",
      "  grestore",
      "} def",
      "/xyshowbbox {",
      "  gsave",
      "  currentpoint newpath moveto",
      # For each char, save currentpoint, flatten charpath,
      # move back to start point and shift across using 'numarray'
      "  dup length -1 add 0 exch 1 exch { dup 2 index exch 1 getinterval currentpoint 3 -1 roll true charpath flattenpath moveto 2 index exch 2 mul 2 getinterval dup 0 get exch 1 get rmoveto } for",
      "  boxinit",
      "  {boxmove} {boxline} {boxcurve} {boxclose} pathforall",
      # Update global picture xmin/xmax/ymin/ymax
      "  boxupdate",
      "  grestore",
      "} def",
      "/widthshowbbox {",
      "  gsave",
      "  currentpoint newpath moveto",
      "  dup length -1 add 0 exch 1 exch { 1 index exch 1 getinterval dup true charpath flattenpath 0 get 2 index eq { 3 index 3 index rmoveto } if } for",
      "  boxinit",
      "  {boxmove} {boxline} {boxcurve} {boxclose} pathforall",
      # Update global picture xmin/xmax/ymin/ymax
      "  boxupdate",
      "  grestore",
      "} def",
      "/ashowbbox {",
      "  gsave",
      "  currentpoint newpath moveto",
      "  dup length -1 add 0 exch 1 exch { 1 index exch 1 getinterval true charpath flattenpath 2 index 2 index rmoveto } for",
      "  boxinit",
      "  {boxmove} {boxline} {boxcurve} {boxclose} pathforall",
      # Update global picture xmin/xmax/ymin/ymax
      "  boxupdate",
      "  grestore",
      "} def",
      "/awidthshowbbox {",
      "  gsave",
      "  currentpoint newpath moveto",
      "  dup length -1 add 0 exch 1 exch { 1 index exch 1 getinterval dup true charpath flattenpath 3 index 3 index rmoveto 0 get 4 index eq { 5 index 5 index rmoveto } if } for",
      "  boxinit",
      "  {boxmove} {boxline} {boxcurve} {boxclose} pathforall",
      # Update global picture xmin/xmax/ymin/ymax
      "  boxupdate",
      "  grestore",
      "} def",
      
      "/showwidthadj {",
      "  convertToR /curadj 0 put",
      "} def",
      "/widthshowwidthadj {",
      "  convertToR /curadj 0 put",
      # Determine width adj in device coords
      "  3 index 0 transwidth",
      # Check each char, if special, add 
      "  1 index { 3 index eq { convertToR /curadj get 1 index add convertToR /curadj 3 -1 roll put } if } forall",
      # Clean up
      "  pop",
      "} def",
      "/ashowwidthadj {",
      # How many chars
      "  dup length",
      # Determine width adj in device coords
      "  3 index 0 transwidth",
      "  mul convertToR /curadj 3 -1 roll put",
      "} def",
      "/awidthshowwidthadj {",
      # How many chars
      "  dup length",
      # Determine width adj in device coords
      "  3 index 0 transwidth",
      "  mul convertToR /curadj 3 -1 roll put",
      # Determine width adj in device coords
      "  5 index 0 transwidth",
      # Check each char, if special, add 
      "  1 index { 5 index eq { convertToR /curadj get 1 index add convertToR /curadj 3 -1 roll put } if } forall",
      # Clean up
      "  pop",
      "} def",
      
      "/show {",
      "  convertToR /charpathfun (showpaths) put",
      "  convertToR /charfun (showchars) put",
      "  convertToR /widthadjfun (showwidthadj) put",
      "  convertToR /bboxfun (showbbox) put",
      "  mytext",
      if (!(charpath || charpos)) {
          # Generate charpath so that currentpoint returns correct location
          # (this also consumes the original text)
          "  true charpath"
      } else {
          "pop"
      },
      "  currentpoint newpath moveto",
      "} def",
      "/xshow {",
      "  convertToR /charpathfun (xshowpaths) put",
      "  convertToR /charfun (xshowchars) put",
      "  convertToR /widthadjfun (showwidthadj) put",
      "  convertToR /bboxfun (xshowbbox) put",
      # Switch 'string' and 'numarray' so common code
      # for show operators can see 'string' on top of stack
      "  exch",
      "  mytext",
      if (!(charpath || charpos)) {
          # Perform the x-shifts from 'numarray'
          "  exch { 0 rmoveto } forall"
      } else {
          # remove 'numarray'
          "  exch pop"
      },
      # Remove 'string'
      "  pop", 
      "  currentpoint newpath moveto",
      "} def",
      "/yshow {",
      "  convertToR /charpathfun (yshowpaths) put",
      "  convertToR /charfun (yshowchars) put",
      "  convertToR /widthadjfun (showwidthadj) put",
      "  convertToR /bboxfun (yshowbbox) put",
      # Switch 'string' and 'numarray' so common code
      # for show operators can see 'string' on top of stack
      "  exch",
      "  mytext",
      if (!(charpath || charpos)) {
          # Perform the y-shifts from 'numarray'
          "  exch { 0 exch rmoveto } forall"
      } else {
          # remove 'numarray'
          "  exch pop"
      },
      # Remove 'string'
      "  pop", 
      "  currentpoint newpath moveto",
      "} def",
      "/xyshow {",
      "  convertToR /charpathfun (xyshowpaths) put",
      "  convertToR /charfun (xyshowchars) put",
      "  convertToR /widthadjfun (showwidthadj) put",
      "  convertToR /bboxfun (xyshowbbox) put",
      # Switch 'string' and 'numarray' so common code
      # for show operators can see 'string' on top of stack
      "  exch",
      "  mytext",
      if (!(charpath || charpos)) {
          # Perform the x/y-shifts from 'numarray' and pop 'numarray'
          "  exch dup length 2 idiv -1 add 0 exch 1 exch { 1 index exch 2 mul 2 getinterval dup 0 get exch 1 get rmoveto } for pop"
      } else {
          # remove 'numarray'
          "  exch pop"
      },
      # Remove 'string'
      "  pop", 
      "  currentpoint newpath moveto",
      "} def",
      "/widthshow {",
      "  convertToR /charpathfun (widthshowpaths) put",
      "  convertToR /charfun (widthshowchars) put",
      "  convertToR /widthadjfun (widthshowwidthadj) put",
      "  convertToR /bboxfun (widthshowbbox) put",
      "  mytext",
      if (!(charpath || charpos)) {
          # Check for any char adjustments
          c("  dup dup { 2 index eq {3 index 3 index rmoveto} if } forall",
          # Generate charpath so that currentpoint returns correct location
          # (this also consumes the original text)
            "  true charpath")
      } else {
          "pop"
      },
      # Remove the cx and cy and char from the /widthshow
      "  pop pop pop", 
      "  currentpoint newpath moveto",
      "} def",
      "/ashow {",
      "  convertToR /charpathfun (ashowpaths) put",
      "  convertToR /charfun (ashowchars) put",
      "  convertToR /widthadjfun (ashowwidthadj) put",
      "  convertToR /bboxfun (ashowbbox) put",
      "  mytext",
      if (!(charpath || charpos)) {
          # Do per-char adjustment
          c("  dup length dup 4 index mul exch 3 index mul rmoveto",
          # Generate charpath so that currentpoint returns correct location
          # (this also consumes the original text)
            "  true charpath")
      } else {
          "pop"
      },
      # Remove ax and ay
      "  pop pop", 
      "  currentpoint newpath moveto",
      "} def",
      "/awidthshow {",
      "  convertToR /charpathfun (awidthshowpaths) put",
      "  convertToR /charfun (awidthshowchars) put",
      "  convertToR /widthadjfun (awidthshowwidthadj) put",
      "  convertToR /bboxfun (awidthshowbbox) put",
      "  mytext",
      if (!(charpath || charpos)) {
          # Do per-char adjustment
          c("  dup length -1 add dup 4 index mul exch 3 index mul rmoveto",
          # Check for any char adjustments
            "  dup { 4 index eq {5 index 5 index rmoveto} if } forall",
          # Generate charpath so that currentpoint returns correct location
          # (this also consumes the original text)
            "  true charpath")
      } else {
          "pop"
      },
      # Remove ax,ay,char,cx,cy
      "  pop pop pop pop pop", 
      "  currentpoint newpath moveto",
      "} def",

      "end",
      # end global settings
      "false setglobal",
      "%%EndProcSet",
      "%% EndProlog",
      "",
      # Put converToR on the dictionary stack
      "convertToR begin",
      # Put userdict on top of dictionary stack for the file that is
      # about to run
      "userdict begin",
      if (!is.null(setflat)) {
          paste(setflat, "setflat")
      },
      "")
}

PScaptureFoot <-
    c(
      # Restore convertToR dict to dictionary stack in case the file
      # that we just ran does something nasty like clearing the dict stack!
      "convertToR begin",
      # XML file footer info
      "(<summary count=') print convertToR /id get 1 sub str cvs print (') print",
      "( ymax=') print convertToR /ymax get str cvs print (') print",
      "( ymin=') print convertToR /ymin get str cvs print (') print",
      "( xmax=') print convertToR /xmax get str cvs print (') print",
      "( xmin=') print convertToR /xmin get str cvs print (') print",
      "(/>\n\n) print",
      "(</picture>\n) print",

      # EOF
      "%% EOF\n"
      )

# Handle a "stringLine" that contains a NUL character
# ONE LINE ONLY at a time
processNUL <- function(filename, linenum) {
    if (length(linenum) != 1)
        stop("Bad computer!")
    # Read in each piece of the element
    bits <- scan(filename, skip=linenum - 1, nlines=1, quote=NULL,
                 sep=" ", quiet=TRUE, what="character")
    # The
    #   string='###TEXT blah blah ###TEXT'
    # part will be incomplete (terminated by the NUL char in blah blah)
    # Replace
    #   string='###TEXT whatever
    # with
    #   string=''
    bits[grep("string='###TEXT", bits)] <-
        paste(bits[grep("string='###TEXT", bits)], "###TEXT'", sep="")
    # OR could be
    #   char='###TEXT^@###TEXT'
    bits[grep("char='###TEXT", bits)] <- "char='###TEXT###TEXT'"
    # Check for ###TEXT' too
    bits[grep("###TEXT' $", bits)] <- " "
    # Put it all back together
    paste(bits, collapse=" ")
}

# Perform some post-processing of the XML file that is generated by the
# above PostScript code.
# (e.g., do some character escaping that would be more painful to do
#  in PostScript code)
postProcess <- function(outfilename, enc, defaultcol) {
    processStringLine <- function(stringLine) {
        paste(stringLine[1],
              gsub("<", "&lt;",
                   gsub(">", "&gt;",
                        gsub("'", "&apos;",
                             gsub("&", "&amp;",
                                  stringLine[2])))),
              stringLine[3],
              sep="")
    }
    # The XML file has been created by ghostscript in ISO-8859-1
    infile <- file(outfilename, "r", encoding="ISO-8859-1")
    lines <- readLines(infile)
    close(infile)
    # All string values have been marked 
    stringLines <- grep("###TEXT", lines)
    if (length(stringLines) > 0) {
        # Check first for NUL-terminated text
        # Test is no tag terminator, ">", at end of line
        nulTerminated <- !grepl(">$", lines[stringLines])
        if (any(nulTerminated)) {
            for (i in stringLines[nulTerminated]) {
                lines[i] <- processNUL(outfilename, i)
            }
        }
        stringLinesBits <- strsplit(lines[stringLines], "###TEXT")
        stringLinesMod <- sapply(stringLinesBits, processStringLine)
        lines[stringLines] <- stringLinesMod
    }
    # Write the file using 'enc' encoding
    outfile <- file(outfilename, "w", encoding=enc)
    writeLines(lines, outfile)
    close(outfile)
}

# Generate RGML file from PostScript file
PostScriptTrace <- function(file, outfilename,
                            charpath=TRUE, charpos=FALSE,
                            setflat=NULL, defaultcol="black",
                            encoding="ISO-8859-1") {
    # Create temporary PostScript file which loads
    # dictionary redefining stroke and fill operators
    # and then runs target PostScript file
    psfilename <- paste("capture", basename(file), sep="")
    psfile <- file(psfilename, "w")
    writeLines(PScaptureHead(file, charpath, charpos,
                             setflat, encoding), psfile)
    # Reconstitute file name here to handle Windows-style paths
    # in the file name
    writeLines(paste("(", file.path(dirname(file), basename(file)),
                     ") run", sep=""), psfile)
    writeLines(PScaptureFoot, psfile)
    close(psfile)

    if (missing(outfilename)) {
        outfilename <- paste(basename(file), ".xml", sep="")
    }
    
    # Run temp file using ghostscript
    # Determination of ghostscript command syntax derived from dev2bitmap()
    gsexe <- Sys.getenv("R_GSCMD")
    if (.Platform$OS.type == "windows") {
        if(!nzchar(gsexe)) gsexe <- Sys.getenv("GSC")
        if(is.null(gsexe) || !nzchar(gsexe)) {
            poss <- Sys.which(c("gswin64c.exe", "gswin32c.exe"))
            poss <- poss[nzchar(poss)]
            gsexe <- if(length(poss)) poss else "gswin32c.exe"
        } else if(grepl(" ", gsexe, fixed = TRUE))
            gsexe <- shortPathName(gsexe)
        outfile <- tempfile()
    } else {
        if (is.null(gsexe) || !nzchar(gsexe)) {
            gsexe <- "gs"
            rc <- system(paste(shQuote(gsexe), "-help > /dev/null"))
            if (rc != 0) 
                stop("sorry, 'gs' cannot be found")
        }
        outfile <- "/dev/null"
    }
    cmd <- paste(gsexe, 
                 " -q -dBATCH -dNOPAUSE -sDEVICE=ps2write -sOutputFile=",
                 outfile, " -sstdout=",
                 outfilename, " ",
                 psfilename, sep="")
    ret <- switch(.Platform$OS.type,
                  unix = system(cmd),
                  windows = system(cmd, invisible = TRUE))
    if(ret != 0) {
        stop(gettextf("status %d in running command '%s'", ret, cmd),
             domain = NA)
    } else {
        postProcess(outfilename, encoding, defaultcol)
    }
    invisible(cmd)
}

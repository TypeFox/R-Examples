#	*	*	*			Package "Rgnuplot"		*	*	*

pausetermPresent <- NA

"%s%" <- function(x,y) paste(x,y,sep='', collapse='') # string concatenation operator



GpplotMap<-function(mapvectfiles=NA, projection='PlateCarree', linetype='l', linestyle=1, plotTitle=NA, maprastfile=NA, maprastpalette=NA, AdditionalCode=NA,projectionInit=NA,returnCode=FALSE)
{
validProjections<-c('Winkeltripel','Robinson','Aitoff','NaturalEarth','WernersEquivalent','HammerWagner','EckertI','AlbersConical','SansonFlamsteed',
'Lambert','PlateCarree','EstereoAzimutal','Ortographic','Mercator')
AddCode<-c(rep('',4),'set xrange [-w: w]\nset yrange [-h-600: h]\nap=7\n',rep('',9))
selProj<-which(validProjections==projection)
if (length(selProj)==0 ) stop('Invalid projection')
Initvalues<-c('.2','15','.2','15','0,300','12','12','-20,0','0','0','0','0,0','0,0,90,90','0')
if (length(linetype)!=length(mapvectfiles)) linetype<-rep(linetype,length.out=length(mapvectfiles))
if (length(linestyle)!=length(mapvectfiles)) linestyle<-rep(linestyle,length.out=length(mapvectfiles))
if (is.na(AdditionalCode)) AdditionalCode<-''
linestyle<-gsub('^(\\d+)$',' ls \\1',linestyle)
#if (is.na(maprastfile) | is.na(maprastpalettes)) { maprastfiles<-'';maprastpalettes<-'' }
if (!returnCode) s<-'load "projections.gnu"\n' else s<-''
s<-s %s% 'p=' %s% projection %s% 'Init('
if (is.na(projectionInit)) s<-s %s% Initvalues[selProj] %s% ')\n' else s<-s %s% projectionInit %s% ')\n'
if (is.na(plotTitle)) s<-s  %s% 'set title p\n' else s<-s %s% 'set title "' %s% plotTitle %s% '"\n'
#if (nchar(AdditionalCode)>0) s<-s %s% AdditionalCode %s% '"\n'
s<-s %s% 'set size ratio -1' %s% '\n'

s<-s %s% AddCode[selProj]
if (!is.na(AdditionalCode)) s<-s %s% AdditionalCode %s% '\n'

if (is.na(maprastfile)){
s<-s  %s% 'plot '
for (n in 1:length(mapvectfiles))
{
s<-s %s% '"'  %s% mapvectfiles[n] %s% '" using (' %s% projection %s% 'YC($2,$1)):(' %s% projection %s% 'XC($2,$1)) notit w ' %s% linetype[n] %s% linestyle[n]
if (n<length(mapvectfiles)) s<-s  %s% ', \\\n'
}
}
if (!is.na(maprastfile)){
paletteRGB<- read.table(maprastpalette,  stringsAsFactors=FALSE)
NpaletteColors<-dim(paletteRGB)[1]-1 # number of palette colors starting from zero
s<-s  %s% 'set pm3d map
unset key;unset tics;unset border;unset colorbox\n'
s<-s  %s% 'set palette model RGB file "' %s% maprastpalette %s% '" u ($1/255):($2/255):($3/255)
set cbrange[0:' %s% NpaletteColors %s% ']
set pm3d corners2color c1
splot "'  %s% maprastfile %s% '" using (' %s% projection %s% 'YC($2,$1)):(' %s% projection %s% 'XC($2,$1)):3 w pm3d notit '
if (!is.na(mapvectfiles)) 
{
s<-s  %s% ', \\\n'
for (n in 1:length(mapvectfiles))
{
s<-s  %s% '"'  %s% mapvectfiles[n] %s% '" using (' %s% projection %s% 'YC($2,$1)):(' %s% projection %s% 'XC($2,$1)):(1) notit w ' %s% linetype[n] %s% linestyle[n]
if (n<length(mapvectfiles)) s<-s  %s% ', \\\n'
}
}}
if (!returnCode) Gprun(s,TRUE) else return(s %s% '\n')
#cat(s)
}


Gpmatrix2PNG<-function(matM,PNGfile)
{# saves a matrix with a mask of a map to a PNG file
#library('png', character.only=TRUE)
Mwidth<-dim(matM)[1]
Mheight<-dim(matM)[2]
matM<-matM[,Mheight:1]
matRGB <- array(0,c(Mheight,Mwidth,3))
matRGB[,,1] <- matrix(matM,Mheight,Mwidth,byrow=TRUE)
matRGB[,,2] <- matrix(matM,Mheight,Mwidth,byrow=TRUE)
matRGB[,,3] <- matrix(matM,Mheight,Mwidth,byrow=TRUE)
png::writePNG(matRGB, target = PNGfile)
}

GpmapPNG2lines<-function(PNGfile, landoutlinefile)
{# draw "squarish" coastlines given a gridded map from a PNG image
#based on code from Prof. Patrick J. Bartlein "The End of the Rainbow? Color Schemes for Improved Data Graphics" 
#library('png', character.only=TRUE)
p<-png::readPNG(PNGfile)
p<-p*255
numrows<-dim(p)[1]
numcols<-dim(p)[2]
p1<-t(p[,,1])
p2<-p1[,numrows:1]
pxl.x<-360/numcols
pxl.y<-180/numrows
l<-seq(-180,180,pxl.x)
l.e<-l[-1]
l.w<-l[-length(l)]
l<-seq(-90,90,pxl.y)
l.n<-l[-1]
l.s<-l[-length(l)]
# loop over lats and lons and get p2 outline segments
# latitudes
XY <- matrix(nrow=numcols*numrows*4,ncol=2)
n <- 0
for (k in 1:numrows) {
  for (j in 1:(numcols-1)) {  
    if(p2[j,k] != p2[j+1,k]) {
      n <- n+1
      XY[n*2,1] <- l.e[j]; XY[n*2,2] <- l.s[k]; XY[n*2+1,1] <- l.e[j]; XY[n*2+1,2] <- l.n[k]
    }
  }
}
# longitudes
for (j in 1:numcols) {
  for (k in 1:(numrows-1)) {  
    if(p2[j,k] != p2[j,k+1]) {
      n <- n+1
      XY[n*2,1] <- l.w[j]; XY[n*2,2] <- l.n[k]; XY[n*2+1,1] <- l.e[j]; XY[n*2+1,2] <- l.n[k]
    }
  }
}
p2.segs <- na.omit(XY)
# save segments for reuse
tmp<-tempfile()
write.table(p2.segs, tmp, row.names=FALSE, sep=" ",col.names=FALSE)
Gprun('!awk \'(FNR+1) % 2 == 0 {print ""} 1\' ' %s% tmp %s% ' > ' %s%  landoutlinefile ,FALSE)
}

GpimageCrop<-function(fileIN, fileOUT,x1,y1,x2,y2, filetype='PNG', terminal=NULL)
{#crop an image
#based on Qing Jie Li http://gnuplot-surprising.blogspot.fi/
if (is.null(terminal)) terminal<-Gpext2terminal(filetype)
#Initialize the gnuplot handle
h1<-Gpinit()
#change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
Gpcmd(h1,'reset
x1=' %s% x1  %s% '
y1=' %s% y1  %s% '
x2=' %s% x2  %s% '
y2=' %s% y2  %s% '
set xrange [x1:x2]
set yrange [y1:y2]
set size ratio -1
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
unset key
unset tics
unset border
set term ' %s% terminal  %s% ' size (x2-x1),(y2-y1)
set output "' %s% fileOUT  %s% '"
plot "' %s% fileIN  %s% '" binary filetype=auto w rgba' )
#close gnuplot handle
h1<-Gpclose(h1)  
}

GpimageRotate<-function(fileIN, fileOUT,degrees, filetype='PNG', terminal=NULL)
{#rotate an image
#based on Qing Jie Li http://gnuplot-surprising.blogspot.fi/
if (is.null(terminal)) terminal<-Gpext2terminal(filetype)
h1<-Gpinit()#Initialize the gnuplot handle
#change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
Gpcmd(h1,'reset
angle=' %s% degrees  %s% '*pi/180
set size ratio -1
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
unset key
unset tics
unset border
plot  "' %s% fileIN  %s% '"  binary filetype=auto rotate=angle w rgbima notitle
xmax=GPVAL_DATA_X_MAX
xmin=GPVAL_DATA_X_MIN
ymin=GPVAL_DATA_Y_MIN
ymax=GPVAL_DATA_Y_MAX
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set term ' %s% terminal  %s% ' size (xmax-xmin),(ymax-ymin)
set output "' %s% fileOUT  %s% '"
plot "' %s% fileIN  %s% '" binary filetype=auto rotate=angle w rgba notitle' )
#close gnuplot handle
h1<-Gpclose(h1)  
}

GpimageFlip<-function(fileIN, fileOUT,flipX=TRUE,flipY=FALSE, filetype='PNG', terminal=NULL)
{#flips a png image on the X, Y or XY axis
#based on Qing Jie Li http://gnuplot-surprising.blogspot.fi/
if (is.null(terminal)) terminal<-Gpext2terminal(filetype)
strOptions<-c('set yrange [] reverse;','set xrange [] reverse;','set xrange [] reverse;set yrange [] reverse;')
Gprun('plot "' %s% fileIN  %s% '" binary filetype=auto w rgbimage
xmax=GPVAL_DATA_X_MAX
xmin=GPVAL_DATA_X_MIN
ymin=GPVAL_DATA_Y_MIN
ymax=GPVAL_DATA_Y_MAX
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set size ratio -1
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
unset key
unset border
unset tics
set term ' %s% terminal  %s% ' size (xmax-xmin), (ymax-ymin)
set output "' %s% fileOUT  %s% '"
' %s% strOptions[flipX*1+flipY*2] %s% 'plot "' %s% fileIN  %s% '" binary filetype=auto w rgba notitle')
}

GpimageResize<-function(fileIN, fileOUT,newWidth, newHeight, filetype='PNG', terminal=NULL)
{#modify an image
#based on Qing Jie Li http://gnuplot-surprising.blogspot.fi/
if (is.null(terminal)) terminal<-Gpext2terminal(filetype)
h1<-Gpinit()#Initialize the gnuplot handle
#change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
Gpcmd(h1,'reset
plot "' %s% fileIN  %s% '" binary filetype=auto w rgbimage
xmax=GPVAL_DATA_X_MAX
xmin=GPVAL_DATA_X_MIN
ymin=GPVAL_DATA_Y_MIN
ymax=GPVAL_DATA_Y_MAX
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set size ratio -1
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
unset key
unset border
unset tics
set term ' %s% terminal  %s% ' size ' %s% newWidth  %s% ',' %s% newHeight  %s% '
set output "' %s% fileOUT  %s% '"
plot "' %s% fileIN  %s% '" binary filetype=auto w rgba' )
#close gnuplot handle
h1<-Gpclose(h1)  
}

GpimageDecimate<-function(fileIN, fileOUT, Xdec=2,Ydec=2, filetype='PNG', terminal=NULL)
{#Decimation of an image
if (is.null(terminal)) terminal<-Gpext2terminal(filetype)
Gprun('plot "' %s% fileIN  %s% '" binary filetype=auto w rgbimage
xmax=GPVAL_DATA_X_MAX
xmin=GPVAL_DATA_X_MIN
ymin=GPVAL_DATA_Y_MIN
ymax=GPVAL_DATA_Y_MAX
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set size ratio -1
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
unset key
unset border
unset tics
set term ' %s% terminal  %s% ' size (xmax-xmin), (ymax-ymin)
set output "' %s% fileOUT  %s% '"
plot "' %s% fileIN  %s% '" binary filetype=auto every ' %s% Xdec  %s% ':' %s% Ydec  %s% ' w rgba notitle')
}

Gpimage2PNG<-function(fileIN, fileOUT, optional256=FALSE,alpha=FALSE, backgroundColor='')
{# convert from any supported format to PNG, the file type is automatically detected
#optional256=TRUE converts a png image to 256 indexed colors, #however the palette is not optimized as in The Gimp
if (backgroundColor!='') backgroundColor <- '\nset object 1 rectangle from GPVAL_DATA_X_MIN,GPVAL_DATA_Y_MIN to GPVAL_DATA_X_MAX,GPVAL_DATA_Y_MAX fillcolor rgb "' %s% backgroundColor %s% '" back\n'
s<-'plot "' %s% fileIN  %s% '" binary filetype=auto w rgbimage
xmax=GPVAL_DATA_X_MAX
xmin=GPVAL_DATA_X_MIN
ymin=GPVAL_DATA_Y_MIN
ymax=GPVAL_DATA_Y_MAX
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set size ratio -1
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
unset key
unset border
unset tics' %s% backgroundColor %s% '\nset term png '
if (optional256==FALSE) s<-s %s% 'truecolor'%s% '  size (xmax-xmin), (ymax-ymin)
set output "' %s% fileOUT  %s% '"
plot "' %s% fileIN  %s% '" binary filetype=auto w '
if (alpha==FALSE) s<-s %s% 'rgbima notitle' else s<-s %s% 'rgba notitle'
Gprun(s)
}

Gpimage2image<-function(fileIN, fileOUT , filetype='PNG', terminal=NULL,alpha=FALSE, backgroundColor='')
{# convert from any supported format to other format, the file type is automatically detected
if (is.null(terminal)) terminal<-Gpext2terminal(filetype)
if (backgroundColor!='') backgroundColor <- '\nset object 1 rectangle from GPVAL_DATA_X_MIN,GPVAL_DATA_Y_MIN to   GPVAL_DATA_X_MAX,GPVAL_DATA_Y_MAX fillcolor rgb "' %s% backgroundColor %s% '" back\n'
s<-'plot "' %s% fileIN  %s% '" binary filetype=auto w rgbimage
xmax=GPVAL_DATA_X_MAX
xmin=GPVAL_DATA_X_MIN
ymin=GPVAL_DATA_Y_MIN
ymax=GPVAL_DATA_Y_MAX
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set size ratio -1
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
unset key
unset border
unset tics' %s% backgroundColor %s% '\nset term ' %s% terminal  %s% ' size (xmax-xmin), (ymax-ymin)
set output "' %s% fileOUT  %s% '"
plot "' %s% fileIN  %s% '" binary filetype=auto w '
if (alpha==FALSE) s<-s %s% 'rgbima notitle' else s<-s %s% 'rgba notitle'
Gprun(s)
}

#sepia algorithm based on http://stackoverflow.com/questions/5132015/how-to-convert-image-to-sepia-in-java
GpimageRgbfiltercolorSepia<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL,sepiaDepth=20,sepiaIntensity=10) GpimageRGBchange(fileIN, fileOUT, filetype, terminal,
'(p=(r+g+b)/3 ,r2 = p + (' %s% sepiaDepth %s% ' * 2), r2 = (r2>255)?255:r2, r2 )',
'(p=(r+g+b)/3 ,g2 = p + ' %s% sepiaDepth %s% ', g2 = (g2>255)?255:g2, g2 )',
'(p=(r+g+b)/3 ,b2 = p - '  %s% sepiaIntensity %s% ', b2 = (b2<0)?0:b2, b2 = (b2>255)?255:b2, b2 )')

GpimageRgbfiltercolorSepia2<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL)
GpimageRGBchange(fileIN, fileOUT,'PNG',NULL,'(p=r * 0.393 + g * 0.769 + b * 0.189,(p>255)?255:p)','(p=r * 0.349 + g * 0.686 + b * 0.168,(p>255)?255:p))','(p=r * 0.272 + g * 0.534 + b * 0.131,(p>255)?255:p))')

GpimageRgbfiltercolorRed<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL) GpimageRGBchange(fileIN, fileOUT, filetype, terminal,'r','0','0')
GpimageRgbfiltercolorGreen<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL) GpimageRGBchange(fileIN, fileOUT, filetype, terminal,'0','g','0')
GpimageRgbfiltercolorBlue<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL) GpimageRGBchange(fileIN, fileOUT, filetype, terminal,'0','0','b')
GpimageRgbfalsecolor<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL) GpimageRGBchange(fileIN, fileOUT, filetype, terminal,'g','b','r')
GpimageRgbgreyscaleY<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL) GpimageRGBchange(fileIN, fileOUT, filetype, terminal,'0.299*r + 0.587*g + 0.114*b') # YIQ/NTSC - RGB colors in a gamma 2.2 color space
GpimageRgbgreyscaleLinear<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL) GpimageRGBchange(fileIN, fileOUT, filetype, terminal,'0.3086*r + 0.6094*g + 0.0820*b') # linear RGB colors
GpimageRgbgreyscaleRMY<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL) GpimageRGBchange(fileIN, fileOUT, filetype, terminal,'0.5*r + 0.419*g + 0.081*b') # RMY
GpimageRgbgreyscaleBT709<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL) GpimageRGBchange(fileIN, fileOUT, filetype, terminal,'0.2125*r + 0.7154*g + 0.0721*b') # BT709
GpimageRgbgreyscaleavg<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL) GpimageRGBchange(fileIN, fileOUT, filetype, terminal,'(r + g + b)/3')
GpimageRgbgreyscaleLuminosity<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL) GpimageRGBchange(fileIN, fileOUT, filetype, terminal,'((R>G)?((R>B)?R:B):((G>B)?G:B) + (R<G)?((R<B)?R:B):((G<B)?G:B)) / 2')
#GpimageRgbBW<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL) GpimageRGBchange(fileIN, fileOUT, filetype, terminal,'(r==250 & g==250 & b==250)?255:0') # 2 colors BW



GpimageRGBchange<-function(fileIN, fileOUT, filetype='PNG', terminal=NULL,rgbformula,rgbformulaG='',rgbformulaB='')
{#modify an image with either a formula for each color component of RGB or for all
#based on Qing Jie Li http://gnuplot-surprising.blogspot.fi/
if (is.null(terminal)) terminal<-Gpext2terminal(filetype)
rflag <- ((rgbformulaG=='') & (rgbformulaB=='')) # rflag TRUE means that there is one formula for all RGB color components

h1<-Gpinit()#Initialize the gnuplot handle
#set gnuplot's additional search directories, to the extdata directory from Rgnuplot (default)
Gpsetloadpath(h1)
#change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
gcmd <- 'reset
plot "' %s% fileIN  %s% '" binary filetype=auto w rgbimage
xmax=GPVAL_DATA_X_MAX
xmin=GPVAL_DATA_X_MIN
ymin=GPVAL_DATA_Y_MIN
ymax=GPVAL_DATA_Y_MAX
set xrange [xmin:xmax]
set yrange [ymin:ymax]
rgbchange(r,g,b)=' %s% rgbformula  %s% '
set size ratio -1
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
unset key
unset border
unset tics
set term ' %s% terminal  %s% ' size (xmax-xmin),(ymax-ymin)
set output "' %s% fileOUT  %s% '"\n'
if (rflag) gcmd <- gcmd %s% 'plot "' %s% fileIN  %s% '" u(rgbchange(column(1),column(2),column(3))):(rgbchange(column(1),column(2),column(3))):(rgbchange(column(1),column(2),column(3))) binary filetype=auto w rgba' 
else gcmd <- gcmd %s% 'rgbchange2(r,g,b)=' %s% rgbformulaG  %s% '
rgbchange3(r,g,b)=' %s% rgbformulaB  %s% '
plot "' %s% fileIN  %s% '" u(rgbchange(column(1),column(2),column(3))):(rgbchange2(column(1),column(2),column(3))):(rgbchange3(column(1),column(2),column(3))) binary filetype=auto w rgba' 
#cat(gcmd)
Gpcmd(h1,gcmd)
#close gnuplot handle
h1<-Gpclose(h1)  
}

GpRGBsample<-function(xRGB='#FFFFFF', optionalTitle='')
{#plot a square filled with an RGB color
#Initialize the gnuplot handle
if (optionalTitle!='') optionalTitle<-'\nset title "' %s% optionalTitle %s% '"\n'
h1<-Gpinit()
#change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
Gpcmd(h1,'reset
set object  1 rect from screen 0, 0 to screen 1, 1 behind fc rgb "' %s% xRGB %s% '"  fillstyle  solid 1.00 noborder
unset xtics 
unset ytics 
unset key
unset border' %s% optionalTitle %s% '
plot 0 lc rgb "' %s% xRGB %s% '" notit')
Gppause()
#close gnuplot handle
h1<-Gpclose(h1)  
}

GpresampleDEM<-function(fileIN, fileOUT, xrange=NA, yrange=NA, interpolationMethod='50,50 gauss 1',XYformat=FALSE)#(xmax-xmin+1),(ymax-ymin+1) qnorm 10
{#resamples DEM data, similarly to what GIS applications do
# interpolationMethod = { splines | qnorm {<norm>} | (gauss | cauchy | exp | box | hann) {<dx>} {,dy} }
h1<-Gpinit()#Initialize the gnuplot handle
#change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
if (XYformat) u<-'u 1:2:3' else u<-'matrix'
gcmd <- 'reset
'
if (is.na(xrange) | is.na(yrange)) gcmd <- gcmd %s% 'splot "' %s% fileIN  %s% '" matrix w l
xmax=GPVAL_DATA_X_MAX
xmin=GPVAL_DATA_X_MIN
ymin=GPVAL_DATA_Y_MIN
ymax=GPVAL_DATA_Y_MAX' else gcmd <- gcmd %s% 'xmin=' %s% xrange[1]  %s% '
xmax=' %s% xrange[2]  %s% '
ymin=' %s% yrange[1]  %s% '
ymax=' %s% yrange[2]  %s% ''
gcmd <- gcmd %s% '
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set dgrid3d ' %s% interpolationMethod  %s% '
set table "' %s% fileOUT  %s% '"
splot "' %s% fileIN  %s% '" ' %s% u %s%'  w l
unset table'
Gpcmd(h1, gcmd)
#cat(gcmd)
#close gnuplot handle
h1<-Gpclose(h1)  
}

GphexRGB<-function(r,g=NA,b=NA)
{#returns an RGB sequence, from decimal to hexadecimal, ready for gnuplot
if (is.na(g) & is.na(b) & length(r)==3) sprintf("#%02X%02X%02X", as.integer(r[1]),as.integer(r[2]),as.integer(r[3]))
else sprintf("#%02X%02X%02X", as.integer(r),as.integer(g),as.integer(b))
}

GpimagePlot<-function(fileIN,alpha=FALSE, backgroundColor='', filetype='', terminal=NULL)
{#display an image file on the screen or to a chosen terminal except to a file
if (is.null(terminal))
{
terminal<-''
if (filetype!='') terminal<-'\nset term ' %s% Gpext2terminal(filetype) %s% '\n'
} 
rgbtype<-c('rgbima','rgbalpha')
if (backgroundColor!='') backgroundColor <- '\nset object 1 rectangle from GPVAL_DATA_X_MIN,GPVAL_DATA_Y_MIN to   GPVAL_DATA_X_MAX,GPVAL_DATA_Y_MAX fillcolor rgb "' %s% backgroundColor %s% '" back\nreplot'
h1<-Gpinit()#Initialize the gnuplot handle
#change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
Gpcmd(h1,terminal %s% '
set size ratio -1
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
unset key
unset tics
unset border
plot "' %s% fileIN  %s% '" binary filetype=auto w ' %s% rgbtype[alpha+1]  %s% ' notitle' %s% backgroundColor)
Gppause()
#close gnuplot handle
h1<-Gpclose(h1)  
}

Gprun <- function(cmd, optPause=FALSE){
 # execute a gnuplot string directly
if (!is.character(cmd)) stop('Argument <<cmd>> should be a string')
  h1<-Gpinit()
  Gpsetloadpath(h1)
  Gpsetwd(h1)
  Gpcmd(h1, cmd)
if (optPause) Gppause()
  h1<-Gpclose(h1)
}

GpRGB2image<-function(RGBfile, imagefile,width,height,filetype='PNG',gpterminal=NULL)
{#converts an RGB file to another image file
if (!is.character(RGBfile)) stop('Argument <<RGBfile>> should be a string')
if (!is.character(imagefile)) stop('Argument <<imagefile>> should be a string')
if (!is.numeric(width)) stop('Argument <<width>> should be integer')
if (!is.numeric(height)) stop('Argument <<height>> should be integer')
if (!file.exists(RGBfile)) stop('Error! File ' %s% RGBfile %s% ' does not exist' )
if (is.null(gpterminal)) gpterminal<-Gpext2terminal(filetype)
s<-'set terminal ' %s% gpterminal %s% ' size ' %s% width %s% ',' %s% height %s% ';set output "' %s% imagefile %s% '"
set size ratio -1
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
unset key
unset tics
unset border
plot "' %s% RGBfile %s% '" binary array=' %s% width %s% 'x' %s% height %s% ' flipy format="%uchar" with rgbimage'
Gprun(s)
#cat(s)
}

GpRGB2DAT<-function(RGBfile, DATfile,width,height)
{#converts an RGB file to a matrix data file
if (!is.character(RGBfile)) stop('Argument <<RGBfile>> should be a string')
if (!is.character(DATfile)) stop('Argument <<DATfile>> should be a string')
if (!is.numeric(width)) stop('Argument <<width>> should be integer')
if (!is.numeric(height)) stop('Argument <<height>> should be integer')
if (!file.exists(RGBfile)) stop('Error! File ' %s% RGBfile %s% ' does not exist' )
RGBsize <- file.info(RGBfile)$size
to.read = file(RGBfile, "rb")
RGBbin <- readBin(to.read, raw(), n = RGBsize, size = 1)
close(to.read)
matRGB <- matrix(RGBbin,width,height,byrow=TRUE)
Gpmatrixr2gnu(matRGB, DATfile)
}

Gpxydata2matrix <- function(fileIN,fileOUT){
# convert a file with splot data to a gnuplot matrix file

#tmpcntlines<-read.table(fileIN, header = F,sep='*',stringsAsFactors=FALSE,colClasses='character',comment.char = "+",blank.lines.skip=FALSE)
#tmpcntlines<-unlist(c(tmpcntlines))
#str(tmpcntlines)
#tcomments<-which(substr(tmpcntlines,1,1)=='#')
#tmpcntlines<-tmpcntlines[-tcomments]
#sum(tmpcntlines=='')

t1<-read.table(fileIN, header = F,strip.white=TRUE,blank.lines.skip=TRUE)#,sep=' '
#t1<-t1[,-c(1,3,5,7)]
t1<-t1[,1:3]
m<-length(unique(t1[,1]))
n<-length(t1[,1])/m
t2<-matrix(t1[  ,3],m,n,byrow=TRUE)
t2<-apply(t2,2,rev)  #flip the matrix
if (file.exists(fileOUT)) file.remove(fileOUT)
Gpmatrixr2gnu((t2),fileOUT)
}

GpmatrixfilePad<-function(fileIN,fileOUT,overwrite=TRUE)
{#add 1 row and 1 column (by duplicating the last row and last column) to a matrix file to display it 
#fileIN="3colors.dat";fileOUT="3colorsb.dat";overwrite=TRUE
if (!file.exists(fileIN)) stop('Error! File ' %s% fileIN %s% ' does not exist' )
if (overwrite==FALSE & file.exists(fileOUT)) stop('Error! File' %s% fileOUT %s% ' already exists' )
if (overwrite==TRUE & file.exists(fileOUT)) file.remove(fileOUT)
z<-read.delim(fileIN,sep=' ',header=FALSE)
numrows<-nrow(z)
numcols<-ncol(z)
z<-matrix(unlist(z),nrow=numrows)
z<-cbind(z,z[,numcols])
z2<-matrix(c(c(t(z)), z[numrows,]),numrows+1,numcols+1,byrow=TRUE)
write.table(z2, file=fileOUT, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

GpshowPaletteColornames<-function()
{#get gnuplot's RGB color names as a dataframe - 'ColorName','ColorHex','R','G','B'
logfile <- tempfile()
handle<-GpinitSaveStderr(logfile)
Gpcmd(handle, 'set print "' %s% logfile %s% '";show palette colornames;set print')
Sys.sleep(1) #wait 1 second
#s<-Gpfile2string(logfile)
handle<-Gpclose(handle)
s<-read.table(logfile,header=FALSE ,skip=1,strip.white=TRUE,comment.char='',stringsAsFactors=FALSE)
s<-s[,-3]
colnames(s)<-c('ColorName','ColorHex','R','G','B')
s
}

GpshowDatafileBinaryFiletypes<-function()
{#get gnuplot's RGB binary types as a vector
logfile <- tempfile()
handle<-GpinitSaveStderr(logfile)
Gpcmd(handle, 'set print "' %s% logfile %s% '";show datafile binary filetypes;set print')
Sys.sleep(1) #wait 1 second
handle<-Gpclose(handle)
s<-read.table(logfile,header=FALSE ,skip=2,strip.white=TRUE,comment.char='',stringsAsFactors=FALSE)
s<-c(unlist(s))
attr(s,'names')<-NULL
s
}

Gpmatrix2palette<-function(paletteData, paletteFileName,paletteIndeces=0,SolidColor=FALSE)
{#save a palette from matrix data with optional palette indeces and the option to modify the palette for solid colors
paletteSize<-length(paletteData)
paletteRows<-paletteSize/3
if (length(paletteIndeces)!=1) if (length(paletteIndeces)!=paletteRows) stop('Error! The number of indeces must be the same as the number of palette colors')
m<-matrix(paletteData,paletteRows,3,byrow=TRUE)
if (SolidColor==FALSE)
{
if (length(paletteIndeces)==1) 
{
Gpmatrixr2gnu(m,paletteFileName) # save palette data
} else {
m<-cbind(m, paletteIndeces)# add indices to the palette data
Gpmatrixr2gnu(m,paletteFileName)# save palette data
}
} else {# SolidColor==TRUE
m2<-matrix(0,paletteRows*2,ncol=3)
for (n in 1:paletteRows*2) { m2[n-1,]<-m[n/2,];m2[n,]<-m[n/2,] }
#m2<-rbind(m, m)
#m2<-matrix(c(m2),paletteRows*2,ncol=3,byrow=TRUE)
if (length(paletteIndeces)==1) 
{
paletteIndeces2<-c(-1,rep(0:(paletteRows-1),each=2,length=(paletteRows*2-1)))
m3<-cbind(m2, paletteIndeces2)# add indices to the palette data
Gpmatrixr2gnu(m3,paletteFileName)# save palette data
} else {
paletteIndeces2<- c(paletteIndeces[1]-1,rep(paletteIndeces,each=2,length=(paletteRows*2-1)))
m3<-cbind(m2, paletteIndeces2)# add indices to the palette data
Gpmatrixr2gnu(m3,paletteFileName)# save palette data
}
}
}

GpRGB1to3channels<-function(RGB1channel=NULL, fileRGB1channel=NULL, fileRGB3channel=NULL)
{#converts a vector or matrix of RGB 1-channel decimal values (24 bit) to 3-channel decimal values (3x8 bit)
#the input can be either a vector or matrix (RGB1channel) or a file (fileRGB1channel) 
#the output can be either a vector (returned value) or a file (fileRGB3channel)
outFile<-(is.null(fileRGB3channel))
if (!outFile) if (file.exists(fileRGB3channel)) file.remove(fileRGB3channel)
#if (outFile) fileRGB3channel
tmpOut<- tempfile()#output vector
if (is.null(fileRGB1channel))#input vector or matrix (RGB1channel)
{
if (is.matrix(RGB1channel)) RGB1channel<-c(matrix(RGB1channel,dim(RGB1channel)[1],dim(RGB1channel)[2],byrow=FALSE))
s<-'set table "' %s% tmpOut %s% '"
splot "-" u (int($1/65536)):((int($1) & 65280)/256):(int($1) & 255)
' %s% gsub('\n+$','',paste(RGB1channel,'\n',sep='',collapse='')) %s% '
e
unset table'
Gprun(s)
r<-read.table(tmpOut,skip=6,stringsAsFactors=FALSE, col.names=c('R','G','B','Z'))
r<-r[,-4]
if (outFile) return(r) #output vector
r<-as.matrix(r,ncol=3)
if (!outFile ) Gpmatrixr2gnu(r,fileRGB3channel) #write.table(r,file=fileRGB3channel, col.names=FALSE,row.names=FALSE)#output file
} else {#input file (fileRGB1channel)
#if (is.matrix(RGB1channel)) RGB1channel<-c(matrix(RGB1channel,dim(RGB1channel)[1],dim(RGB1channel)[2],byrow=FALSE))
Gprun('set table "' %s% tmpOut %s% '"
splot "' %s% fileRGB1channel %s% '" u (int($1/65536)):((int($1) & 65280)/256):(int($1) & 255)
unset table')
r<-read.table(tmpOut,skip=6,stringsAsFactors=FALSE, col.names=c('R','G','B','Z'))
r<-r[,-4]
if (outFile) return(r) #output vector
if (!outFile) write.table(r,file=fileRGB3channel, col.names=FALSE,row.names=FALSE)#output file
}
}

GpcreateIndexFromMatrixAndPalette<-function(matrixRGB, paletteRGB)
{#from a matrix with RGB colors (decimal 24bit) from an image file and its palette as a vector, create a matrix with indices
indexedVector<-matrix(0,dim(matrixRGB)[1],dim(matrixRGB)[2])
for (m in 1:dim(matrixRGB)[1])
for (n in 1:dim(matrixRGB)[2])
indexedVector[m,n]<-which(matrixRGB[m,n]==paletteRGB)-1
indexedVector
}

GpcreatePaletteFromMatrix<-function(matrixRGB, sortType='')
{#from a matrix with RGB colors (decimal 24bit) from an image file create a palette of 256 colors (decimal 24bit), as a vector
paletteRGB<-unique(c(matrixRGB))#get all the different colors
if (length(paletteRGB)>256) warning('The palette has ' %s% length(paletteRGB) %s% ' colors') #else {
#paletteRGB<-c(paletteRGB,rep(0,256-length(paletteRGB)))
#}
if (sortType=='') return(paletteRGB)
if (sortType=='desc')
{
#sort by descending number of occurrences, then by position on the image
uVectorSum<-sapply(paletteRGB, function(x) sum(x==matrixRGB))
uVectorPos<-sapply(paletteRGB, function(x) min(which(x==matrixRGB)))
uVector2<-cbind(paletteRGB, uVectorSum,uVectorPos)
uVector3<-uVector2[order(uVector2[,2], decreasing=TRUE),]
return(uVector3[,1])
}
if (sortType=='GIMP')
{
stmp<-sprintf('%x',paletteRGB)
uVectorSum<- as.numeric(paste('0x',sep='', substr(stmp,1,2))) + as.numeric(paste('0x',sep='', substr(stmp,3,4))) + as.numeric(paste('0x',sep='', substr(stmp,5,6))) 
uVector2<-cbind(paletteRGB, uVectorSum)
uVector3<-uVector2[order(uVector2[,2], decreasing=FALSE),]
return(uVector3[,1])
}
}

GppalettePlot<-function(filepal, sortType='', TheGimp=FALSE)
{#plots a palette from an indexed PNG file
PNGdata2<-GpPNG2color(filepal)#get the color matrix from an indexed PNG file
paletteRGB<-GpcreatePaletteFromMatrix(PNGdata2, sortType)#create a palette
tmppal<-tempfile()
GpRGB1to3channels(paletteRGB,fileRGB3channel=tmppal)#save the palette to a file with separated RGB components
if (!TheGimp) r<-read.table(tmppal, stringsAsFactors=FALSE) else {
r<-read.table(tmppal, skip=6, stringsAsFactors=FALSE)
r<-r[,-4]# get rid of "index"
}
r<-as.matrix(r)
tmppalsolid<-tempfile()
Gpmatrix2palette(c(t(r)), tmppalsolid,0,TRUE)#save the indexed color matrix to a file
Gprun('set view map
set size square
set yrange [] reverse
set cbrange [0:255]
set samples 16
set isosamples 16
set palette model RGB file "' %s% tmppalsolid %s% '" u ($1/255):($2/255):($3/255)
set pm3d corners2color c1
splot [1:16][1:16] "++" u 1:2:($1+(16*($2-1))) w pm notit',TRUE)
}

GpgimpPalette2matrix<-function(paletteGimp,returnIndex=FALSE)
{#reads a Gimp palette into a matrix, optionally the index can be returned too
r<-read.table(paletteGimp,skip=4,stringsAsFactors=FALSE, col.names=c('R','G','B','Z','index'))
r<-r[,-4]# get rid of "index"
if (!returnIndex) r<-r[,-4] # get rid of the index value
as.matrix(r)#data.frame to matrix
}

Gpmatrix2GimpPalette<-function(paletteMatrix, gplFileName, GimpColumns=16)
{#saves a matrix into a Gimp palette file (gpl)
Nrows <- dim(paletteMatrix)[1]
Ncols <- dim(paletteMatrix)[2]
if (Ncols==3) paletteMatrix <- cbind(paletteMatrix,0:(Nrows-1))
s <- 'GIMP Palette
Name: ' %s% gplFileName %s% '
Columns: ' %s% GimpColumns %s% '
#
'
paletteMatrix<-apply(paletteMatrix,1:2,as.character)
paletteMatrix<-apply(paletteMatrix,1:2,function(x) paste(substr('   ',1,3-nchar(x, type = "chars")),x,sep='',collapse=''))
paletteMatrix[,4]<-sub(' *([0-9]+)','	Index \\1\n',paletteMatrix[,4],perl=TRUE)
cat(file=gplFileName %s% '.gpl',s %s% paste(t(paletteMatrix),sep='',collapse=''))
}

Gpmapsr2gnu <- function(mapsObj,mapFileName)
{#saves the lines and points from a maps object to a data file readable by gnuplot
tmpMatrix<-cbind(mapsObj$x,mapsObj$y)
Gpmatrixr2gnu(tmpMatrix,mapFileName)
s<-Gpfile2string(mapFileName)
s<-gsub('NA NA', '', s, ignore.case = FALSE, perl=TRUE)
cat(s,file=mapFileName)
}

GpboxXY<-function(fileName)
{#returns the box around a location
m<-read.table(fileName, colClasses='numeric', stringsAsFactors=FALSE)
matrix(unlist(m),ncol=2,byrow=FALSE)
maxValues<-apply(m,2,max)
minValues<-apply(m,2,min)
r<-rbind(minValues, maxValues)
colnames(r)<-c('X','Y')
r
}

GpSHP2gnu<-function(SHPfilename, SHPlayername,gnufilename,toCRS='+init=epsg:4326')
{#given a shapefile (SHP) with full path and the shapefile layer name, the coordinates are saved to a text file readable by gnuplot
#library('rgdal', character.only=TRUE)
vecpxl = rgdal::readOGR(SHPfilename, SHPlayername)#load shapefile
fromCRS<-vecpxl@proj4string@projargs
isLines<-(is(vecpxl)[1]=='SpatialLinesDataFrame')
if (isLines) z<-vecpxl@lines else z<-vecpxl@polygons # get the lines or polygons
if (file.exists(gnufilename)) file.remove(gnufilename)
for (n in 1:length(z))
{
if (isLines) z1<-z[[n]]@Lines else z1<-z[[n]]@Polygons
for (m in 1:length(z1))
{
if (isLines) z2<-z[[n]]@Lines[[m]]@coords else z2<-z[[n]]@Polygons[[m]]@coords
if (toCRS==fromCRS) write(t(z2),gnufilename,ncolumns = 2,append=TRUE) else {
b_sp <- sp::SpatialPoints(z2, proj4string=sp::CRS(fromCRS))
z3<-sp::spTransform(b_sp, sp::CRS(toCRS)) # rgdal::spTransform
write(t(z3@coords),gnufilename,ncolumns = 2,append=TRUE)
}
cat(file=gnufilename,'\n',append=TRUE)
}
}
}

Gpext2terminal<-function(filetype='PNG')
{# determines a suitable terminal from a file extension
if (filetype=='JPEG') terminal<-'jpeg'
if (filetype=='PNG') terminal<-'pngcairo'
if (filetype=='SVG') terminal<-'svg'
if (filetype=='EPS') terminal<-'postscript eps enhanced color'
if (filetype=='PS') terminal<-'postscript'
terminal
}

GpimageTile<-function(fileOUT, matrixFilenamesIn, vectorWidths, vectorHeights, matrixXscale=NULL, matrixYscale=NULL,alpha=FALSE, filetype='PNG', terminal=NULL)
{#tiles multiple image files together into one
totalX<-sum(vectorWidths)
totalY<-sum(vectorHeights)
if (alpha) alphu <- ' u 1:2:3:4 w rgba' else alphu <- ' u 1:2:3 w rgbima'
if (is.null(terminal)) terminal<-Gpext2terminal(filetype)
s<-'set term ' %s% terminal %s% ' size ' %s% totalX %s% ',' %s% totalY %s% '
set output "' %s% fileOUT %s% '"
set xrange [0:' %s% (totalX-1) %s% '] noreverse nowriteback
set yrange [0:' %s% (totalY-1) %s% '] noreverse nowriteback
unset key;unset xtics; unset ytics; unset x2tics; unset y2tics;set format ""
unset border;unset title;unset label;unset legend
set lmargin 0
set rmargin 0
set bmargin 0
set tmargin 0
plot '
if (is.null(matrixXscale)) matrixXscale<-matrix(1,dim(matrixFilenamesIn)[1],dim(matrixFilenamesIn)[2],byrow=TRUE)
if (is.null(matrixYscale)) matrixYscale<-matrix(1,dim(matrixFilenamesIn)[1],dim(matrixFilenamesIn)[2],byrow=TRUE)
for (x in 1:dim(matrixFilenamesIn)[1])
for (y in 1:dim(matrixFilenamesIn)[2])
{
originX<-(totalX-sum(vectorWidths[1:x]))-1
if (originX<0) originX<-0
originY<-(totalY-sum(vectorHeights[2:y]))
#print(originY)
if (is.na(originY)) originY<-1
if (originY==0) originY<-1
s1<-'"' %s% matrixFilenamesIn[x,y] %s% '" binary filetype=auto origin=(' %s% originX %s% ',' %s% originY %s% ') dx=' %s% (matrixXscale[x,y]) %s% ' dy=' %s% (matrixYscale[x,y]) %s% alphu
if ((x != dim(matrixFilenamesIn)[1]) | (y != dim(matrixFilenamesIn)[2]))  s1<-s1 %s% ', '
s<-s %s% s1 
}
Gprun(s,FALSE)
#cat(s)
}

GpXYcoords2shpere<-function(xyfile,newxyfile)
{#read a XY data file and change the coords to fit a sphere
z<-read.table(xyfile, header = F,sep=' ',strip.white=TRUE,blank.lines.skip=TRUE)
z<-z[,c(1,3,5)]
Mwidth=length(unique(z[,2]))
Mheight=length(unique(z[,1]))
z[1]<-(-((z[1]-1)-Mheight/2)/Mheight*180)
z[2]<-(-((z[2]-1)-Mwidth/2)/Mwidth*360)
z2 <-data.matrix(z)
uniqz2<-unique(z2[,1])
if (file.exists(newxyfile)) file.remove(newxyfile)
for (uz in uniqz2) {write(t(z2[which(z2[,1]==uz),]),newxyfile,ncolumns = 3,append=TRUE);cat('\n',file=newxyfile,append=TRUE)}
}

GpsaveXYfile<-function(Rmatrix,xyfile)
{#saves an R matrix with coords as a XY data file
Mwidth=length(unique(Rmatrix[,2]))
Mheight=length(unique(Rmatrix[,1]))
z2 <-data.matrix(Rmatrix)
uniqz2<-unique(z2[,1])
if (file.exists(xyfile)) file.remove(xyfile)
for (uz in uniqz2) {write(t(z2[which(z2[,1]==uz),]),xyfile,ncolumns = 3,append=TRUE);cat('\n',file=xyfile,append=TRUE)}
}

Gpcolorhistogram<-function(fileIN, Ncols=1, vectCols=NULL, gtitle='',gxlabel='',gylabel='')
{#plots a color histogram
#based on http://gnuplot.sourceforge.net/demo_cvs/histograms.html
#example: file.copy('/usr/share/doc/gnuplot-doc/examples/immigration.dat',getwd() %s% '/immigration.dat')
#Gpcolorhistogram('immigration.dat', 4, c(6,12:14),'Immigration from Northern Europe\n(columstacked histogram)','Immigration by decade','Country of Origin')
if (is.null(vectCols)) vectCols<-1:Ncols
s<-'set border 3 front linetype -1 linewidth 1.000
set boxwidth 0.75 absolute
set style fill   solid 1.00 border lt -1
set grid nopolar
set grid noxtics nomxtics ytics nomytics noztics nomztics nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000
set key outside right top vertical Left reverse noenhanced autotitles columnhead box linetype -1 linewidth 1.000
set style histogram columnstacked title  offset character 0, 0, 0
set datafile missing "-"
set style data histograms
set xtics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set xtics norangelimit
set xtics ()
set ytics border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify
set ztics border in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
set cbtics border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify
set rtics axis in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
set title "' %s% gtitle %s% '" 
set xlabel "' %s% gxlabel %s% '" 
set ylabel "' %s% gylabel %s% '" 
set yrange [ 0.00000 : * ] noreverse nowriteback
i = 23
plot "' %s% fileIN %s% '" using ' %s% vectCols[1] %s% ' ti col'
for (n in 2:Ncols) s<-s %s% ', "" using ' %s% vectCols[n] %s% ' ti col'
Gprun(s,TRUE)
}

Gpmath3dPlot<-function(foo='(1+sgn(x-0.5))*(sin(4*pi*x) + sin(4*pi*y)+2)/2', xrange=c(0,1), yrange=c(0,1))
{# 3d plots like Mathematica
# from: Ben Ruijl Ben's blog http://negativeprobability.blogspot.fi/2012/07/making-gnuplot-plots-look-like.html
s<-'set palette defined ( 0 "#bd2c29", 2 "#ffd35a", 6 "white")
set pm3d implicit at s hidden3d 100
set style line 100 lc rgb "#000000" lt 1 lw 0.6
unset hidden3d
unset surface
set border 4095 front linetype -1 linewidth 0.8
set ticslevel 0
set xrange [' %s% as.character(xrange[1]) %s% ':' %s% as.character(xrange[2]) %s% ']
set yrange [' %s% as.character(yrange[1]) %s% ':' %s% as.character(yrange[2]) %s% ']
set ztics 1
set isosamples 40,40
set samples 20,20
unset colorbox
set view 32,28
splot ' %s% foo %s%  ' notit'
Gprun(s,TRUE)
}

GpmapMerpar<-function(parfile='worldpar.dat',merfile='worldmer.dat',pardeg=10,merdeg=10,splot=FALSE)
{# create a data file with meridians or parallels
if (splot) seprtr <- '\n\n' else seprtr <- '\n'
if (is.character(parfile))
{
if (file.exists(parfile)) file.remove(parfile)
for (m in seq(-90,90,pardeg))
{
for (n in seq(-180,180,merdeg))
cat(file=parfile,n,' ',m,seprtr,append=TRUE)
cat(file=parfile,'\n',append=TRUE)
}
}
if (is.character(merfile))
{
if (file.exists(merfile)) file.remove(merfile)
for (m in seq(-180,180,merdeg))
{
for (n in seq(-90,90,pardeg))
cat(file=merfile,m,' ',n,seprtr,append=TRUE)
cat(file=merfile,'\n',append=TRUE)
}
}
}

GpdivergingColormap<-function(s,rgb1,rgb2, outColorspace='sRGB')
{# This function is based on Kenneth Moreland's code for creating Diverging Colormaps.
# Matlab code created by Andy Stein. Translated to R by Jose Gama.
# s is a vector that goes between zero and one
# rgb1,rgb2 are objects from the colorspace package
# RGB, sRGB, HLS, HSV, LAB, LUV, PolarLAB, PolarLUV, XYZ
# outColorspace is the color space for the output
#library('colorspace', character.only=TRUE)
LabToMsh<-function(Lab)
{
L<-Lab@coords[1]
a<-Lab@coords[2]
b<-Lab@coords[3]
M <- sqrt(L*L + a*a + b*b)
s <- (M > 0.001) * acos(L/M)
h <- (s > 0.001) * atan2(b,a)
c(M,s,h)
}
MshToLab<-function(Msh)
{
M<-Msh[1]
s<-Msh[2]
h<-Msh[3]
L <- M*cos(s)
a <- M*sin(s)*cos(h)
b <- M*sin(s)*sin(h)
colorspace::LAB(L,a,b)
}

AngleDiff<-function(a1, a2)
{# Given two angular orientations, returns the smallest angle between the two.
v1<-matrix(c(cos(a1), sin(a1)),1,2,byrow=TRUE)
v2<-matrix(c(cos(a2), sin(a2)),1,2,byrow=TRUE)
x<-acos(sum(v1 * v2))
x
}
AdjustHue<-function(msh, unsatM)
{# For the case when interpolating from a saturated color to an unsaturated
# color, find a hue for the unsaturated color that makes sense.
        if (msh[1] >= unsatM-0.1  )                  
            # The best we can do is hold hue constant.
            h <- msh[3]
        else {
            # This equation is designed to make the perceptual change of the interpolation to be close to constant.
            hueSpin <- (msh[2]*sqrt(unsatM^2 - msh[1]^2)/(msh[1]*sin(msh[2])))
            # Spin hue away from 0 except in purple hues.
            if (msh[3] > -0.3*pi) h <- msh[3] + hueSpin else h <- msh[3] - hueSpin
            }
            h
}
divergingMap1val<-function(s, rgb1, rgb2, outColorspace='sRGB')
{
# Interpolate a diverging color map
# s is a number between 0 and 1
msh1 <- LabToMsh(as(rgb1, "LAB"))
msh2 <- LabToMsh(as(rgb2, "LAB"))
# If the endpoints are distinct saturated colors, then place white in between them
if (msh1[2] > 0.05 & msh2[2] > 0.05 & AngleDiff(msh1[3],msh2[3]) > pi/3)
{
# Insert the white midpoint by setting one end to white and adjusting the scalar value.
Mmid <- max(88.0, msh1[1], msh2[1])
#Mmid <- max(Mmid)
if (s < 0.5)
{
msh2[1] <- Mmid;  msh2[2] <- 0.0;  msh2[3] <- 0.0;s <- 2.0*s
} else {
msh1[1] <- Mmid;  msh1[2] <- 0.0;  msh1[3] <- 0.0; s <- 2.0*s - 1.0
}
}
# If one color has no saturation, then its hue value is invalid.  In this
# case, we want to set it to something logical so that the interpolation of hue makes sense.
        if ((msh1[2] < 0.05) & (msh2[2] > 0.05))
            msh1[3] <- AdjustHue(msh2, msh1[1]) else if ((msh2[2] < 0.05) & (msh1[2] > 0.05)) msh2[3] <- AdjustHue(msh1, msh2[1])
	mshTmp<-msh1
        mshTmp[1] <- (1-s)*msh1[1] + s*msh2[1]
        mshTmp[2] <- (1-s)*msh1[2] + s*msh2[2]
        mshTmp[3]<- (1-s)*msh1[3] + s*msh2[3]
        # Now convert back to the desired color space
 as(MshToLab(mshTmp),outColorspace)
}
dvmap<-matrix(0,length(s),3)
for (n in 1:length(s)) dvmap[n,]<-divergingMap1val(s[n], rgb1, rgb2, outColorspace)@coords
dvmap
}

GpXYcoordsConvertFun<-function(xyfile,newxyfile,fun1,fun2,swapXY=FALSE)
{#read a XY data file and change the coords according to 2 functions
z<-read.table(xyfile, header = F,sep=' ',strip.white=TRUE,blank.lines.skip=TRUE)
z<-z[,c(1,3,5)]
z[1]<-fun1(z[1])
z[2]<-fun2(z[2])
if (swapXY) z2 <-data.matrix(z[,c(2,1,3)]) else z2 <-data.matrix(z)
uniqz2<-unique(z2[,1])
if (file.exists(newxyfile)) file.remove(newxyfile)
for (uz in uniqz2) {write(t(z2[which(z2[,1]==uz),]),newxyfile,ncolumns = 3,append=TRUE);cat('\n',file=newxyfile,append=TRUE)}
}

GpplotFunction <- function(x, xlab='x', ylab='y', main='', type='l',...){
#coded by Oscar Perpinan Lamigueiro
if (!is.character(x)) stop('Argument <<x>> should be a string')
if (!is.character(xlab)) stop('Argument <<xlab>> should be a string')
if (!is.character(ylab)) stop('Argument <<ylab>> should be a string')
if (!is.character(main)) stop('Argument <<main>> should be a string')
if (!is.character(type)) stop('Argument <<type>> should be a string')
  h1<-Gpinit()
  GpsetXlabel(h1, xlab)
  GpsetYlabel(h1, ylab)
  style <- switch(type,
                  l = 'lines',
                  p = 'points'
                  )
  Gpsetstyle(h1, style)
  GpplotEquation(h1, as.character(x), main)
  Gppause()
  h1<-Gpclose(h1)
}

Gpsplot <- function(x, type = c('hidden3d', 'pm3d', 'map', 'contour')){
  #coded by Oscar Perpinan Lamigueiro
if (!is.numeric(x)) stop('Argument <<x>> should be numeric')
if (!is.numeric(x)) stop('Argument <<x>> should be numeric')

  h1<-Gpinit()

  old <- setwd(tempdir())

  Gpsetwd(h1)
  write.table(x, file='data.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)

  type <- match.arg(type)

  setCMD <- switch(type,
                   hidden3d = 'set hidden3d',
                   pm3d = 'set pm3d',
                   map = 'set pm3d map',
                   contour = paste('unset surface',
                     'set view map',
                     'set contour base', sep='\n')
                   )
  
  plotCMD <- switch(type,
                    hidden3d = 'splot "data.txt" matrix w l',
                    pm3d = 'splot "data.txt" matrix w pm3d',
                    map = 'splot "data.txt" matrix w pm3d',
                    contour = 'splot "data.txt" u 2:1:3 matrix w l notitle')

  cmd <- paste('reset', setCMD, 'unset key', plotCMD, sep = '\n')
  Gpcmd(h1, cmd)

  setwd(old)
}

GpplotPolyFit <- function(x, y, order){
#coded by Oscar Perpinan Lamigueiro
if (!is.numeric(x)) stop('Argument <<x>> should be numeric')
if (!is.numeric(y)) stop('Argument <<y>> should be numeric')
if (!is.numeric(order)) stop('Argument <<order>> should be numeric')

  old <- setwd(tempdir())
                                        #write the data to a text file
  write(t(cbind(x,y)),'data.txt',ncolumns =2)
                                        #initial values of the parameters
  initvalues <- rep(1, order+1)
  coefs <- letters[seq_len(order+1)]
  write(t(cbind(coefs,'=',initvalues)),'guess.txt',
        ncolumns = 3, sep='')
        
  h1<-Gpinit()
  Gpsetwd(h1, getwd())

  xs <- c('x', paste('x', 2:order, sep='**'))
  equation <- paste(coefs[1],
                    paste(coefs[-1], xs, sep='*', collapse='+'),
                    sep='+')
  equation <- paste('y(x)', equation, sep='=')
  cmd <- paste(equation,
               'fit y(x) "data.txt" via "guess.txt"',
               'set xlabel "x"',
               'set ylabel "y"',
               'unset key',
               'plot y(x), "data.txt"',
               sep='\n')

  Gpcmd(h1, cmd)
  setwd(old)
  }








Gpmandel <- function(x,y,z,n,maxiterations)
{#recursive implementation of the Mandelbrot set function

if (!is.numeric(x)) stop('Argument <<x>> should be numeric')
if (!is.numeric(y)) stop('Argument <<y>> should be numeric')
if (!is.numeric(z)) stop('Argument <<z>> should be numeric')
if (!is.numeric(n)) stop('Argument <<n>> should be numeric')
if (!is.numeric(maxiterations)) stop('Argument <<maxiterations>> should be numeric')

retmandel <- 0
ret <- .C("Rgnuplot_mandel",as.numeric(x),as.numeric(y),as.numeric(z),as.integer(n),as.integer(maxiterations),retmandel=as.integer(retmandel),DUP = TRUE,PACKAGE="Rgnuplot" )
ret$retmandel
}

Gpcols2rows <- function(filename,newfile,filecolseparator=' ',fileheader = FALSE,newfilerowseparator='\n\n')
{#convert a file with columns yyyy,m1,m2,m3,...,m12 to rows yyyy,m1 yyyy,m3 ... yyyy,m12 readable by gnuplot

if (!is.character(filename)) stop('Argument <<filename>> should be a string')
if (!is.character(newfile)) stop('Argument <<newfile>> should be a string')
if (!is.character(filecolseparator)) stop('Argument <<filecolseparator>> should be a string')
if (!is.logical(fileheader)) stop('Argument <<fileheader>> should be a boolean')
if (!is.character(newfilerowseparator)) stop('Argument <<newfilerowseparator>> should be a string')

# read table data, TAB separated
RTable<-read.table(filename, header = fileheader, sep = filecolseparator)
# store the number of rows and columns
iNrows<-dim(RTable)[1]
iNcols<-dim(RTable)[2]
# sort the repeated year values
yearvalues<-sort(rep(RTable[,1],iNcols-1))
# month values are needed for each year
month12<-rep(1:12,iNrows)
# convert to vector, by columns
Mdata<-c(t(as.matrix(RTable[-1])))
#create a matrix with the new data
newdata<-cbind(month12,yearvalues,Mdata)
newdata2<-newdata[order(month12),]
# save the new data
for (m in unique(newdata[,1]))
{
write.table(newdata2[which(newdata2[,1]==m),], file = newfile, sep = filecolseparator,row.names =FALSE, col.names =FALSE, quote =FALSE,append=TRUE)
z <- file(newfile, "at")
writeLines(newfilerowseparator,z, sep='')
close(z)
}
}

GpR2plot <- function(Rfunction,filename,gnuxrange,gnuyrange,gnusamples)
{#saves the output from a 2D R function to a file readable by gnuplot
if (!is.function(Rfunction)) stop('Argument <<Rfunction>> should be a function')
if (!is.character(filename)) stop('Argument <<filename>> should be a string')
if (!is.numeric(gnuxrange)) stop('Argument <<gnuxrange>> should be numeric')
if (!is.numeric(gnuyrange)) stop('Argument <<gnuyrange>> should be numeric')
if (!is.numeric(gnusamples)) stop('Argument <<gnusamples>> should be numeric')
col1 <- seq(from = gnuxrange[1], to = gnuxrange[2], length.out = gnusamples[1])
col2 <- Rfunction(col1)
col3 <- rep('i',times =length(col1))
col3[which((col2<gnuyrange[1]) | (col2>gnuyrange[2])) ] <- 'o'
rmatrix<-cbind(col1,col2,col3)
write(t(rmatrix),filename,ncolumns = 3)
}

GpR2splot <- function(Rfunction,filename,gnuxrange,gnuyrange,gnusamples,gnuisosamples,hidden3d=FALSE)
{#saves the output from a 3D R function to a file readable by gnuplot
if (!is.function(Rfunction)) stop('Argument <<Rfunction>> should be a function')
if (!is.character(filename)) stop('Argument <<filename>> should be a string')
if (!is.numeric(gnuxrange)) stop('Argument <<gnuxrange>> should be numeric')
if (!is.numeric(gnuyrange)) stop('Argument <<gnuyrange>> should be numeric')
if (!is.numeric(gnuisosamples)) stop('Argument <<gnuisosamples>> should be numeric')
if (!is.numeric(gnusamples)) stop('Argument <<gnusamples>> should be numeric')
if (length(gnuisosamples)==1) gnuisosamples <- c(gnuisosamples, gnuisosamples)
if (length(gnusamples)==1) gnusamples <- c(gnusamples, gnusamples)
if (hidden3d) filerows <- gnuisosamples[1] else filerows <- gnusamples[1]
col1 <- rep(seq(from = gnuxrange[1], to = gnuxrange[2], length.out = filerows),times=gnuisosamples[2])
col2 <- rep(seq(from = gnuyrange[1], to = gnuyrange[2], length.out = gnuisosamples[2]),each=filerows)
col3 <- vector(mode = "numeric", length =length(col1))
for (coloop in 1:length(col1)) col3[coloop] <- Rfunction(col1[coloop],col2[coloop])
col4 <- rep('i',times =length(col1))
col4[which((col2<gnuyrange[1]) | (col2>gnuyrange[2])) ] <- 'o'
rmatrix <- cbind(col1,col2,col3,col4)
n<-dim(rmatrix)[1]
for (eoloop in 1:(n / filerows-1)) rmatrix<-rbind(rmatrix[1:(eoloop*filerows+eoloop-1),],rep('',3),rmatrix[-1:-(eoloop*filerows+eoloop-1),])
write(t(rmatrix),filename,ncolumns = 4)
}

Gph <- function(handle, gnustring)
{#shows the output from gnuplot's help command
GpcheckHandle(handle)
if (!is.character(gnustring)) stop('Argument <<gnustring>> should be a string')
logfile <-tempfile()
Gpcmd(handle, 'set print "' %s% logfile %s% '";?' %s% gnustring %s% ';set print')
Gpfile2string(logfile)
}

GpPNG2RGB<-function(PNGfile, RGBfile,forceRGB=FALSE)
{# converts a PNG file to an RGB or RGBA file
#if forceRGB is TRUE then the alpha channel is ignored
#library('png', character.only=TRUE) #
if (!is.character(PNGfile)) stop('Argument <<PNGfile>> should be a string')
if (!is.character(RGBfile)) stop('Argument <<RGBfile>> should be a string')
if (!is.logical(forceRGB)) stop('Argument <<forceRGB>> should be logical')
if (!file.exists(PNGfile)) stop('Error! File ' %s% PNGfile %s% ' does not exist' )
p<-png::readPNG(PNGfile)
p<-p*255
numrows<-dim(p)[1]
numcols<-dim(p)[2]
#cat(dim(p)[3])
to.write = file(RGBfile, "wb")
if (forceRGB) rgbsize <- 3 else rgbsize <- 4
for (m in 1:numrows)
for (n in 1:numcols)
{
temp <- as.vector(as.integer(p[m,n,1:rgbsize]), mode = 'raw')
  writeBin(temp, to.write)
}
close( to.write )
}



Gpmatrixr2gnu <- function(rmatrix, gnufile)
{#saves an R matrix in a format that can be read by gnuplot
if (!is.matrix(rmatrix)) stop('Argument <<rmatrix>> should be a matrix')
if (!is.character(gnufile)) stop('Argument <<gnufile>> should be a string')
Ncol <- dim(rmatrix)[2]
write(t(rmatrix),gnufile,ncolumns = Ncol)
}

GpfitProgress <- function( fitprogressfile, fitparameter)
{#returns the value of a parameter from a logged stderr fit file
if (!is.character(fitprogressfile)) stop('Argument <<fitprogressfile>> should be a string')
if (!is.character(fitparameter)) stop('Argument <<fitparameter>> should be a string')
if (!file.exists(fitprogressfile)) stop('File ' %s% fitprogressfile %s% ' does not exist')
fitprogress <- Gpfile2string(fitprogressfile)
r <- regexpr('initial set of free parameter values[^\\/]+', fitprogress, perl=TRUE, useBytes = FALSE)
params <- substr(fitprogress,r[1],r[1]+attr(r,"match.length"))
params <- gsub('\\s*=.*?\n','\n', params)
params <- unlist(strsplit(params, '\n'))
params <- params[-c(1,2,length(params))]
if (!(fitparameter %in% c('WSSR','deltaWSSR','lambda',params))) stop('Parameter ' %s% fitparameter %s% ' does not exist')
if (fitparameter %in% params) grx <- gregexpr('\\s' %s% fitparameter %s% '\\s*\\=\\s*([0-9\\.\\-\\+e]+)',fitprogress, perl=TRUE) else {
if ( fitparameter == 'WSSR' ) grx <- gregexpr('\\sWSSR\\s*\\:\\s*([0-9\\.\\-\\+e]+)',fitprogress, perl=TRUE)
if ( fitparameter == 'deltaWSSR' ) grx <- gregexpr('delta\\(WSSR\\)\\s*\\:\\s*([0-9\\.\\-\\+e]+)',fitprogress, perl=TRUE)
if ( fitparameter == 'lambda' ) grx <- gregexpr('lambda\\s*\\:\\s*([0-9\\.\\-\\+e]+)',fitprogress, perl=TRUE)
}
x<-''
for (a in 1:length(grx[[1]])) x<-x %s%
substr(fitprogress,attr(grx[[1]],"capture.start")[a],attr(grx[[1]],"capture.start")[a]+attr(grx[[1]],"capture.length")[a])
x<-gsub('\\s+','\n',x)
z<-unlist(strsplit(x,'\n'))
z<-as.numeric(z)
r <- regexpr('After (\\d+) iterations the fit converged', fitprogress, perl=TRUE, useBytes = FALSE)
niterations <- substr(fitprogress,attr(r,"capture.start"),attr(r,"capture.start")+attr(r,"capture.length"))
niterations <- as.numeric(niterations)+1
z[1:niterations]
}

GpfitAllprogress <- function( fitprogressfile, boolScreenOut = TRUE)
{#returns the value of all parameters from a logged stderr fit file
if (!is.character(fitprogressfile)) stop('Argument <<fitprogressfile>> should be a string')
if (!is.logical(boolScreenOut)) stop('Argument <<boolScreenOut>> should be logical')
if (!file.exists(fitprogressfile)) stop('File ' %s% fitprogressfile %s% ' does not exist')
fitprogress <- Gpfile2string(fitprogressfile)
r <- regexpr('initial set of free parameter values[^\\/]+', fitprogress, perl=TRUE, useBytes = FALSE)
params <- substr(fitprogress,r[1],r[1]+attr(r,"match.length"))
params <- gsub('\\s*=.*?\n','\n', params)
params <- unlist(strsplit(params, '\n'))
params <- params[-c(1,2,length(params))]
r <- regexpr('After (\\d+) iterations the fit converged', fitprogress, perl=TRUE, useBytes = FALSE)
niterations <- substr(fitprogress,attr(r,"capture.start"),attr(r,"capture.start")+attr(r,"capture.length"))
niterations <- as.numeric(niterations)+1
matrixresult <- matrix( 0, nrow = niterations, ncol = 3 + length( params ) )
paramloop <- 1
for (fitparameter in c('WSSR','deltaWSSR','lambda',params))
{
if (fitparameter %in% params) grx <- gregexpr('\\s' %s% fitparameter %s% '\\s*\\=\\s*([0-9\\.\\-\\+e]+)',fitprogress, perl=TRUE) else {
if ( fitparameter == 'WSSR' ) grx <- gregexpr('\\sWSSR\\s*\\:\\s*([0-9\\.\\-\\+e]+)',fitprogress, perl=TRUE)
if ( fitparameter == 'deltaWSSR' ) grx <- gregexpr('delta\\(WSSR\\)\\s*\\:\\s*([0-9\\.\\-\\+e]+)',fitprogress, perl=TRUE)
if ( fitparameter == 'lambda' ) grx <- gregexpr('lambda\\s*\\:\\s*([0-9\\.\\-\\+e]+)',fitprogress, perl=TRUE)
}
x<-''
for (a in 1:length(grx[[1]])) x<-x %s%
substr(fitprogress,attr(grx[[1]],"capture.start")[a],attr(grx[[1]],"capture.start")[a]+attr(grx[[1]],"capture.length")[a])
x<-gsub('\\s+','\n',x)
z<-unlist(strsplit(x,'\n'))
z<-as.numeric(z)
z <- z[1:niterations]
matrixresult[,paramloop] <- z
paramloop <- paramloop + 1
}
if (boolScreenOut == TRUE)
{
titletext <- c('WSSR','deltaWSSR','lambda',params)
rettext <- paste(c('Iter. ',format(titletext, width=13,justify='right'),'\n'), sep = '', collapse = '') %s% 
paste(t(cbind(format(as.character(0:(niterations-1)),justify='left', width=max(6,nchar(niterations-1)+1)),format(matrixresult,scientific=TRUE, width=13), '\n')), sep = '', collapse = '')
cat(rettext)
} else return ( matrixresult )
}

GpinitSaveStderr <- function(logfile)
{#initialize gnuplot, save stderr to a log file
if (!is.character(logfile)) stop('Argument <<logfile>> should be a string')
Gpinit('gnuplot 2>' %s% logfile)
}

Gpinit <- function(optcmd='gnuplot')
{#initialize gnuplot
handle<-0#as.integer(0)
ret <- .C("Rgnuplot_init",rtrnvalue=handle,as.character(optcmd),DUP = TRUE,PACKAGE="Rgnuplot" )
ret$rtrnvalue
}

Gpgetfontpath<-function(handle)
{#get gnuplot's additional directories, for fonts
GpcheckHandle(handle)
options(warn=-1)
tmpfile<-tempfile()
Gpcmd(handle, 'save set "' %s% tmpfile %s% '"' )
Sys.sleep(1)
s<-Gpfile2string(tmpfile)
rx<-regexpr('\nset fontpath.*?\n',s)
s<-substr(s,rx[1],rx[1]+ attr(rx, "match.length")[1]-1)
s<-gsub("\\n", "", s, perl=TRUE)
s<-gsub("set fontpath", "", s, perl=TRUE)
s
}

Gpsetfontpath<-function(handle,fontpath=system.file(package='Rgnuplot') %s% '/extdata')
{#set gnuplot's additional directories, for fonts, default path = extdata directory from Rgnuplot
GpcheckHandle(handle)
options(warn=-1)
Gpcmd(handle, 'set fontpath "' %s% fontpath %s% '"')
}

Gpsetvariable<-function(handle,variablename,variabledata)
{#set a system or environment variable "variablename", with value "variabledata"
GpcheckHandle(handle)
options(warn=-1)
if (is.numeric(variabledata)) Gpcmd(handle, 'set ' %s% variablename %s% ' ' %s% variabledata ) else  Gpcmd(handle, 'set ' %s% variablename %s% ' "' %s% variabledata %s% '"')
}

Gpgetvariable<-function(handle,variablename)
{# returns the value of a system or environment variable "variablename"
GpcheckHandle(handle)
options(warn=-1)
#redirect output to a temporary file and save the variable's value
tmpfile<-tempfile()
Gpcmd(handle, 'set print "' %s% tmpfile %s% '";print ' %s% variablename %s% ';set print' )
Sys.sleep(1)
s<-Gpfile2string(tmpfile)
s # return the variable's value
}

Gpgetloadpath<-function(handle)
{#get gnuplot's additional directories, for data and scripts
GpcheckHandle(handle)
options(warn=-1)
tmpfile<-tempfile()
Gpcmd(handle, 'save set "' %s% tmpfile %s% '"' )
Sys.sleep(1)
s<-Gpfile2string(tmpfile)
rx<-regexpr('\nset loadpath.*?\n',s)
s<-substr(s,rx[1],rx[1]+ attr(rx, "match.length")[1]-1)
s<-gsub("\\n", "", s, perl=TRUE)
s<-gsub("set loadpath", "", s, perl=TRUE)
s
}

Gpsetloadpath<-function(handle,loadpath=system.file(package='Rgnuplot') %s% '/extdata')
{#set gnuplot's additional directories, for data and scripts, default path = extdata directory from Rgnuplot
GpcheckHandle(handle)
options(warn=-1)
Gpcmd(handle, 'set loadpath "' %s% loadpath %s% '"')
}

Gpversion<-function(handle)
{# returns the gnuplot version
GpcheckHandle(handle)
options(warn=-1)
#redirect output to a temporary file and save the working directory's path
tmpfile<-tempfile()
Gpcmd(handle, 'set print "' %s% tmpfile %s% '";print GPVAL_VERSION;set print' )
Sys.sleep(1)
s<-Gpfile2string(tmpfile)
s # return the gnuplot version
}

Gperrmsg<-function(handle)
{#get gnuplot's error messages
GpcheckHandle(handle)
options(warn=-1)
#redirect output to a temporary file and save the working directory's path
tmpfile<-tempfile()
Gpcmd(handle, 'set print "' %s% tmpfile %s% '";print GPVAL_ERRMSG;set print' )
Sys.sleep(1)
s<-Gpfile2string(tmpfile)
s # return the error message
}

Gpgetwd<-function(handle)
{#get gnuplot working directory
GpcheckHandle(handle)
options(warn=-1)
#redirect output to a temporary file and save the working directory's path
tmpfile<-tempfile()
Gpcmd(handle, 'set print "' %s% tmpfile %s% '";print GPVAL_PWD;set print' )
Sys.sleep(1)
s<-Gpfile2string(tmpfile)
s # return the working directory
}

Gpsetwd<-function(handle,wd=getwd())
{#set gnuplot working directory, default path = R's working directory
GpcheckHandle(handle)
options(warn=-1)
#print(wd)
Gpcmd(handle, 'cd \"' %s% wd %s% '\"')
}

GpPNG4DEM<-function(filePNG, fileDEMtab, file3Ddat)
{# converts a PNG file and a DEM tab separated to a text format readable by gnuplot 
p<-GpPNG2color(filePNG)
z<-read.delim(fileDEMtab,header=FALSE, stringsAsFactors =FALSE,comment.char = "#",sep=' ')#
numrows<-nrow(z)
numcols<-ncol(z)
z2<-matrix(unlist(z),nrow=numrows)
file.remove(file3Ddat)
for (m in 1:numrows)
{
for (n in 1:numcols) cat(file=file3Ddat, m, '\t', n, '\t', z[m,n], '\t', p[m,n], '\n',append=TRUE)
cat(file=file3Ddat,'\n',append=TRUE)
}
}

GpPNG2color<-function(fileName)
{# converts a PNG file to a text format readable by gnuplot 
#library('png', character.only=TRUE) #
if (!file.exists(fileName)) stop('Error! File does not exist' )
p<-png::readPNG(fileName)
p<-p*255
numrows<-dim(p)[1]
numcols<-dim(p)[2]
p1<-matrix(0,nrow=numrows, ncol=numcols)
for (m in 1:numrows)
for (n in 1:numcols)
#p1[m,n]<-p[numrows-m+1,numcols-n+1,1]*65536+p[numrows-m+1,numcols-n+1,2]*256+p[numrows-m+1,numcols-n+1,3]
p1[m,n]<-p[m,n,1]*65536+p[m,n,2]*256+p[m,n,3]
p1
}

Gpmatrix2XYdata<-function(fileName1,fileName2,optMatrix=NA,surfacegrid=FALSE,overwrite=TRUE)
{# converts a file with data in matrix format to an X, Y format readable by gnuplot
if (!file.exists(fileName1)) stop('Error! File does not exist' )
if (overwrite==FALSE & file.exists(fileName2)) stop('Error! File already exists' )
if (overwrite==TRUE & file.exists(fileName2)) file.remove(fileName2)
z<-read.table(fileName1,sep=' ')
numrows<-nrow(z)
numcols<-ncol(z)
z<-matrix(unlist(z),nrow=numrows)
#cat((dim(optMatrix)[1]),' ',(numrows*numcols),'\n')
if(class(optMatrix)=='matrix') if(dim(optMatrix)[1]!=numrows*numcols) stop('Error! The optional matrix has incorrect dimensions' )
#optMatrix<-matrix(unlist(optMatrix),nrow=numrows)
c=1
for (m in 1:numrows)
{
if(class(optMatrix)!='matrix') for (n in 1:numcols) cat(file=fileName2,m,'\t',n,'\t',z[m,n],'\n',append=TRUE)
else for (n in 1:numcols) 
{
cat(file=fileName2,m,'\t',n,'\t',z[m,n],'\t', optMatrix[c,1],'\n',append=TRUE)
c=c+1
}
if(surfacegrid) cat(file=fileName2,'\n\n',append=TRUE) else cat(file=fileName2,'\n',append=TRUE)
}
}

GppauseX<-function(delaySecs=2)
{# pauses the system for a number of seconds and then waits for the user to press a key, X11 only
if (!interactive()) return#return if not in interactive mode
print('Pause - close window to continue')
scounter<-0
while ((!GpisWindowOpen(0)) & (scounter<50)) { Sys.sleep(0.1);scounter<-scounter+1 } # wait for the graphics window to show up, for very slow calculations
repeat
{
Sys.sleep(delaySecs)
if(!GpisWindowOpen(0)) break
}
}

Gppause<-function(delaySecs=2)
{# pauses the system for a number of seconds and then waits for the user to press a key, detecs X11 beforehand
if (!interactive()) return#return if not in interactive mode
if (is.na(pausetermPresent)) pausetermPresent<-Gppausableterm()
if (pausetermPresent) GppauseX(delaySecs)
else {
Sys.sleep(delaySecs)
readline(prompt = "Pause. Press <Enter> to continue...")
}
}

Gpshowterm<-function()
{#get gnuplot's current terminal
logfile <- tempfile()
handle<-GpinitSaveStderr(logfile)
Gpcmd(handle, 'set print "' %s% logfile %s% '";show term;set print')
Sys.sleep(1) #wait 1 second
handle<-Gpclose(handle)
s<-Gpfile2string(logfile)
s
}

Gppausableterm<-function()
{#returns TRUE if the current terminal can be paused
if (is.na(pausetermPresent))
{
options(warn=-1)
curTerm <- Gpshowterm()
curTerm <- sub('\\s+terminal type is ', '', curTerm)
curTerm <- sub('\\s+\\d+\n$', '', curTerm)
if (curTerm %in% c('wxt','x11','qt','Windows')) pausetermPresent <-TRUE  else pausetermPresent <- FALSE
}
pausetermPresent
}

GpsetTerm<-function(optionalNewTerminal='')
{#set gnuplot's current terminal or get a list of all available terminals
logfile <- tempfile()
handle<-GpinitSaveStderr(logfile)
Gpcmd(handle, 'set print "' %s% logfile %s% '";set term ' %s% optionalNewTerminal %s% ';set print')
Sys.sleep(1) #wait 1 second
handle<-Gpclose(handle)
s<-Gpfile2string(logfile)
s
}

GpisWindowOpen<-function(windowid)
{# returns true if an X11 window is opened
#based on Joel VanderWerf http://www.ruby-forum.com/topic/144886
options(warn=-1)
rtxt<-system('xprop -name "Gnuplot (window id : "' %s% windowid %s% '")" WM_STATE 1', intern =TRUE)
if (length(rtxt)==0) return (FALSE)
return(grepl('Normal',rtxt[2]) | grepl('Iconic',rtxt[2]))#return true if the window is open
}

GpX11Present<-function()
{# returns true if X11 is present in the system
options(warn=-1)
rtxt<-system('echo $DISPLAY')
if (length(rtxt)==0)  return (TRUE) else return (FALSE)
}

GppidX11<-function(windowid=0)
{# returns the pid from an X11 window (gnuplot)
#based on Joel VanderWerf http://www.ruby-forum.com/topic/144886
options(warn=-1)
rtxt<-system('xprop -name "Gnuplot (window id : "' %s% windowid %s% '")" _NET_WM_PID', intern =TRUE)
if (length(rtxt)==0) return (NaN)
as.numeric(sub('\\D+','',rtxt))
}

GpkillpidX11<-function(windowid=0)
{# kills X11 windows (gnuplot)
#useful for testing purposes only- experimental, use at your own risk!
while (GpisWindowOpen(windowid))
{
pid <- GppidX11(windowid)
if (is.numeric(pid)) system('kill ' %s% pid, intern =TRUE)
Sys.sleep(1)
}
}

GploadDemo<-function(handle, mfile)
{# load a .dem gnuplot file and execute it, allowing pause statements
GpcheckHandle(handle)
m<-read.delim(mfile,sep='\n',stringsAsFactors=FALSE,header=FALSE,colClasses='character')
m<-unlist(m)
tmpfile<-tempfile()
Gpcmd(handle, 'set print "' %s% tmpfile %s% '"' )
for(lcode in m)
{
if (grepl('pause',lcode)) 
{
Gppause(2)
}
else {
Gpcmd(handle, lcode)
s<-Gpfile2string(tmpfile)
if (s !='') print(s)
}
}
}

GpWindowStatus<-function(windowid)
{# returns the status of an X11 window
#based on Joel VanderWerf http://www.ruby-forum.com/topic/144886
options(warn=-1)
rtxt<-system('xprop -name "Gnuplot (window id : "' %s% windowid %s% '")" WM_STATE 1', intern =TRUE)
if (length(rtxt)==0) return ('')
if (grepl('Withdrawn',rtxt[2])) return('Withdrawn')
if (grepl('Iconic',rtxt[2])) return('Iconic')
if (grepl('Normal',rtxt[2])) return('Normal')
}

Gpfile2string<-function(mfile)
{# read a text file to a string
if (!file.exists(mfile)) return('')
if (file.info(mfile)$size==0) return('')
#m<-read.delim(mfile,sep='\n',stringsAsFactors=FALSE,header=FALSE,colClasses='character')
m<-readLines(mfile)
m<-paste(m,sep='\n', collapse='\n')
#m<-unlist(m)
m
}

GpURL2string<-function(mURL)
{# read a text file from the web to a string
con <- url(mURL) # open a connection
mytxt <- readLines(con) # read the file
close(con) # close the connection
mytxt<-unlist(mytxt)
mytxt<-paste(mytxt,sep='\n', collapse='\n')
mytxt
}

GpcheckHandle <- function(handle)
{
#if (class(handle) != 'integer') stop('Argument <<handle>> is not a valid pointer')
if (handle==0) stop('Argument <<handle>> is not a valid pointer')
}

#	*	*	*	Function "Gpclose" Package "Rgnuplot"		*	*	*
Gpclose <- function(handle)
{
GpcheckHandle(handle)
ret <- .C("Rgnuplot_close",handle,DUP = TRUE,PACKAGE="Rgnuplot" )
as.integer(0)
}

#	*	*	*	Function "Gpcmd" Package "Rgnuplot"		*	*	*
Gpcmd <- function(handle,cmd, ...)
{
GpcheckHandle(handle)
if (!is.character(cmd)) stop('Argument <<cmd>> should be a string')
c<-list(...)
if(length(c)) ret <- .C("Rgnuplot_cmd",handle,as.character(cmd),...,DUP = TRUE,PACKAGE="Rgnuplot" ) else 
ret <- .C("Rgnuplot_send",handle,as.character(cmd),DUP = TRUE,PACKAGE="Rgnuplot" )
}

#	*	*	*	Function "Gpsetstyle" Package "Rgnuplot"		*	*	*
Gpsetstyle <- function(handle,plot.style)
{
GpcheckHandle(handle)
if (!is.character(plot.style)) stop('Argument <<plot.style>> should be a string')
ret <- .C("Rgnuplot_setstyle",handle,as.character(plot.style),DUP = TRUE,PACKAGE="Rgnuplot" )
}

#	*	*	*	Function "GpsetXlabel" Package "Rgnuplot"		*	*	*
GpsetXlabel <- function(handle,label)
{
GpcheckHandle(handle)
if (!is.character(label)) stop('Argument <<label>> should be a string')
ret <- .C("Rgnuplot_set_xlabel",handle,as.character(label),DUP = TRUE,PACKAGE="Rgnuplot" )
}

#	*	*	*	Function "GpsetYlabel" Package "Rgnuplot"		*	*	*
GpsetYlabel <- function(handle,label)
{
GpcheckHandle(handle)
if (!is.character(label)) stop('Argument <<label>> should be a string')
ret <- .C("Rgnuplot_set_ylabel",handle,as.character(label),DUP = TRUE,PACKAGE="Rgnuplot" )
}

#	*	*	*	Function "Gpresetplot" Package "Rgnuplot"		*	*	*
Gpresetplot <- function(handle)
{
GpcheckHandle(handle)
ret <- .C("Rgnuplot_resetplot",handle,DUP = TRUE,PACKAGE="Rgnuplot" )
}

#	*	*	*	Function "GpplotX" Package "Rgnuplot"		*	*	*
GpplotX <- function(handle,d,n,title)
{
GpcheckHandle(handle)
if (!is.numeric(d)) stop('Argument <<d>> should be a number')
if (!is.numeric(n)) stop('Argument <<n>> should be a number')
if (!is.character(title)) stop('Argument <<title>> should be a string')
ret <- .C("Rgnuplot_plot_x",handle,as.numeric(d),as.integer(n),as.character(title),DUP = TRUE,PACKAGE="Rgnuplot" )
}

#	*	*	*	Function "GpplotXY" Package "Rgnuplot"		*	*	*
GpplotXY <- function(handle,x,y,n,title)
{
GpcheckHandle(handle)
if (!is.numeric(x)) stop('Argument <<x>> should be a number')
if (!is.numeric(y)) stop('Argument <<y>> should be a number')
if (!is.numeric(n)) stop('Argument <<n>> should be a number')
if (!is.character(title)) stop('Argument <<title>> should be a string')
ret <- .C("Rgnuplot_plot_xy",handle,as.numeric(x),as.numeric(y),as.integer(n),as.character(title),DUP = TRUE,PACKAGE="Rgnuplot" )
}

#	*	*	*	Function "GpplotOnce" Package "Rgnuplot"		*	*	*
GpplotOnce <- function(title,style,label.x,label.y,x,y,n)
{
if (!is.character(title)) stop('Argument <<title>> should be a string')
if (!is.character(style)) stop('Argument <<style>> should be a string')
if (!is.character(label.x)) stop('Argument <<label.x>> should be a string')
if (!is.character(label.y)) stop('Argument <<label.y>> should be a string')
if (!is.numeric(x)) stop('Argument <<x>> should be a number')
if (!is.numeric(y)) stop('Argument <<y>> should be a number')
if (!is.numeric(n)) stop('Argument <<n>> should be a number')
ret <- .C("Rgnuplot_plot_once",as.character(title),as.character(style),as.character(label.x),as.character(label.y),as.numeric(x),as.numeric(y),as.integer(n),DUP = TRUE,PACKAGE="Rgnuplot" )
}

#	*	*	*	Function "GpplotSlope" Package "Rgnuplot"		*	*	*
GpplotSlope <- function(handle,a,b,title)
{
GpcheckHandle(handle)
if (!is.numeric(a)) stop('Argument <<a>> should be a number')
if (!is.numeric(b)) stop('Argument <<b>> should be a number')
if (!is.character(title)) stop('Argument <<title>> should be a string')
ret <- .C("Rgnuplot_plot_slope",handle,as.numeric(a),as.numeric(b),as.character(title),DUP = TRUE,PACKAGE="Rgnuplot" )
}

#	*	*	*	Function "GpplotEquation" Package "Rgnuplot"		*	*	*
GpplotEquation <- function(handle,equation,title)
{
GpcheckHandle(handle)
if (!is.character(equation)) stop('Argument <<equation>> should be a string')
if (!is.character(title)) stop('Argument <<title>> should be a string')
ret <- .C("Rgnuplot_plot_equation",handle,as.character(equation),as.character(title),DUP = TRUE,PACKAGE="Rgnuplot" )
}

#	*	*	*	Function "GpwriteXcsv" Package "Rgnuplot"		*	*	*
GpwriteXcsv <- function(fileName,x,n,title)
{
if (!is.character(fileName)) stop('Argument <<fileName>> should be a string')
if (!is.numeric(x)) stop('Argument <<x>> should be a number')
if (!is.numeric(n)) stop('Argument <<n>> should be a number')
if (!is.character(title)) stop('Argument <<title>> should be a string')
rtrnvalue<-0
ret <- .C("Rgnuplot_write_x_csv",as.character(fileName),as.numeric(x),as.integer(n),as.character(title),rtrnvalue=as.integer(rtrnvalue),DUP = TRUE,PACKAGE="Rgnuplot" )
ret$rtrnvalue
}

#	*	*	*	Function "GpwriteXYcsv" Package "Rgnuplot"		*	*	*
GpwriteXYcsv <- function(fileName,x,y,n,title)
{
if (!is.character(fileName)) stop('Argument <<fileName>> should be a string')
if (!is.numeric(x)) stop('Argument <<x>> should be a number')
if (!is.numeric(y)) stop('Argument <<y>> should be a number')
if (!is.numeric(n)) stop('Argument <<n>> should be a number')
if (!is.character(title)) stop('Argument <<title>> should be a string')
rtrnvalue<-0
ret <- .C("Rgnuplot_write_xy_csv",as.character(fileName),as.numeric(x),as.numeric(y),as.integer(n),as.character(title),rtrnvalue=as.integer(rtrnvalue),DUP = TRUE,PACKAGE="Rgnuplot" )
ret$rtrnvalue
}

#	*	*	*	Function "GpwriteMultiCsv" Package "Rgnuplot"		*	*	*
GpwriteMultiCsv <- function(fileName,xListPtr,n,numColumns,title)
{
if (!is.character(fileName)) stop('Argument <<fileName>> should be a string')
if (!is.numeric(xListPtr)) stop('Argument <<xListPtr>> should be a number')
if (!is.numeric(n)) stop('Argument <<n>> should be a number')
if (!is.numeric(numColumns)) stop('Argument <<numColumns>> should be a number')
if (!is.character(title)) stop('Argument <<title>> should be a string')
rtrnvalue<-0
ret <- .C("Rgnuplot_write_multi_csv",as.character(fileName),as.numeric(xListPtr),as.integer(n),as.integer(numColumns),as.character(title),rtrnvalue=as.integer(rtrnvalue),DUP = TRUE,PACKAGE="Rgnuplot" )
ret$rtrnvalue
}


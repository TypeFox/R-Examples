.onLoad <- function(lib, pkg)
{
  library.dynam('Rgnuplot', pkg, lib)
}

.onAttach <- function(lib, pkg)
{
checkG <- try(system('gnuplot -V', intern = TRUE, ignore.stderr = T),silent=TRUE)
checkG <- unlist(strsplit(checkG,' '))
notInstG <- TRUE
V <- ''
 if (length(checkG)>=2) if (checkG[1] == 'gnuplot') {
  V <- as.numeric(checkG[2])
  if (is.numeric(V)) if (V >= 4.6) notInstG <- FALSE
}

if (notInstG)
{
if (is.numeric(V)) if (V < 4.6)
{
packageStartupMessage('Rgnuplot works with gnuplot version 4.6 or higher\nThe current version of gnuplot in this system is less than 4.6 \nPlease reinstall it.\n')
} else packageStartupMessage('R cannot locate gnuplot, it is either not installed or there is no path to it.\nPlease install it.\n' )
mySystem <- .Platform$OS.type # unix windows
if (mySystem == 'unix') if (Sys.info()['sysname']=='Darwin') mySystem <- 'mac'
if (mySystem == 'unix')
{
packageStartupMessage('Installing gnuplot on Linux:

On Debian Linux derivatives:
sudo apt-get install gnuplot

On RedHat Linux derivatives:
yum install gnuplot*

How to install GNUPLOT+TikZ on Ubuntu:
Check GNUPLOT+TikZ For Stunning LaTeX Plots by Tim Teatro
http://www.timteatro.net/2010/07/18/gnuplottikz-for-stunning-latex-plots/

Building from CVS and more:
http://www.gnuplot.info/development/
')
} else if (mySystem == 'windows')
{
packageStartupMessage('Installing gnuplot on Windows:

Download gnuplot for Windows:
http://www.gnuplot.info/download.html

Choose the MinGW Windows binaries built by Tatsuro Matsuoka. Install it.

On Windows 7:

Control Panel/System and Security/System
System Properties/Environment Variables/Edit System Variable
Add to variable "Path" the path to gnuplot, by default:
C:\\Program Files\\gnuplot\\bin\\

Building from CVS and more:
http://www.gnuplot.info/development/
')
} else if (mySystem == 'mac')
{
packageStartupMessage('Installing gnuplot on MacOSX:

Download the gnuplot source:
http://sourceforge.net/projects/gnuplot/files/gnuplot/

Unzip to a folder created for gnuplot.

Edit Makefile (sub folder shlib in the ReadLine source folder), search for -dynamic and change it to -dynamiclib

Type in a console:
make install

Change the directory to the gnuplot source folder, type:
./configure --with-readline=< path of the gnuplot folder>
make install

Source of this info:
How to install gnuplot in Mac OS X lion
http://bhou.wordpress.com/2011/09/13/how-to-install-gnuplot-in-mac-os-x-lion/

Building from CVS and more:
http://www.gnuplot.info/development/
')
}

}
}

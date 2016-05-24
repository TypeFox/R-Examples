# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"change.subset"  <- function(object, classic = FALSE)
{
  data <- object

  print.flag.help <- function() {
cat("
Subset selection help
---------------------
To specify a subset of the data to process, you use the varable names
and the regular R \"selection\" operators. To combine a subset over two
or more variables, the selection expressions for the two variables are
combined using R\'s unary logical operators.

The variable names are those that are specified in the NONMEM table
files (e.g. PRED, TIME, SEX).

The \"selection\" operators are:
==   (equal)
!=   (not equal)
||   (or)
>    (greater than)
<    (less than)
For example, to specify that TIME less than 24 should be processed, you
type the expression: TIME < 24.

The unary logical operators are:
&    (and)
|    (or)
For example, to specify TIME less than 24 and males (SEX equal to 1), you
type:TIME < 24 & SEX == 1

This subset selection scheme works on all variables, including ID numbers.

The subset selection is not entirely stable. For example, there is no
check that the user enters a valid expression, nor that the user specifies 
existing variable names. An erroneous expression will not become evident 
until a plot is attempted and the expression takes effect. 
")

}
cat("The current subset expression is:\n")
  cat(object@Prefs@Subset,"\n\n")

  cat("Type the expression that will evaluate to the logical vector which\n")
  cat("indicates the subset you want to use. Type h to get help, q to \n")
  cat("leave it as it is or NULL to unspecify:\n")
  

  ans <- readline()
  if(ans == "q")
    return(cat(""))
  if(ans == "n" || ans == "NULL") {
    ans <- NULL
    data@Prefs@Subset <- ans
        if (classic==TRUE) {
          c1<-call("assign",paste("xpdb", object@Runno, sep = ""), data, immediate=T, envir = .GlobalEnv)
          eval(c1)
          c2<-call("assign",pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))
          eval(c2)
          return(cat(""))
          
        } else {
          return(data)
        }
  } else if(ans == "h") {
    print.flag.help()
  } else {
    data@Prefs@Subset <- ans
        if (classic==TRUE) {
          c1<-call("assign",paste("xpdb", object@Runno, sep = ""), data, immediate=T, envir = .GlobalEnv)
          eval(c1)
          c2<-call("assign",pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))
          eval(c2)
          return(cat(""))
        } else {
          return(data)
        }
  }
  
}



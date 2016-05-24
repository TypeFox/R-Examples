### entropy.NSB.R  (2008-08-21)
###
###    R interface to the nsb-entropy estimator
###
### Copyright 2006-2008 Jean Hausser
###
###
### This file is part of the `entropy' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA




##################### public function ##########################################

entropy.NSB = function(y, unit=c("log", "log2", "log10"), CMD="nsb-entropy")
{
  unit = match.arg(unit)

  tmpfile = tempfile()
  tmpfile.txt = paste(tmpfile, "txt", sep=".")
  nsboutfile = nsbsave(filename=tmpfile.txt, y)
  system(paste(CMD, "-dnum -s1 -e1 -cY -v0", tmpfile))
  H = nsbload(filename=nsboutfile)
  unlink(nsboutfile)
  unlink(tmpfile.txt)

  if (unit == "log2")  H = H/log(2)  # change from log to log2 scale
  if (unit == "log10") H = H/log(10) # change from log to log10 scale

  return(H)
}
##################### private function #########################################

# saves a vector of counts (theta) such that it can be read by nsb-entropy
nsbsave = function(filename="samples.txt", theta) 
{
  datastr = character()
  for (i in 1:length(theta))
  {
    if ( theta[i] == 0 ) next
    for (j in 1:theta[i])
    {
      datastr = paste(datastr, i-1)
    }
  }
  #Removing leading space
  datastr = substr(datastr, 2, nchar(datastr))
  
  fileHandle = file(filename, "w")
  cat(
    paste(
      "# name: ASize",
      paste("# type: scalar", length(theta)),
      "# name: data",
      "# type: matrix",
      "# rows: 1",
      paste("# columns:", sum(theta)),
      datastr, 
      sep="\n"),
    file=fileHandle)
  close(fileHandle)

  # return the name nsb-entropy will store the results in
  nsboutfile = character()  
  splitname = strsplit(filename,'.',fixed=TRUE)[[1]]
  return(paste(splitname[1], "_uni_num", length(theta),
                 "_mf1f0_1_entr.", splitname[2], sep=''))
}

# read text file output by nsb-entropy
nsbload = function(filename) 
{
  nsbout = readLines(filename)
  # Seek to the Snsb section
  for (i in 1:length(nsbout)) {
    if ( nsbout[i] == "# name: Snsb" ) break
  }
  #Snsb is at i+4
  return(as.double(nsbout[i+4]))
}

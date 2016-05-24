#Function: ghap.makefile
#License: GPLv3 or later
#Modification date: 2 Feb 2016
#Written by: Marco Milanesi
#Contact: marco.milanesi.mm@gmail.com
#Description: Create a copy of the example file in the working directory

ghap.makefile<-function(){
  
  #Copy files to current directory
  file.copy(system.file("extdata","HapMap3_chr2.tar.bz2", package = "GHap"), "HapMap3_chr2.tar.bz2")
  untar("HapMap3_chr2.tar.bz2",compressed = "bzip2")
  file.remove("HapMap3_chr2.tar.bz2")
  
}
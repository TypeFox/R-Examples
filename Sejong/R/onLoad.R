#Copyright 2012 Heewon Jeon(madjakarta@gmail.com)
#
#This file is part of Sejong.
#
#Sejong is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#Sejong is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.



#.SejongEnv <- new.env()



.onAttach <- function(libname, pkgname){
  packageStartupMessage("Successfully Loaded Sejong Package.")
}



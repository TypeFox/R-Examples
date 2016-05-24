
#Copyright (c) 2009-2014 Sebastien Bihorel
#All rights reserved.
#
#This file is part of scaRabee.
#
#    scaRabee is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    scaRabee is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with scaRabee.  If not, see <http://www.gnu.org/licenses/>.
#

scarabee.clean <- function(files=NULL, analysis=NULL){
  
  options(warn=-1)
  if (file.exists(files$data)) file.remove(files$data)
  if (file.exists(files$param)) file.remove(files$param)
  
  file <- paste(analysis,'.R',sep='')
  if (file.exists(file)) file.remove(file) 
  
  options(warn=0)
  
  setwd('../')
  
}
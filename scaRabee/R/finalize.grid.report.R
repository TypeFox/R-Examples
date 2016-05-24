
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

finalize.grid.report <- function(problem=NULL,
                                 fgrid=NULL,
                                 files=NULL){
  
  # Format fgrid to 6 significant digits and convert to character
  formatgrid <- apply(fgrid,2,function(x) as.character(signif(x,6)))
  
  # Determine the maximum number of character in each column 
  collen <- apply(formatgrid,2,function(x) max(nchar(as.character(x))))
 
  # Determine the number of character in the header of fgrid
  nameslen <- nchar(names(fgrid))
  
  # Determine the width of each column to be printed
  maxlen <- apply(rbind(collen,nameslen),2,max)
  
  # Create format string to print grid
  format <- paste('%',maxlen,'s',sep='')
  
  # Print grid to report
  write('\nGrid summary table\n',file=files$report,append=TRUE,sep='')
  tmp <- paste(sprintf(format,names(maxlen)),collapse=' ')
  write(tmp,file=files$report,append=TRUE,sep='\n')
  
  apply(formatgrid,
        1,
        function(x,...){
          tmp <- paste(sprintf(format,x),collapse=' ')
          write(tmp,file=files$report,append=TRUE,sep='\n')},
        files,format)
  
  write(sprintf('\nGrid search completed at: %s\n', Sys.time()),
        file=files$report,append=TRUE,sep='\n')
      
}

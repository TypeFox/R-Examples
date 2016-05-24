
###################################################################
#
# This function is part of WACSgen V1.0
# Copyright © 2013,2014,2015, D. Allard, BioSP,
# and Ronan Trépos MIA-T, INRA
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. http://www.gnu.org
#
###################################################################




#' Function that reads a file of format Agroclim
#' 
#' @export
#' 
#' @param filename  path to the filemane to read
#' @return A list of 2 objects:
#'         \itemize{
#'           \item{data}{The table data}
#'           \item{mapping}{The names of special var: year, month, day rain,
#'                          tmin and tmax}
#'         }
#'


WACSreadAgroclim = function(filename) 
{
  if (!is.character(filename) || ! file.exists(filename)) {
    stop ("[WACSreadAgroclim] the file does not exist")
  }
  if (! file.exists(filename)) {
    stop (paste("[WACSreadAgroclim] file not found:", filename))
  }
  
  tmplines = readLines(filename, n=100);
  if (length(tmplines) < 1) {
    stop (paste("[WACSreadAgroclim] file bad agroclim format:", filename))
  }
  if (length(grep("Nom de la demande", tmplines[1])) == 0) {
    stop (paste("[WACSreadAgroclim] file bad agroclim format:", filename));
  }
  index_line = grep ("NUM_POSTE", tmplines);
  if ((length(index_line) != 1) || (index_line[1] == 0)) {
    stop (paste("[WACSreadAgroclim] file bad agroclim format:", filename));
  }
  index_line = index_line[1];
  table = read.table(filename, header=TRUE, skip=index_line[1]-1,
                     dec=".", sep=";")
  
  name_year = NULL
  name_month = NULL
  name_day = NULL
  name_tmin = NULL
  name_tmax = NULL
  name_rain = NULL
  for (i in 1:ncol(table)) {
    colName = names(table)[i];
    
    if (colName == "AN") {
      name_year = "AN";
    }
    if (colName == "MOIS") {
      name_month = "MOIS";
    }
    if (colName == "JOUR") {
      name_day = "JOUR";
    }
    if (colName == "TN") {
      name_tmin = "TN";
    }
    if (colName == "TX") {
      name_tmax = "TX";
    }
    if (colName == "RR") {
      name_rain = "RR";
    }
  }
  if (is.null(name_year)){
    stop ("[WACSreadAgroclim] stop year column not found");
  }
  if (is.null(name_month)){
    stop ("[WACSreadAgroclim] stop month column not found");
  }
  if (is.null(name_day)){
    stop ("[WACSreadAgroclim] stop day column not found");
  }
  if (is.null(name_tmin)){
    stop ("[WACSreadAgroclim] stop tmin column not found");
  }
  if (is.null(name_tmax)){
    stop ("[WACSreadAgroclim] stop tmax column not found");
  }
  if (is.null(name_rain)){
    stop ("[WACSreadAgroclim] stop rain column not found");
  }
  res = list();
  res[[name_year]] = "year";
  res[[name_month]] = "month";
  res[[name_day]] = "day";
  res[[name_rain]] = "rain";
  res[[name_tmin]] = "tmin";
  res[[name_tmax]] = "tmax";
  return (list(data=table, mapping=res))
}

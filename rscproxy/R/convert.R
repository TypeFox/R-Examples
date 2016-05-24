# 
#  RSCProxy : Portable, easy-to-use implementation of C-style interface to R
#             (StatConnector)
#  Copyright (C) 2003-2009 Thomas Baier
#
#  convert.R: public conversion functions (exported)
# 
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; version 2 of the License.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# 
as.matrix.rcomdata<-function(x,...)
{
  dims<-dim(x)
  data<-unlist(x)
  dim(data)<-dims
  data
}

as.data.frame.rcomdata<-function(x,row.names = NULL,optional = FALSE,...){
  arglist<-list(...)
  if (!is.null(arglist$col.names)) col.names <- arglist$col.names
  else col.names <- NULL
  if (!is.null(arglist$header)) header <- arglist$header
  else header <- TRUE
  if (header) {
    myColNames <- unlist(x[1,])
    myData <- x[-1,]
  } else {
    myData <- x
  }
  myData <- as.data.frame(lapply(1:dim(myData)[2],function(i)unlist(myData[,i])))
  if (header) names(myData) <- myColNames
  if (!is.null(row.names))  rownames(myData) <- row.names
  if (!is.null(col.names))  colnames(myData) <- col.names
  myData
} 

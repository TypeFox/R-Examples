#
# file: plainNames2data.frameNames.R
#
# This file is part of the R-package decisionSupport
#
# Authors:
#   Lutz Göhring <lutz.goehring@gmx.de>
#   Eike Luedeling (ICRAF) <eike@eikeluedeling.com>
#
# Copyright (C) 2015 World Agroforestry Centre (ICRAF)
#	http://www.worldagroforestry.org
#
# The R-package decisionSupport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The R-package decisionSupport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R-package decisionSupport.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################################
#' Transform model function variable names: plain to data.frame names.
#'
#' The variable names of a function are transformed from plain variable names to data.frame names
#' of the form \code{x$<globalName>}.
#' @param modelFunction a function whose body contains variables with plain names. The 
#'   function must not contain any arguments.
#' @param plainNames a \code{character} vector containing the names of the  variables that
#'   shall be transformed.
#' @details
#'  The input function must be of the form:
#'  \preformatted{
#'    modelFunction<-function(){
#'      ...
#'      <expression with variable1>
#'      ...
#'    }
#'  }
#' @return The transformed function which is of the form:
#'  \preformatted{
#'    function(x){
#'      ...
#'      <expression with x$variable1>
#'      ...
#'    }
#'  }
#' @section Warning:
#'   If there are local functions within the function \code{modelFunction} defined, whose arguments
#'   have identical names to any of the \code{plainNames} the function fails!
#' @examples
#'  profit1<-function(){
#'    list(Profit=revenue-costs)
#'  }
#'  profit2<-plainNames2data.frameNames(modelFunction=profit1, 
#'                                                plainNames=c("revenue", "costs"))
#   #CAVE: here is a problem with the environment: ToDo check!!! (cf. above where eval is used)
#'  print(profit2)
#'  is.function(profit2)
#'  profit2(data.frame("revenue"=10,"costs"=2))
#' @seealso \code{\link{mcSimulation}}, \code{\link{estimate}}
#' @export
plainNames2data.frameNames<-function(modelFunction, plainNames){
	modelFunctionString<-deparse(modelFunction)
	# Replace all occurences of the global variable names:
	for (i in plainNames){
		modelFunctionString<-gsub(pattern=i,
															replacement=paste("x$",i,sep=""),
															x=modelFunctionString,
															fixed=TRUE)
	}
	# Replace only the first occurrence of "()" ; ToDo: regular expression s.t. any "(    )" is replace, i.e. not depending on the
	# count of whitespaces:
	modelFunctionString<-sub(pattern=c("()"),replacement=c("(x)"),x=modelFunctionString, fixed=TRUE)
	#print(modelFunctionString)
	# ToDo: Probably here some problem with the environment emerges:
	modelFunctionData.frameNames<-eval(parse(text=modelFunctionString))
	# Return the transformed function:
	modelFunctionData.frameNames
}
if(0){
	##############################################################################################
	# Test it:
	profit1<-function(){
		list(Profit=revenue-costs)
	}
	profit2<-plainNames2data.frameNames(modelFunction=profit1, plainNames=c("revenue", "costs"))
	# CAVE: here is a problem with the environment: ToDo check!!! (cf. above where eval is used)
	print(profit2)
	is.function(profit2)
	profit2(data.frame("revenue"=10,"costs"=2))
	##############################################################################################
	# Scratch:
	profit1String<-deparse(profit1)
	print(profit1String)
	#profit2String<-cat(paste(sub(pattern=c("revenue"),replacement=c("x$revenue"),x=profit1String),collapse="\n"))
	# Replace all occurences of "revenue":
	profit2String<-gsub(pattern=c("revenue"),replacement=c("x$revenue"),x=profit1String, fixed=TRUE)
	print(profit2String)
	# Replace all occurences of "costs":
	profit2String<-gsub(pattern=c("costs"),replacement=c("x$costs"),x=profit2String, fixed=TRUE)
	print(profit2String)
	# Replace only the first occurrence of "()" ; ToDo: regular expression s.t. any "(    )" is replace, i.e. not depending on the
	# count of whitespaces:
	profit2String<-sub(pattern=c("()"),replacement=c("(x)"),x=profit2String, fixed=TRUE)
	print(profit2String)
	profit2<-eval(parse(text=profit2String))
	print(profit2)
	is.function(profit2)
	profit2(data.frame("revenue"=10,"costs"=2))
}

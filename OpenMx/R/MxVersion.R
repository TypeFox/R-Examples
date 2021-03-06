#
#   Copyright 2007-2015 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

gitVersion <- "v2.5.2"

mxVersion <- function (model=NULL, verbose=TRUE) {
    pvers <- try(packageVersion("OpenMx"))
    if ("try-error" %in% class(pvers)) {
        pvers = gitVersion
    }
	if(verbose){
		msg = paste("OpenMx version: ", pvers, " [GIT ", gitVersion, "]", sep="")
		msg = paste(msg, "\nR version: ", version$version.string, sep="")
		msg = paste(msg, "\nPlatform: ", version$platform, sep="")
		if ("Darwin" ==Sys.info()["sysname"]){
			msg = paste(msg, "\nMacOS:", system("sw_vers -productVersion", intern=TRUE))
		}
		msg = paste(msg, "\nDefault optimiser: ", mxOption(NULL, "Default optimizer"), sep="")
		if (!is.null(model)) {
			thisModelsOptimiser = mxOption(model, "Default optimizer")
		    if(!is.null(thisModelsOptimiser)){
				msg = paste(msg, "(optimizer for this model is ", thisModelsOptimiser, ")", sep="")
		    }
		}
	    message(msg)
	}
	invisible(pvers)
}

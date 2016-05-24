##
#  Copyright (c) 2008-2010 Fred Hutchinson Cancer Research Center
# 
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
##

parseHeader <- function (lines)
{
    if (length(lines) < 1)
        return(NULL)
    if (length(lines) == 1)
        lines = strsplit(lines, "\r\n")[[1]]
    status = lines[1]
    lines = lines[-c(1, length(lines))]
    lines = gsub("\r\n", "", lines)
    if (FALSE) {
        header = lines[-1]
        header <- read.dcf(textConnection(header))
    }
    else {
        els <- sapply(lines, function(x) strsplit(x, ":[ ]*"))
        header <- lapply(els, function(x) x[2])
        names(header) <- sapply(els, function(x) x[1])
    }
    els <- strsplit(status, " ")[[1]]
    header[["status"]] <- as.integer(els[2])
	hstring <- NULL
	for(i in 3:length(els)) hstring <- paste(hstring," ",els[i],sep="")
	hstring <- substr(hstring,2,nchar(hstring))
    header[["statusMessage"]] <- hstring
    header
}


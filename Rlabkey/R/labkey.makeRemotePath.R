##
# Copyright (c) 2012 LabKey Corporation
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##

labkey.makeRemotePath <- function(localRoot, remoteRoot, fullPath)
{

## we need a full path to get started
if (nzchar(fullPath) == FALSE)
{
    stop(paste("A non-empty value must be specified for fullPath"));
}

## normalize the fullPath right now to ensure consistent output
fullPath <- gsub("[\\]", "/", fullPath);

## if we have an empty remote path then just return what the user passed in
if (nzchar(remoteRoot) == FALSE)
{
    return(fullPath);
}


## normalize the roots
localRoot <- gsub("[\\]", "/", localRoot);
remoteRoot <- gsub("[\\]", "/", remoteRoot);

## ensure roots have trailing "/"
if(substr(localRoot, nchar(localRoot), nchar(localRoot))!="/")
{
    localRoot <- paste(localRoot,"/",sep="")
}
if(substr(remoteRoot, nchar(remoteRoot), nchar(remoteRoot))!="/")
{
    remoteRoot <- paste(remoteRoot,"/",sep="")
}

## if the local path and remote path are the same then just return the fullPath
if ( localRoot %in% remoteRoot )
{
    return(fullPath);
}


## do the replacement only if the localRoot is part of the filePath
if (localRoot %in% substr(fullPath, 1, nchar(localRoot)))
{
    return(paste(remoteRoot, substr(fullPath, nchar(localRoot) + 1, nchar(fullPath)), sep=""));
}

## otherwise just return what was passed in
return(fullPath);
}
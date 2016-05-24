##
# Copyright (c) 2010 LabKey Corporation
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

fromJSON2 <-
    function( json_str )
{
    .Call(".JSON_to_R", json_str, PACKAGE="Rlabkey")
}

newJSONParser2 <-
    function()
{
    p <- .Call(".parser_new", PACKAGE="Rlabkey")
    class(p) <- "jsonParser2"
    p
}   

`$.jsonParser2` <-
    function( x, name )
{
    switch(name,
           addObject=function( json_str ) {
               invisible(.Call("parser_add", x, json_str, TRUE,
                               PACKAGE="Rlabkey"))
           },
           getObject=function() {
               .Call("parser_finalize", x, PACKAGE="Rlabkey")
           },
           reset=function() {
               invisible(.Call("parser_delete", x, PACKAGE="Rlabkey"))
           },
           stop("unknown parser action '", name, "'"))
}

print.jsonParser2 <-
    function( x, ... )
{
    cat("json parser2\n")
}

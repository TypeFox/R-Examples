## Copyright 2013-2014 Stefan Widgren and Maria Noremark,
## National Veterinary Institute, Sweden
##
## Licensed under the EUPL, Version 1.1 or - as soon they
## will be approved by the European Commission - subsequent
## versions of the EUPL (the "Licence");
## You may not use this work except in compliance with the
## Licence.
## You may obtain a copy of the Licence at:
##
## http://ec.europa.eu/idabc/eupl
##
## Unless required by applicable law or agreed to in
## writing, software distributed under the Licence is
## distributed on an "AS IS" basis,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
## express or implied.
## See the Licence for the specific language governing
## permissions and limitations under the Licence.

library(EpiContactTrace)

##
## Node position checking
##

##
## Case 1
##
tree <- structure(list(
    node = c("O", "E", "F", "N", "A", "D", "G", "M", "B", "C", "H", "I", "J", "K", "L"),
    parent = c(NA, "O", "O", "O", "E", "E", "N", "N", "D", "D", "M", "M", "M", "M", "M"),
    level = c(0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3)),
                  .Names = c("node", "parent", "level"),
                  row.names = c(NA, -15L), class = "data.frame")
tree.exp <- structure(list(
    node = c("O", "E",  "F", "N", "A", "D", "G", "M", "B", "C", "H", "I", "J", "K", "L"),
    parent = c(NA, "O", "O", "O", "E", "E", "N", "N", "D", "D", "M", "M", "M", "M", "M"),
    level = c(0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3),
    x =  c(0, -10.5, 0, 10.5, -13.5, -7.5, 7.5, 13.5, -10.5, -4.5, 1.5, 7.5, 13.5, 19.5, 25.5),
    y = c(0, -1, -1, -1, -2, -2, -2, -2, -3, -3, -3, -3, -3, -3, -3)),
    .Names = c("node", "parent", "level", "x", "y"), row.names = c(NA, -15L), class = "data.frame")
tree.obs <- EpiContactTrace:::position_tree(tree)
tree.obs$level <- as.numeric(tree.obs$level)
str(tree.obs)
stopifnot(identical(tree.obs, tree.exp))

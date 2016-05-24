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
data(transfers)

##
## Check NetworkSummary
##

##
## Case 1
##
load(file=system.file("extdata", "ns.rda", package="EpiContactTrace"))
root <- sort(unique(c(transfers$source, transfers$destination)))
result <- NetworkSummary(transfers, root=root, tEnd='2005-10-31', days=90)
stopifnot(identical(result, ns))

##
## Case 2
##
ns <- NetworkSummary(transfers, root=584, tEnd='2005-10-31', days=91)
ns.trace <- NetworkSummary(Trace(transfers,
                                 root=584,
                                 tEnd='2005-10-31',
                                 days=91))
stopifnot(identical(ns, ns.trace))

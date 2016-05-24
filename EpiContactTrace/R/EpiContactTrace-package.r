## Copyright 2013 Stefan Widgren and Maria Noremark,
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

##' Epidemiological tool for contact tracing.
##'
##' Routines for epidemiological contact tracing and visualisation of network
##' of contacts.
##'
##'
##' In many countries, livestock movement data are
##' collected with the major objective to enable contact tracing during disease
##' outbreaks. Livestock movement data can also be of relevance for risk based
##' surveillance - both during outbreak or when investigating if a disease is
##' present in the population. However, the livestock movement databases are
##' not always structured in such a way that relevant information for contact
##' tracing or surveillance design is easily retrieved. EpiContactTrace uses
##' the network parameters in-degree, out-degree, ingoing contact-chain and
##' outgoing contact-chain, which are relevant for forward- and backward
##' contact-tracing respectively. The measures can also be used for identifying
##' herds with many contacts, which can be used in risk based disease
##' surveillance. Different time periods for ingoing and outgoing contacts can
##' be of interest in the contact tracing, based on possible window of
##' introduction, and this can be adjusted in the tool. The output from the
##' analysis is available as a dataset, but moreover, the tool automatically
##' generates a report on farm level. The report both contains an overview of
##' the situation on the farm, including a graph, as well as detailed
##' information including dates of movements on group or individual level on
##' all contacts.
##'
##' @name EpiContactTrace-package
##' @aliases EpiContactTrace-package EpiContactTrace
##' @docType package
##' @author Stefan Widgren Maria Noremark
##' @section Maintainer:
##' Stefan Widgren <stefan.widgren@@sva.se>
##' @references \itemize{
##'   \item Dube, C., et al., A review of network analysis terminology
##'     and its application to foot-and-mouth disease modelling and policy
##'     development. Transbound Emerg Dis 56 (2009) 73-85, doi:
##'     10.1111/j.1865-1682.2008.01064.x
##'
##'   \item Noremark, M., et al., Network analysis
##'     of cattle and pig movements in Sweden: Measures relevant for
##'     disease control and riskbased surveillance.  Preventive Veterinary
##'     Medicine 99 (2011) 78-90, doi: 10.1016/j.prevetmed.2010.12.009
##' }
##' @note Plots are not supported in version 0.8.6 since igraph0 has
##' been archived. We intend to resolve the issue in a future
##' version. Install version 0.8.5 and igraph0 manually from the
##' archive if plots are required. See section 6.3 in 'R Installation
##' and Administration' on how to install packages from source.
##' @keywords package
##' @import methods
##' @useDynLib EpiContactTrace
##' @examples
##'
##' ## Load data
##' data(transfers)
##'
##' ## Perform contact tracing
##' contactTrace <- Trace(movements=transfers,
##'                       root=2645,
##'                       tEnd='2005-10-31',
##'                       days=90)
##'
##' show(contactTrace)
##'
##' \dontrun{
##' ## Generate an html report showing details of the contact tracing for
##' ## root 2645.
##' ## Note: Creates the files 2645.html and 2645.png in the working
##' ## directory.
##' Report(contactTrace)
##' }
##'
NULL

##' Movement Example Data
##'
##' Movement data included in the package. The data contains fictitious
##' example data of cattle movements during the period 2005-08-01 -- 2005-10-31.
##'
##'
##' @name transfers
##' @docType data
##' @usage data(transfers)
##' @format A data frame with 70190 observations on the following 6 variables.
##' \describe{
##'   \item{source}{
##'     a numeric vector with the holding identifier of the source
##'   }
##'
##'   \item{destination}{
##'     a numeric vector with holding identifier of the destination.
##'   }
##'
##'   \item{id}{
##'     a character vector with the identity of the animal.
##'     In this dataset an 5 character hexadecimal vector.
##'   }
##'
##'   \item{t}{
##'     a Date of the transfers
##'   }
##'
##'   \item{n}{
##'     a numeric vector with the number of animals moved.
##'     Always 1 in this dataset.
##'   }
##'
##'   \item{category}{
##'     a factor describing the category of the animal.
##'     Always \code{Cattle} in this dataset.
##'   }
##' }
##' @keywords datasets
##' @examples
##'
##' \dontrun{
##' data(transfers)
##'
##' contactTrace <- Trace(movements=transfers,
##'                       root=2645,
##'                       tEnd='2005-10-31',
##'                       days=90)
##'
##' show(contactTrace)
##' }
##'
NULL

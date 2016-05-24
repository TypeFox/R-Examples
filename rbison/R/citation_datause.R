#' Get BISON data use agreement details and examples for how to cite data.
#' 
#' @export
#' @references \url{http://bison.usgs.ornl.gov/doc/api.jsp#data}
#' @rdname bison_datause
bison_datause <- function(){
  cat('
  Any use of BISON web services automatically commits the user to acceptance of the Data 
    Use Agreement.
  
  Adapted from the Data Use Agreement of the Global Biodiversity Information Facility (GBIF).
  
  Background
  
  Biodiversity Information Serving Our Nation (BISON) is committed to providing free and open 
  access to species occurrence data. Data currently available through BISON are provided through 
  the Global Biodiversity Information Facility (GBIF). GBIF Participants who have signed the GBIF 
  Memorandum of Understanding have expressed their willingness to make biodiversity data available 
  through their nodes to foster scientific research development internationally and to support the 
  public use of these data (see http://data.gbif.org/tutorial/datasharingagreement). GBIF data 
  sharing should take place within a framework of due attribution. Therefore, using data available 
  through BISON requires agreement with the following:
    
    1. Data Use Agreements
  
  The quality and completeness of data cannot be guaranteed. Users employ these data at their 
  own risk. Users shall respect restrictions of access to sensitive data. In order to make 
  attribution of use for owners of the data possible, the identifier of ownership of data must be 
  retained with every data record. Users must publicly acknowledge, in conjunction with the 
  use of the data, the data publishers whose biodiversity data they have used. Data publishers 
  may require additional attribution of specific collections within their institution. Users must 
  comply with additional terms and conditions of use set by the data publisher. Where these exist 
  they will be available through the metadata associated with the data 
  (see http://data.gbif.org/datasets/). 
  ')
}

#' @export
#' @rdname bison_datause
bison_citation <- function(){
  cat('
  Use the following format to cite data retrieved from BISON:
    [Data Provider or owner name]. [Resource or Dataset Name] published by [Data Provider 
    name, address or affiliation(s)] (Accessed through Biodiversity Information Serving Our 
    Nation (BISON), http://bison.usgs.ornl.gov, YYYY-MM-DD)
  
  For example: 
    Field Museum of Natural History. U.S. Bird Occurrences published by Field Museum of 
    Natural History, Museum of Vertebrate Zoology, University of Washington Burke Museum, 
    and University of Turku (Accessed through Biodiversity Information Serving Our Nation 
    (BISON), http://bison.usgs.ornl.gov, 2007-02-22)
  
  Or
  
    Gordon, J. U.S. Bird Occurrences published by Field Museum of Natural History, Museum 
    of Vertebrate Zoology, University of Washington Burke Museum, and University of Turku 
    (Accessed through Biodiversity Information Serving Our Nation (BISON), 
    http://bison.usgs.ornl.gov, 2007-02-22)'
  )
}

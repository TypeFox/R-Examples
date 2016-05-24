#############################################################################
#
# Copyright Institut Télécom-Télécom Bretagne, 2012
#
# Contributors:
#   Patrick Meyer <patrick.meyer@telecom-bretagne.eu>
#   Sebastien Bigaret <sebastien.bigaret@telecom-bretagne.eu>
#		
# This software, RXMCDA, is a library to for the R statistical software which 
# allows you to handle XMCDA tags and transform them into R variables. 
# 
# This software is governed by the CeCILL license (v2) under French law
# and abiding by the rules of distribution of free software. You can
# use, modify and/ or redistribute the software under the terms of the
# CeCILL license as circulated by CEA, CNRS and INRIA at the following
# URL "http://www.cecill.info".
# 
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#		
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#		
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
##############################################################################


checkXSD <- function(tree){
  
  # get the namespaces defined in tree
  
  namespaces<-xmlNamespaces(tree, simplify=TRUE)
  
  # search for namespaces containing XMCDA-2.*
  
  i <- grep("XMCDA-2.",namespaces)	
  
  # xsdLocations <- c("http://www.decision-deck.org/2009/XMCDA-2.0.0" = "http://www.decision-deck.org/xmcda/_downloads/XMCDA-2.0.0.xsd", "http://www.decision-deck.org/2009/XMCDA-2.1.0" = "http://www.decision-deck.org/xmcda/_downloads/XMCDA-2.1.0.xsd")
  
  # the schema are locally stored
  
  xsdLocations <- c("http://www.decision-deck.org/2009/XMCDA-2.0.0" = "XMCDA-2.0.0.xsd", 
                    "http://www.decision-deck.org/2009/XMCDA-2.1.0" = "XMCDA-2.1.0.xsd", 
                    "http://www.decision-deck.org/2012/XMCDA-2.2.0" = "XMCDA-2.2.0.xsd", 
                    "http://www.decision-deck.org/2012/XMCDA-2.2.1" = "XMCDA-2.2.1.xsd")
  
  if (!is.na(xsdLocations[namespaces[i[1]]])){
    # xsd <- xmlTreeParse(xsdLocations[namespaces[i[1]]], isSchema =TRUE, useInternalNodes = TRUE)
    xsd <- xmlTreeParse(system.file("extdata",xsdLocations[namespaces[i[1]]],package="RXMCDA"), isSchema =TRUE, useInternalNodes = TRUE)
  }else{
    # xsd <- xmlTreeParse("http://www.decision-deck.org/xmcda/_downloads/XMCDA-2.1.0.xsd", isSchema =TRUE, useInternalNodes = TRUE)
    xsd <- xmlTreeParse(system.file("extdata","XMCDA-2.2.1.xsd",package="RXMCDA"), isSchema =TRUE, useInternalNodes = TRUE)
  }
  
  if (xmlSchemaValidate(xsd,tree)$status != 0)
    return(0)
  else
    return(1)
}

# getNumericValue returns the first numeric value read in an XML tree under the value tag (real, integer, rational).
# If none is found, NA is returned.

getNumericValue <- function(tree){
  
  out<-NA
  
  if (names(xmlChildren(tree[[1]]))[1] == "real"){
    out<-as.double(xmlValue(getNodeSet(tree[[1]], "real")[[1]]))
  }
  else if (names(xmlChildren(tree[[1]]))[1] == "integer"){
    out<-as.integer(xmlValue(getNodeSet(tree[[1]], "integer")[[1]]))
  }
  else if (names(xmlChildren(tree[[1]]))[1] == "rational"){
    num<-as.integer(xmlValue(getNodeSet(tree[[1]], "rational/numerator")[[1]]))
    den<-as.integer(xmlValue(getNodeSet(tree[[1]], "rational/denominator")[[1]]))
    out<-num/den
  }
  
  return(out)
}

# getNumberOfCriteria returns a list containing the number of criteria defined in each <criteria> tag. 
# Possibility to specify which mcdaConcept should be searched. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getNumberOfCriteria <- function(tree, mcdaConcept = NULL){
  
  err<-NULL
  
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the criteria tag(s) from the tree (according to the mcdaConcept if necessary)
  criteria <- getNodeSet(tree, paste("//criteria",specification,sep=""))
  
  # create the empty output list
  out<-list()
  
  # count the number of <criterion> in each <criteria> and add it to the out list
  if (length(criteria)>0){
    for (i in 1:length(criteria)){
      elements <- getNodeSet(criteria[[i]], "criterion")
      inactive <- getNodeSet(criteria[[i]], "criterion[active='false']")
      out <- c(out,list(length(elements)-length(inactive)))
      names(out)[length(out)]<-toString(xmlGetAttr(criteria[[i]],"mcdaConcept"))
    }
  }
  else { #if (length(criteria)>0){
    err<-"No <criteria> found."
  }
  # if there is no error, print status = OK, else print status = description of the error
  if (!is.null(err)){
    out<-c(out,list(status=err))
  }
  else{ # if (!is.null(err)){ 
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# getCriteriaIDs returns a list containing the ids of the criteria in each <criteria> tag. 
# Possibility to specify which mcdaConcept should be searched.
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getCriteriaIDs <- function(tree, mcdaConcept = NULL){
  
  err<-NULL
  
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the criteria tag(s) from the tree (according to the mcdaConcept if necessary)
  criteria <- getNodeSet(tree, paste("//criteria",specification,sep=""))
  
  # create the empty output list
  out<-list()
  
  # extract the @id of each criterion
  if (length(criteria)>0){
    for (i in 1:length(criteria)){
      elements <- getNodeSet(criteria[[i]], "criterion")
      criteriaIDs <- c()
      if (length(elements)>0){
        for (j in 1:length(elements)){
          # filter for inactive criteria
          act<-getNodeSet(elements[[j]], "active")
          if (length(act)==0){
            # no actvie tag found, therefore supposed to be active
            criteriaIDs<-c(criteriaIDs,xmlGetAttr(elements[[j]], "id"))
          } else {
            if (xmlValue(act[[1]])=="true"){
              # an active tag found which is true
              criteriaIDs<-c(criteriaIDs,xmlGetAttr(elements[[j]], "id"))
            }
          }
        }
      }
      out <- c(out,list(criteriaIDs))
      names(out)[length(out)]<-toString(xmlGetAttr(criteria[[i]],"mcdaConcept"))
    }
  }
  else { #if (length(criteria)>0){
    err<-"No <criteria> found."
  }
  # if there is no error, print status = OK, else print status = description of the error
  if (!is.null(err)){
    out<-c(out,list(status=err))
  }
  else{ #if (!is.null(err)){
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# getNumberOfAlternatives returns a list containing the number of alternatives defined in each <alternatives> tag.
# Possibility to specify which mcdaConcept should be searched.
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getNumberOfAlternatives <- function(tree, mcdaConcept = NULL){
  
  err<-NULL
  
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <alternatives> from the tree (according to the mcdaConcept if necessary)
  alternatives <- getNodeSet(tree, paste("//alternatives",specification,sep=""))
  
  # create the empty output list
  out<-list()
  
  # count the number of <alternative> in each <alternatives> and add it to the out list
  if (length(alternatives)>0){
    for (i in 1:length(alternatives)){
      elements <- getNodeSet(alternatives[[i]], "alternative")
      inactive <- getNodeSet(alternatives[[i]], "alternative[active='false']")
      out <- c(out,list(length(elements)-length(inactive)))
      names(out)[length(out)]<-toString(xmlGetAttr(alternatives[[i]],"mcdaConcept"))
    }
  }
  else { #if (length(alternatives)>0){
    err<-"No <alternatives> found."
  }
  # if there is no error, print status = OK, else print status = description of the error
  if (!is.null(err)){
    out<-c(out,list(status=err))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  
  return(out)
}

# getAlternativesIDs returns a list containing the ids of the alternatives in each <alternatives> tag.
# Possibility to specify which mcdaConcept should be searched. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getAlternativesIDs <- function(tree, mcdaConcept = NULL){
  
  err <-NULL
  
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <alternatives> from the tree (according to the mcdaConcept if necessary)
  alternatives <- getNodeSet(tree, paste("//alternatives",specification,sep=""))
  
  # create the empty output list
  out<-list()
  
  # extract the @id of each alternative
  if (length(alternatives)>0){
    for (i in 1:length(alternatives)){
      elements <- getNodeSet(alternatives[[i]], "alternative")
      alternativesIDs <- c()
      if (length(elements)>0){
        for (j in 1:length(elements)){
          # don't consider inactive alternatives
          act<-getNodeSet(elements[[j]], "active")
          if (length(act)==0){
            # no actvie tag found, therefore supposed to be active
            alternativesIDs<-c(alternativesIDs,xmlGetAttr(elements[[j]], "id"))
          } else {
            if (xmlValue(act[[1]])=="true"){
              # an active tag found which is true
              alternativesIDs<-c(alternativesIDs,xmlGetAttr(elements[[j]], "id"))
            }
          }
        }
      }
      out <- c(out,list(alternativesIDs))
      names(out)[length(out)]<-toString(xmlGetAttr(alternatives[[i]],"mcdaConcept"))
    }
  }
  else { #if (length(alternatives)>0){
    err<-"No <alternatives> found."
  }
  # if there is no error, print status = OK, else print status = description of the error
  if (!is.null(err)){
    out<-c(out,list(status=err))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# getNumberOfCategories returns a list containing the number of categories defined in each <categories> tag.
# Possibility to specify which mcdaConcept should be searched.
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getNumberOfCategories <- function (tree, mcdaConcept = NULL) 
{
  err <- NULL
  specification = ""
  if (!is.null(mcdaConcept)) 
    specification <- paste("[@mcdaConcept='", mcdaConcept, 
                           "']", sep = "")
  categories <- getNodeSet(tree, paste("//categories", specification, 
                                       sep = ""))
  out <- list()
  if (length(categories) > 0) {
    for (i in 1:length(categories)) {
      elements <- getNodeSet(categories[[i]], "category")
      inactive <- getNodeSet(categories[[i]], "category[active='false']")
      out <- c(out, list(length(elements) - length(inactive)))
      names(out)[length(out)] <- toString(xmlGetAttr(categories[[i]], "mcdaConcept"))
    }
  }
  else {
    err <- "No <categories> found."
  }
  if (!is.null(err)) {
    out <- c(out, list(status = err))
  }
  else {
    out <- c(out, list(status = "OK"))
  }
  return(out)
}

# getCategoriesIDs returns a list containing the ids of the categories in each <categories> tag.
# Possibility to specify which mcdaConcept should be searched. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getCategoriesIDs <- function (tree, mcdaConcept = NULL) 
{
  err <- NULL
  specification = ""
  if (!is.null(mcdaConcept)) 
    specification <- paste("[@mcdaConcept='", mcdaConcept, 
                           "']", sep = "")
  categories <- getNodeSet(tree, paste("//categories", specification, 
                                       sep = ""))
  out <- list()
  if (length(categories) > 0) {
    for (i in 1:length(categories)) {
      elements <- getNodeSet(categories[[i]], "category")
      categoriesIDs <- c()
      if (length(elements) > 0) {
        for (j in 1:length(elements)) {
          act <- getNodeSet(elements[[j]], "active")
          if (length(act) == 0) {
            categoriesIDs <- c(categoriesIDs, xmlGetAttr(elements[[j]], "id"))
          }
          else {
            if (xmlValue(act[[1]]) == "true") {
              categoriesIDs <- c(categoriesIDs, xmlGetAttr(elements[[j]], "id"))
            }
          }
        }
      }
      out <- c(out, list(categoriesIDs))
      names(out)[length(out)] <- toString(xmlGetAttr(categories[[i]], "mcdaConcept"))
    }
  }
  else {
    err <- "No <criteria> found."
  }
  if (!is.null(err)) {
    out <- c(out, list(status = err))
  }
  else {
    out <- c(out, list(status = "OK"))
  }
  return(out)
}



# getParameters returns a list containing the values of the parameters under <methodParameters>. 
# Possibility to specify which parameter name should be searched. 
# The elements of the list are named according to the name attribute, if it has been defined. 
# Currently, an parameter without a value below is not read.

getParameters <- function(tree, name = NULL){
  
  err<-NULL
  err1<-NULL
  # if an name has been specified, search according to this attribute
  specification = ""
  if (!is.null(name)) specification <- paste("[@name='",name,"']",sep="")	
  
  # extract the <parameter> from the tree (according to the mcdaConcept if necessary)
  options <- getNodeSet(tree, paste("//methodParameters/parameter",specification,sep=""))
  
  # create the empty output list
  out<-list()
  
  # extract the "value" of each parameter
  if (length(options)>0){
    for (i in 1:length(options)){
      if (length(xmlChildren(options[[i]]))>0){
        
        value <- getNodeSet(options[[i]], "value")
        if (names(xmlChildren(value[[1]]))[1] == "label"){
          out<-c(out,list(xmlValue(getNodeSet(value[[1]], "label")[[1]])))
          names(out)[length(out)]<-toString(xmlGetAttr(options[[i]],"name"))
        }
        else if (names(xmlChildren(value[[1]]))[1] == "real"){
          out<-c(out,list(as.double(xmlValue(getNodeSet(value[[1]], "real")[[1]]))))
          names(out)[length(out)]<-toString(xmlGetAttr(options[[i]],"name"))
        }
        else if (names(xmlChildren(value[[1]]))[1] == "integer"){
          out<-c(out,list(as.integer(xmlValue(getNodeSet(value[[1]], "integer")[[1]]))))
          names(out)[length(out)]<-toString(xmlGetAttr(options[[i]],"name"))
        }
        else if (names(xmlChildren(value[[1]]))[1] == "boolean") {
          val <- xmlValue(getNodeSet(value[[1]], "boolean")[[1]])
          if ((val=="true")||(val=="1")) {
            out <- c(out, TRUE)
          }
          else if  ((val=="false")||(val=="0")) {
            out <- c(out, FALSE)
          } else {
            out <- c(out, list(NA))
          }
          names(out)[length(out)] <- toString(xmlGetAttr(options[[i]], "name"))
        }
        # If the value is neither of types label, 
        # real, integer or boolean, we do not read it. 
        # TODO: add the other types for more universality!
        
      }
    }
  }
  else { #if (length(options)>0){
    err<-"No <methodParameters> found."
  }
  # if there is no error, print status = OK, else print status = description of the error
  if (!is.null(err)){
    out<-c(out,list(status=err))
  }
  else{
    
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# getPerformanceTables returns a list containing the performance tables. 
# Possibility to specify which mcdaConcept should be searched. 
# If altIDs or critIDs are specified, only those rows or columns are extracted.
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getPerformanceTables <- function(tree, altIDs = NULL, critIDs = NULL, mcdaConcept = NULL){
  
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <performanceTable> from the tree (according to the mcdaConcept if necessary)
  performanceTables <- getNodeSet(tree, paste("//performanceTable",specification,sep=""))
  
  # create the empty output list and the errors
  out<-list()
  err1<-NULL
  err2<-NULL
  err3<-NULL
  err4<-NULL
  err5<-NULL
  err6<-NULL
  performance.table<-NULL
  
  if (length(performanceTables)>0){
    for (i in 1:length(performanceTables)){
      alternativePerformances <- getNodeSet(performanceTables[[i]], "alternativePerformances")
      if (length(alternativePerformances)>0){
        
        # We construct the empty performance.table on basis of the number
        # of number of alternatives and the number of criteria declared for the 
        # first alternative. This might be questionnable, but we had to make a choice.
        
        performance <- getNodeSet(alternativePerformances[[1]], "performance")
        
        if (length(performance)>0){
          
          # numAlt and numCrit contain the supposed number of alternatives and criteria
          
          alternativesIDs<-c()
          for (j in 1:length(alternativePerformances)){
            tmp<-getNodeSet(alternativePerformances[[j]], "alternativeID")
            tmpErr<-try(
              alternativesIDs<-c(alternativesIDs,xmlValue(tmp[[1]]))
            )
            if (inherits(tmpErr, 'try-error')){
              err1<-"At least one <alternativePerformances> contains no <alternativeID>."
            }
          }
          criteriaIDs<-c()
          for (j in 1:length(performance)){
            tmp<-getNodeSet(performance[[j]], "criterionID")
            tmpErr<-try(
              criteriaIDs<-c(criteriaIDs,xmlValue(tmp[[1]]))
            )
            if (inherits(tmpErr, 'try-error')){
              err1<-"The first <alternativePerformances> of a <performanceTable> contains no <criterionID>."
            }
          }
          
          tmpErr<-try(
            performance.table<-matrix(nrow=length(alternativesIDs),ncol=length(criteriaIDs),dimnames = list(alternativesIDs,criteriaIDs))
          )
          if (inherits(tmpErr, 'try-error')){
            err2<-"Impossible to create a performance table."
          }
          
          for (j in 1:length(alternativePerformances)){
            alt<-getNodeSet(alternativePerformances[[j]], "alternativeID")
            perf <- getNodeSet(alternativePerformances[[j]], "performance")
            for (j in 1:length(perf)){
              tmpErr<-try(
{
  crit <- getNodeSet(perf[[j]], "criterionID")
  #										value <- getNodeSet(perf[[j]], "value/real")
  #										performance.table[xmlValue(alt[[1]]),xmlValue(crit[[1]])] <- as.numeric(xmlValue(value[[1]]))
  value <- getNodeSet(perf[[j]], "value")
  performance.table[xmlValue(alt[[1]]),xmlValue(crit[[1]])] <- getNumericValue(value)[[1]]
}
              )
              if (inherits(tmpErr, 'try-error')){
                err2<-"Impossible to read (a) value(s) in a <performanceTable>."
              }
            }	
          }
        }
        else{ #if (length(performance)>0){
          err3<-"The first <alternativePerformances> of a <performanceTable> contains no <performance>."
        }
      }
      else{ #if (length(alternativePerformances)>0){
        err4<-"A <performanceTable> contains no <alternativePerformances>."
      }
      # if criteriaIDs and alternativesIDs have been specified, we want to filter the performance
      # table and retain only those lines and columns which have been specified in criteriaIDs and
      # alternativesIDs.
      
      if ((!is.null(critIDs))|(!is.null(altIDs))){
        # print(altIDs)
        # print(rownames(performance.table))
        # print(rownames(performance.table)%in%altIDs)
        # if ((TRUE%in%(rownames(performance.table)%in%altIDs))&(TRUE%in%(colnames(performance.table)%in%critIDs)))
        performance.table<-performance.table[rownames(performance.table)%in%altIDs,colnames(performance.table)%in%critIDs]
      }		
      
      if (length(performance.table)>0){
        out<-c(out,list(performance.table))
        names(out)[length(out)]<-toString(xmlGetAttr(performanceTables[[i]],"mcdaConcept"))	
      }
      else {
        err6<-"No performance table could be read with constraints from criteria IDs and alternatives IDs."
      }
      
    }
    
  }
  else {#if (length(performanceTables)>0){
    err5<-"No <performanceTable> found."
  }
  if (!is.null(err1)|(!is.null(err2))|(!is.null(err3))|(!is.null(err4))|(!is.null(err5))|(!is.null(err6))){
    out<-c(out,list(status=c(err1,err2,err3,err4,err5,err6)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# Returns a list containing the comparisons of single criteria. 
# Possibility to specify which mcdaConcept should be searched. 
# criteriaIDs is used to filter on which criteria ids the extraction is to be done. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getCriteriaComparisons <- function(tree, criteriaIDs, mcdaConcept = NULL){
  
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <performanceTable> from the tree (according to the mcdaConcept if necessary)
  criteriaComparisons <- getNodeSet(tree, paste("//criteriaComparisons",specification,sep=""))
  
  # create the empty output list and the errors
  out<-list()
  err1<-NULL
  err2<-NULL
  
  
  if (length(criteriaComparisons)>0){
    for (i in 1:length(criteriaComparisons)){
      
      # check whether we only have <criteriaID> under <initial> and <terminal>
      test1<-getNodeSet(criteriaComparisons[[i]], "pairs/pair/initial")
      test2<-getNodeSet(criteriaComparisons[[i]], "pairs/pair/terminal")
      test1.names<-NULL
      test2.names<-NULL
      tmpErr<-try(
{
  for (k in 1:length(test1))
    test1.names<-c(test1.names,names(xmlChildren(test1[[k]])))
  for (k in 1:length(test2))
    test2.names<-c(test2.names,names(xmlChildren(test2[[k]])))
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to read (a) value(s) in a <criteriaComparisons>."
      }
      if (!("criteriaSet" %in% test1.names) & !("criteriaSet" %in% test2.names)){
        
        # if there are no <criteriaSet> under <initial> and <terminal>, 
        # we can suppose that this <criteriaComparison> only contains
        # comparisons of single alternatives
        
        criteriaComp <- matrix(nrow=0,ncol=3)
        pairs <- getNodeSet(criteriaComparisons[[i]], "pairs/pair")
        
        if (length(pairs)>0){
          for (j in 1:length(pairs)){
            
            head<-NULL
            tail<-NULL
            val<-NULL
            
            initial <- getNodeSet(pairs[[j]],"initial")
            terminal <- getNodeSet(pairs[[j]],"terminal")
            
            if (names(xmlChildren(initial[[1]]))[1] == "criterionID"){
              # the comparisons are on single criteria
              tmpErr<-try(
{
  head <- getNodeSet(pairs[[j]], "initial/criterionID")
  tail <- getNodeSet(pairs[[j]], "terminal/criterionID")
}
              )
              if (inherits(tmpErr, 'try-error')){
                err2<-"Impossible to read (a) value(s) in a <criteriaComparisons>."
              }
              tmpErr<-try(
{
  #										val <- getNodeSet(pairs[[j]], "value/real")
  #										if ((length(which(criteriaIDs==xmlValue(head[[1]])))>0) & length(which(criteriaIDs==xmlValue(tail[[1]]))>0)) 
  #											criteriaComp <-  rbind(criteriaComp,c(which(criteriaIDs==xmlValue(head[[1]])),which(criteriaIDs==xmlValue(tail[[1]])),as.numeric(xmlValue(val[[1]]))))
  val <- getNodeSet(pairs[[j]], "value")
  if ((length(which(criteriaIDs==xmlValue(head[[1]])))>0) & length(which(criteriaIDs==xmlValue(tail[[1]]))>0)) 
    criteriaComp <-  rbind(criteriaComp,c(which(criteriaIDs==xmlValue(head[[1]])),which(criteriaIDs==xmlValue(tail[[1]])),getNumericValue(val)))
}
              )
              if (inherits(tmpErr, 'try-error')){
                err2<-"Impossible to read (a) value(s) in a <criteriaComparisons>."
              }
            }		
          }
        } #if (length(pairs)>0){
        if (dim(criteriaComp)[1] == 0)
          criteriaComp <- NULL
        out<-c(out,list(criteriaComp))
        names(out)[length(out)]<-toString(xmlGetAttr(criteriaComparisons[[i]],"mcdaConcept"))
      } # if (!("criteriaSet" %in% test1.names) & !("criteriaSet" %in% test2.names)){
    } # for (i in 1:length(criteriaComparisons)){
  }
  else {#if (length(criteriaComparisons)>0){
    err1<-"No <criteriaComparisons> found."
  }
  # In case there are <criteriaComparison> for sets of criteria, and none for
  # single criteria, out might be empty, but no error is detected. Therefore
  # we add this supplementary control. 
  if (length(out)==0)
    err1<-"No <criteriaComparisons> found."
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# Returns a list containing the comparisons of pairs of criteria.
# This restriction is done because we only use this structure for the comparison of interaction indexes
# on two criteria. 
# criteriaIDs is used to filter on which criteria ids the extraction is to be done.
# Possibility to specify which mcdaConcept should be searched. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getCriteriaPairsComparisons <- function(tree, criteriaIDs, mcdaConcept = NULL){
  
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <performanceTable> from the tree (according to the mcdaConcept if necessary)
  criteriaComparisons <- getNodeSet(tree, paste("//criteriaComparisons",specification,sep=""))
  
  # create the empty output list and the errors
  out<-list()
  err1<-NULL
  err2<-NULL
  
  
  if (length(criteriaComparisons)>0){
    for (i in 1:length(criteriaComparisons)){
      
      # check whether we only have <criteriaSet> under <initial> and <terminal>
      test1<-getNodeSet(criteriaComparisons[[i]], "pairs/pair/initial")
      test2<-getNodeSet(criteriaComparisons[[i]], "pairs/pair/terminal")
      test1.names<-NULL
      test2.names<-NULL
      tmpErr<-try(
{
  for (k in 1:length(test1))
    test1.names<-c(test1.names,names(xmlChildren(test1[[k]])))
  for (k in 1:length(test2))
    test2.names<-c(test2.names,names(xmlChildren(test2[[k]])))
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to read (a) value(s) in a <criteriaComparisons>."
      }
      
      if (!("criterionID" %in% test1.names) & !("criterionID" %in% test2.names)){
        
        criteriaComp <- matrix(nrow=0,ncol=5)
        
        pairs <- getNodeSet(criteriaComparisons[[i]], "pairs/pair")
        
        if (length(pairs)>0){
          for (j in 1:length(pairs)){
            
            head1<-NULL
            head2<-NULL
            tail1<-NULL
            tail2<-NULL
            val<-NULL
            
            initial <- getNodeSet(pairs[[j]],"initial")
            terminal <- getNodeSet(pairs[[j]],"terminal")
            
            
            
            if (names(xmlChildren(initial[[1]]))[1] == "criteriaSet"){
              # the comparisons are on sets of criteria (and we will consider here 
              # only sets of two elements, for the interaction)
              tmpErr<-try(
{	
  elements <- getNodeSet(pairs[[j]], "initial/criteriaSet/element")
  head1<- getNodeSet(elements[[1]], "criterionID")
  head2<- getNodeSet(elements[[2]], "criterionID")
}
              )
              if (inherits(tmpErr, 'try-error')){
                err2<-"Impossible to read (a) value(s) in a <criteriaComparisons>."
              }
              tmpErr<-try(
{	
  elements <- getNodeSet(pairs[[j]], "terminal/criteriaSet/element")
  tail1<- getNodeSet(elements[[1]], "criterionID")
  tail2<- getNodeSet(elements[[2]], "criterionID")
}
              )
              if (inherits(tmpErr, 'try-error')){
                err2<-"Impossible to read (a) value(s) in a <criteriaComparisons>."
              }
              
              
              
              tmpErr<-try(
{
  #										val <- getNodeSet(pairs[[j]], "value/real")
  #										if ((length(which(criteriaIDs==xmlValue(head1[[1]])))>0)&(length(which(criteriaIDs==xmlValue(head2[[1]])))>0)&(length(which(criteriaIDs==xmlValue(tail1[[1]])))>0)&(length(which(criteriaIDs==xmlValue(tail2[[1]])))>0))
  #											criteriaComp <-  rbind(criteriaComp,c(which(criteriaIDs==xmlValue(head1[[1]])),which(criteriaIDs==xmlValue(head2[[1]])),which(criteriaIDs==xmlValue(tail1[[1]])),which(criteriaIDs==xmlValue(tail2[[1]])),as.numeric(xmlValue(val[[1]]))))
  val <- getNodeSet(pairs[[j]], "value")
  if ((length(which(criteriaIDs==xmlValue(head1[[1]])))>0)&(length(which(criteriaIDs==xmlValue(head2[[1]])))>0)&(length(which(criteriaIDs==xmlValue(tail1[[1]])))>0)&(length(which(criteriaIDs==xmlValue(tail2[[1]])))>0))
    criteriaComp <-  rbind(criteriaComp,c(which(criteriaIDs==xmlValue(head1[[1]])),which(criteriaIDs==xmlValue(head2[[1]])),which(criteriaIDs==xmlValue(tail1[[1]])),which(criteriaIDs==xmlValue(tail2[[1]])),getNumericValue(val)))
}
              )
              if (inherits(tmpErr, 'try-error')){
                err2<-"Impossible to read (a) value(s) in a <criteriaComparisons>."
              }
            }
          }
        } #if (length(pairs)>0){
        if (dim(criteriaComp)[1] == 0)
          criteriaComp <- NULL
        out<-c(out,list(criteriaComp))
        names(out)[length(out)]<-toString(xmlGetAttr(criteriaComparisons[[i]],"mcdaConcept"))
      } # if (!("criterionID" %in% test1.names) & !("criterionID" %in% test2.names)){
    }
  }
  else {#if (length(criteriaComparisons)>0){
    err1<-"No <criteriaComparisons> found."
  }
  # In case there are <criteriaComparison> for single criteria, and none for
  # sets of criteria, out might be empty, but no error is detected. Therefore
  # we add this supplementary control. 
  if (length(out)==0)
    err1<-"No <criteriaComparisons> found."
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# Returns a list containing the comparisons of single alternatives. 
# Possibility to specify which mcdaConcept should be searched. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.
# The performance table must also be given as an argument because the matrix which is returned contains
# the evaluation of the alternatives. 

getAlternativesComparisons <- function(tree, performanceTable, mcdaConcept = NULL){
  
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <alternativesComparisons> from the tree (according to the mcdaConcept if necessary)
  alternativesComparisons <- getNodeSet(tree, paste("//alternativesComparisons",specification,sep=""))
  
  # create the empty output list and the errors
  out<-list()
  err1<-NULL
  err2<-NULL
  
  if (length(alternativesComparisons)>0){
    for (i in 1:length(alternativesComparisons)){
      alternativesComp <- matrix(nrow=0,ncol=(2*dim(performanceTable)[2]+1))
      pairs <- getNodeSet(alternativesComparisons[[i]], "pairs/pair")
      
      if (length(pairs)>0){
        for (j in 1:length(pairs)){
          
          head<-NULL
          tail<-NULL
          val<-NULL
          
          noPairs<-FALSE
          tmpErr<-try(
{
  head <- getNodeSet(pairs[[j]], "initial/alternativeID")
  tail <- getNodeSet(pairs[[j]], "terminal/alternativeID")
}
          )
          if (inherits(tmpErr, 'try-error')){
            err2<-"Impossible to read (a) value(s) in a <alternativesComparisons>."
            noPairs<-TRUE
          }
          
          tmpErr2<-try(
{
  #								val <- getNodeSet(pairs[[j]], "value/real")
  val <- getNodeSet(pairs[[j]], "value")
}
          )
          
          noVal<-FALSE
          if (inherits(tmpErr2, 'try-error')||(length(val)==0)){
            noVal<-TRUE
          }
          
          if ((noPairs == FALSE)&(noVal == FALSE)){
            #						val <- getNodeSet(pairs[[j]], "value/real")
            val <- getNodeSet(pairs[[j]], "value")
            try(
{
  #									alternativesComp <-  rbind(alternativesComp,c(performanceTable[xmlValue(head[[1]]),],performanceTable[xmlValue(tail[[1]]),],as.numeric(xmlValue(val[[1]]))))
  alternativesComp <-  rbind(alternativesComp,c(performanceTable[xmlValue(head[[1]]),],performanceTable[xmlValue(tail[[1]]),],getNumericValue(val)))
}
            )	
          }
          else if ((noPairs == FALSE)&(noVal == TRUE)){
            try(
{
  alternativesComp <-  rbind(alternativesComp,c(performanceTable[xmlValue(head[[1]]),],performanceTable[xmlValue(tail[[1]]),],NA))
}
            )
            
          }
          
        } # for (j in 1:length(pairs)){
      } # if (length(pairs)>0){
      if (dim(alternativesComp)[1] == 0)
        alternativesComp <- NULL
      out<-c(out,list(alternativesComp))
      names(out)[length(out)]<-toString(xmlGetAttr(alternativesComparisons[[i]],"mcdaConcept"))
    } # for (i in 1:length(alternativesComparisons)){
  } 
  else{ # if (length(alternativesComparisons)>0){
    err1<-"No <alternativesComparisons> found."
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# Returns a list of alternatives labels containing the comparisons of single alternatives. 
# Possibility to specify which mcdaConcept should be searched. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.
# The performance table must also be given as an argument because the matrix which is returned contains
# the evaluation of the alternatives. 

getAlternativesComparisonsLabels <- function(tree, altIDs=NULL, mcdaConcept = NULL){
  
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <alternativesComparisons> from the tree (according to the mcdaConcept if necessary)
  alternativesComparisons <- getNodeSet(tree, paste("//alternativesComparisons",specification,sep=""))
  
  # create the empty output list and the errors
  out<-list()
  err1<-NULL
  err2<-NULL
  
  if (length(alternativesComparisons)>0){
    for (i in 1:length(alternativesComparisons)){
      alternativesComp <- matrix(nrow=0,ncol=3)
      pairs <- getNodeSet(alternativesComparisons[[i]], "pairs/pair")
      
      if (length(pairs)>0){
        for (j in 1:length(pairs)){
          
          head<-NULL
          tail<-NULL
          val<-NULL
          
          noPairs<-FALSE
          tmpErr<-try(
{
  head <- getNodeSet(pairs[[j]], "initial/alternativeID")
  tail <- getNodeSet(pairs[[j]], "terminal/alternativeID")
}
          )
          if (inherits(tmpErr, 'try-error')){
            err2<-"Impossible to read (a) value(s) in a <alternativesComparisons>."
            noPairs<-TRUE
          }
          
          tmpErr2<-try(
{
  #								val <- getNodeSet(pairs[[j]], "value/real")
  val <- getNodeSet(pairs[[j]], "value")
}
          )
          
          noVal<-FALSE
          if (inherits(tmpErr2, 'try-error')||(length(val)==0)){
            noVal<-TRUE
          }
          
          if ((noPairs == FALSE)&(noVal == FALSE)){
            #						val <- getNodeSet(pairs[[j]], "value/real")
            val <- getNodeSet(pairs[[j]], "value")
            #						if (((xmlValue(head[[1]])%in%altIDs)&(xmlValue(tail[[1]])%in%altIDs))|(is.null(altIDs)))
            #							alternativesComp <-  rbind(alternativesComp,c(xmlValue(head[[1]]),xmlValue(tail[[1]]),as.numeric(xmlValue(val[[1]]))))	
            if (((xmlValue(head[[1]])%in%altIDs)&(xmlValue(tail[[1]])%in%altIDs))|(is.null(altIDs)))
              alternativesComp <-  rbind(alternativesComp,c(xmlValue(head[[1]]),xmlValue(tail[[1]]),getNumericValue(val)))
          }
          else if ((noPairs == FALSE)&(noVal == TRUE)){
            if (((xmlValue(head[[1]])%in%altIDs)&(xmlValue(tail[[1]])%in%altIDs))|(is.null(altIDs)))
              alternativesComp <-  rbind(alternativesComp,c(xmlValue(head[[1]]),xmlValue(tail[[1]]),NA))
          }
          
        } # for (j in 1:length(pairs))
      } # if (length(pairs)>0){
      if (dim(alternativesComp)[1] == 0)
        alternativesComp <- NULL
      out<-c(out,list(alternativesComp))
      names(out)[length(out)]<-toString(xmlGetAttr(alternativesComparisons[[i]],"mcdaConcept"))
    } # for (i in 1:length(alternativesComparisons)){
  } 
  else{ # if (length(alternativesComparisons)>0){
    err1<-"No <alternativesComparisons> found."
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

getAlternativesComparisonsValues <- function (tree, alternativesIDs = NULL, mcdaConcept = NULL) {
  specification = ""
  if (!is.null(mcdaConcept)) 
    specification <- paste("[@mcdaConcept='", mcdaConcept, "']", sep = "")
  alternativesComparisons <- getNodeSet(tree, paste("//alternativesComparisons", 
                                                    specification, sep = ""))
  out <- list()
  err1 <- NULL
  err2 <- NULL
  if (length(alternativesComparisons) > 0) {
    for (i in 1:length(alternativesComparisons)) {
      alternativesComp <- matrix(nrow = 0, ncol = 3)
      pairs <- getNodeSet(alternativesComparisons[[i]], "pairs/pair")
      if (length(pairs) > 0) {
        for (j in 1:length(pairs)) {
          head <- NULL
          tail <- NULL
          val <- NULL
          noPairs <- FALSE
          tmpErr <- try({
            head <- getNodeSet(pairs[[j]], "initial/alternativeID")
            tail <- getNodeSet(pairs[[j]], "terminal/alternativeID")
          })
          if (inherits(tmpErr, "try-error")) {
            err2 <- "Impossible to read (a) value(s) in a <alternativesComparisons>."
            noPairs <- TRUE
          }
          tmpErr2 <- try({
            val <- getNodeSet(pairs[[j]], "value")
          })
          noVal <- FALSE
          if (inherits(tmpErr2, "try-error") || (length(val) == 0)) {
            noVal <- TRUE
          }
          if ((noPairs == FALSE) & (noVal == FALSE)) {
            val <- getNodeSet(pairs[[j]], "value")
            if (((xmlValue(head[[1]]) %in% alternativesIDs) & 
                   (xmlValue(tail[[1]]) %in% alternativesIDs)) | (is.null(alternativesIDs))) 
              alternativesComp <- rbind(alternativesComp, 
                                        c(which(alternativesIDs == xmlValue(head[[1]])),
                                          which(alternativesIDs == xmlValue(tail[[1]])), 
                                          getNumericValue(val)))
          }
          else if ((noPairs == FALSE) & (noVal == TRUE)) {
            if (((xmlValue(head[[1]]) %in% alternativesIDs) & 
                   (xmlValue(tail[[1]]) %in% alternativesIDs)) | (is.null(alternativesIDs))) 
              alternativesComp <- rbind(alternativesComp, 
                                        c(which(alternativesIDs == xmlValue(head[[1]])),
                                          which(alternativesIDs == xmlValue(tail[[1]])), 
                                          NA))
          }
        }
      }
      if (dim(alternativesComp)[1] == 0) 
        alternativesComp <- NULL
      out <- c(out, list(alternativesComp))
      names(out)[length(out)] <- toString(xmlGetAttr(alternativesComparisons[[i]], 
                                                     "mcdaConcept"))
    }
  }
  else {
    err1 <- "No <alternativesComparisons> found."
  }
  if (!is.null(err1) | (!is.null(err2))) {
    out <- c(out, list(status = c(err1, err2)))
  }
  else {
    out <- c(out, list(status = "OK"))
  }
  return(out)
}

getCriteriaComparisonsLabels <- function(tree, critIDs=NULL, mcdaConcept = NULL){
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <criteriaComparisons> from the tree (according to the mcdaConcept if necessary)
  criteriaComparisons <- getNodeSet(tree, paste("//criteriaComparisons",specification,sep=""))
  
  # create the empty output list and the errors
  out<-list()
  err1<-NULL
  err2<-NULL
  
  if (length(criteriaComparisons)>0){
    for (i in 1:length(criteriaComparisons)){
      # check whether we only have <criteriaSet> under <initial> and <terminal>
      test1<-getNodeSet(criteriaComparisons[[i]], "pairs/pair/initial")
      test2<-getNodeSet(criteriaComparisons[[i]], "pairs/pair/terminal")
      test1.names<-NULL
      test2.names<-NULL
      tmpErr<-try(
{
  for (k in 1:length(test1))
    test1.names<-c(test1.names,names(xmlChildren(test1[[k]])))
  for (k in 1:length(test2))
    test2.names<-c(test2.names,names(xmlChildren(test2[[k]])))
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to read (a) value(s) in a <criteriaComparisons>."
      }
      
      if (!("criteriaSet" %in% test1.names) & !("criteriaSet" %in% test2.names)){
        criteriaComp <- matrix(nrow=0,ncol=3)
        pairs <- getNodeSet(criteriaComparisons[[i]], "pairs/pair")
        
        if (length(pairs)>0){
          for (j in 1:length(pairs)){
            
            head<-NULL
            tail<-NULL
            val<-NULL
            
            noPairs<-FALSE
            tmpErr<-try(
{
  head <- getNodeSet(pairs[[j]], "initial/criterionID")
  tail <- getNodeSet(pairs[[j]], "terminal/criterionID")
}
            )
            if (inherits(tmpErr, 'try-error')){
              err2<-"Impossible to read (a) value(s) in a <criteriaComparisons>."
              noPairs<-TRUE
            }
            
            tmpErr2<-try(
{
  #								val <- getNodeSet(pairs[[j]], "value/real")
  val <- getNodeSet(pairs[[j]], "value")
}
            )
            
            noVal<-FALSE
            if (inherits(tmpErr2, 'try-error')||(length(val)==0)){
              noVal<-TRUE
            }
            
            if ((noPairs == FALSE)&(noVal == FALSE)){
              #						val <- getNodeSet(pairs[[j]], "value/real")
              val <- getNodeSet(pairs[[j]], "value")
              #						if (((xmlValue(head[[1]])%in%critIDs)&(xmlValue(tail[[1]])%in%critIDs))|(is.null(critIDs)))
              #							criteriaComp <-  rbind(criteriaComp,c(xmlValue(head[[1]]),xmlValue(tail[[1]]),as.numeric(xmlValue(val[[1]]))))	
              if (((xmlValue(head[[1]])%in%critIDs)&(xmlValue(tail[[1]])%in%critIDs))|(is.null(critIDs)))
                criteriaComp <-  rbind(criteriaComp,c(xmlValue(head[[1]]),xmlValue(tail[[1]]),getNumericValue(val)))
            }
            else if ((noPairs == FALSE)&(noVal == TRUE)){
              if (((xmlValue(head[[1]])%in%critIDs)&(xmlValue(tail[[1]])%in%critIDs))|(is.null(critIDs)))
                criteriaComp <-  rbind(criteriaComp,c(xmlValue(head[[1]]),xmlValue(tail[[1]]),NA))
            }
            
          } # for (j in 1:length(pairs))
        } # if (length(pairs)>0){
        if (dim(criteriaComp)[1] == 0)
          criteriaComp <- NULL
        out<-c(out,list(criteriaComp))
        names(out)[length(out)]<-toString(xmlGetAttr(criteriaComparisons[[i]],"mcdaConcept"))
      }
    } # for (i in 1:length(criteriaComparisons)){
  } 
  else{ # if (length(criteriaComparisons)>0){
    err1<-"No <criteriaComparisons> found."
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}


# Returns a list containing the values of single criteria. 
# Possibility to specify which mcdaConcept should be searched. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getCriteriaValues <- function(tree, criteriaIDs, mcdaConcept = NULL){
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <performanceTable> from the tree (according to the mcdaConcept if necessary)
  criteriaValues <- getNodeSet(tree, paste("//criteriaValues",specification,sep=""))
  
  # create the empty output list and the errors
  out<-list()
  err1<-NULL
  err2<-NULL
  
  
  if (length(criteriaValues)>0){
    for (i in 1:length(criteriaValues)){
      
      # check whether we only have <criterionID> under <criterionValue>
      # and <real> under <value> (and not an interval)
      test1<-getNodeSet(criteriaValues[[i]], "criterionValue")
      test1.names<-NULL
      test2<-getNodeSet(criteriaValues[[i]], "criterionValue/value")
      test2.names<-NULL
      tmpErr<-try(
{
  for (k in 1:length(test1))
    test1.names<-c(test1.names,names(xmlChildren(test1[[k]])))
  for (k in 1:length(test2))
    test2.names<-c(test2.names,names(xmlChildren(test2[[k]])))
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to read (a) value(s) in a <criteriaValues>."
      }
      if (!("criteriaSet" %in% test1.names)& !("interval" %in% test2.names)){			
        criteriaVal <- matrix(nrow=0,ncol=2)
        
        vals <- getNodeSet(criteriaValues[[i]], "criterionValue")
        
        if (length(vals)>0){
          for (j in 1:length(vals)){
            tmpErr<-try(
{
  criterionID <- getNodeSet(vals[[j]], "criterionID")
  #									val <- getNodeSet(vals[[j]], "value/real")
  val <- getNodeSet(vals[[j]], "value")
  #									if (length(which(criteriaIDs==xmlValue(criterionID[[1]])))>0)
  #										criteriaVal <-rbind(criteriaVal,c(which(criteriaIDs==xmlValue(criterionID[[1]])),as.numeric(xmlValue(val[[1]]))))
  if (length(which(criteriaIDs==xmlValue(criterionID[[1]])))>0)
    criteriaVal <-rbind(criteriaVal,c(which(criteriaIDs==xmlValue(criterionID[[1]])),getNumericValue(val)))
}
            )
            if (inherits(tmpErr, 'try-error')){
              err2<-"Impossible to read (a) value(s) in a <criteriaValues>."
            }
          }
        } #if (length(vals)>0){
        if (dim(criteriaVal)[1] == 0)
          criteriaVal <- NULL
        out<-c(out,list(criteriaVal))
        names(out)[length(out)]<-toString(xmlGetAttr(criteriaValues[[i]],"mcdaConcept"))
      }
    }
  }
  else {#if (length(criteriaValues)>0){
    err1<-"No <criteriaValues> found."
  }
  # In case there are <criteriaValues> for sets of criteria, and none for
  # single criteria, out might be empty, but no error is detected. Therefore
  # we add this supplementary control. 
  if (length(out)==0)
    err1<-"No <criteriaValues> found."
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# Returns a list containing the intervals around the value of single criteria. 
# Possibility to specify which mcdaConcept should be searched. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getCriteriaIntervalValues <- function(tree, criteriaIDs, mcdaConcept = NULL){
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <criteriaValues> from the tree (according to the mcdaConcept if necessary)
  criteriaValues <- getNodeSet(tree, paste("//criteriaValues",specification,sep=""))
  
  # create the empty output list and the errors
  out<-list()
  err1<-NULL
  err2<-NULL
  
  
  if (length(criteriaValues)>0){
    for (i in 1:length(criteriaValues)){
      
      # check whether we only have <criterionID> under <criterionValue>
      # and <interval> under <value> (and not a real)
      test1<-getNodeSet(criteriaValues[[i]], "criterionValue")
      test1.names<-NULL
      test2<-getNodeSet(criteriaValues[[i]], "criterionValue/value")
      test2.names<-NULL
      tmpErr<-try(
{
  for (k in 1:length(test1))
    test1.names<-c(test1.names,names(xmlChildren(test1[[k]])))
  for (k in 1:length(test2))
    test2.names<-c(test2.names,names(xmlChildren(test2[[k]])))
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to read (a) value(s) in a <criteriaValues>."
      }
      #			if (!("criteriaSet" %in% test1.names)& !("real" %in% test2.names)){
      if (!("criteriaSet" %in% test1.names)& ("interval" %in% test2.names)){
        criteriaVal <- matrix(nrow=0,ncol=3)
        
        vals <- getNodeSet(criteriaValues[[i]], "criterionValue")
        
        if (length(vals)>0){
          for (j in 1:length(vals)){
            tmpErr<-try(
{
  criterionID <- getNodeSet(vals[[j]], "criterionID")
  #									lb <- getNodeSet(vals[[j]], "value/interval/lowerBound/real")
  #									ub <- getNodeSet(vals[[j]], "value/interval/upperBound/real")
  #									if (length(which(criteriaIDs==xmlValue(criterionID[[1]])))>0)
  #										criteriaVal <-rbind(criteriaVal,c(which(criteriaIDs==xmlValue(criterionID[[1]])),as.numeric(xmlValue(lb[[1]])),as.numeric(xmlValue(ub[[1]]))))
  lb <- getNodeSet(vals[[j]], "value/interval/lowerBound")
  ub <- getNodeSet(vals[[j]], "value/interval/upperBound")
  if (length(which(criteriaIDs==xmlValue(criterionID[[1]])))>0)
    criteriaVal <-rbind(criteriaVal,c(which(criteriaIDs==xmlValue(criterionID[[1]])),getNumericValue(lb),getNumericValue(ub)))						
}
            )
            if (inherits(tmpErr, 'try-error')){
              err2<-"Impossible to read (a) value(s) in a <criteriaValues>."
            }
          }
        } #if (length(vals)>0){
        if (dim(criteriaVal)[1] == 0)
          criteriaVal <- NULL
        out<-c(out,list(criteriaVal))
        names(out)[length(out)]<-toString(xmlGetAttr(criteriaValues[[i]],"mcdaConcept"))
      }
    }
  }
  else {#if (length(criteriaValues)>0){
    err1<-"No <criteriaValues> found."
  }
  # In case there are <criteriaValues> for sets of criteria, and none for
  # single criteria, out might be empty, but no error is detected. Therefore
  # we add this supplementary control. 
  if (length(out)==0)
    err1<-"No <criteriaValues> found."
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# Returns a list containing the values of sets of alternatives. 
# Possibility to specify which mcdaConcept should be searched. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getCriteriaPairsValues <- function(tree, criteriaIDs, mcdaConcept = NULL){
  
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <criteriaValues> from the tree (according to the mcdaConcept if necessary)
  criteriaValues <- getNodeSet(tree, paste("//criteriaValues",specification,sep=""))
  
  # create the empty output list and the errors
  out<-list()
  err1<-NULL
  err2<-NULL
  
  
  if (length(criteriaValues)>0){
    for (i in 1:length(criteriaValues)){
      
      # check whether we only have <criteriaSet> under <criterionValue>
      # and <real> under <value> (and not an interval)
      test1<-getNodeSet(criteriaValues[[i]], "criterionValue")
      test1.names<-NULL
      test2<-getNodeSet(criteriaValues[[i]], "criterionValue/value")
      test2.names<-NULL
      tmpErr<-try(
{
  for (k in 1:length(test1))
    test1.names<-c(test1.names,names(xmlChildren(test1[[k]])))
  for (k in 1:length(test2))
    test2.names<-c(test2.names,names(xmlChildren(test2[[k]])))
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to read (a) value(s) in a <criteriaValues>."
      }
      if (!("criterionID" %in% test1.names)& !("interval" %in% test2.names)){			
        criteriaVal <- matrix(nrow=0,ncol=3)
        
        vals <- getNodeSet(criteriaValues[[i]], "criterionValue")
        
        if (length(vals)>0){
          for (j in 1:length(vals)){
            
            tmpErr<-try(
{	
  elements <- getNodeSet(vals[[j]], "criteriaSet/element")
  head1<- getNodeSet(elements[[1]], "criterionID")
  head2<- getNodeSet(elements[[2]], "criterionID")
  #									val <- getNodeSet(vals[[j]], "value/real")
  val <- getNodeSet(vals[[j]], "value")
  #									if ((length(which(criteriaIDs==xmlValue(head1[[1]])))>0)&(length(which(criteriaIDs==xmlValue(head2[[1]])))>0))
  #										criteriaVal <-rbind(criteriaVal,c(which(criteriaIDs==xmlValue(head1[[1]])),which(criteriaIDs==xmlValue(head2[[1]])),as.numeric(xmlValue(val[[1]]))))
  if ((length(which(criteriaIDs==xmlValue(head1[[1]])))>0)&(length(which(criteriaIDs==xmlValue(head2[[1]])))>0))
    criteriaVal <-rbind(criteriaVal,c(which(criteriaIDs==xmlValue(head1[[1]])),which(criteriaIDs==xmlValue(head2[[1]])),getNumericValue(val)))								
}
            )
            if (inherits(tmpErr, 'try-error')){
              err2<-"Impossible to read (a) value(s) in a <criteriaValues>."
            }
          }
        } #if (length(vals)>0){
        if (dim(criteriaVal)[1] == 0)
          criteriaVal <- NULL
        out<-c(out,list(criteriaVal))
        names(out)[length(out)]<-toString(xmlGetAttr(criteriaValues[[i]],"mcdaConcept"))
      }
    }
  }
  else {#if (length(criteriaValues)>0){
    err1<-"No <criteriaValues> found."
  }
  # In case there are <criteriaValues> for sets of criteria, and none for
  # single criteria, out might be empty, but no error is detected. Therefore
  # we add this supplementary control. 
  if (length(out)==0)
    err1<-"No <criteriaValues> found."
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# Returns a list containing the intervals around the value of sets of criteria. 
# Possibility to specify which mcdaConcept should be searched. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getCriteriaPairsIntervalValues <- function(tree, criteriaIDs, mcdaConcept = NULL){
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <criteriaValues> from the tree (according to the mcdaConcept if necessary)
  criteriaValues <- getNodeSet(tree, paste("//criteriaValues",specification,sep=""))
  
  # create the empty output list and the errors
  out<-list()
  err1<-NULL
  err2<-NULL
  
  
  if (length(criteriaValues)>0){
    for (i in 1:length(criteriaValues)){
      
      # check whether we only have <criteriaSet> under <criterionValue>
      # and <interval> under <value> (and not a real)
      test1<-getNodeSet(criteriaValues[[i]], "criterionValue")
      test1.names<-NULL
      test2<-getNodeSet(criteriaValues[[i]], "criterionValue/value")
      test2.names<-NULL
      tmpErr<-try(
{
  for (k in 1:length(test1))
    test1.names<-c(test1.names,names(xmlChildren(test1[[k]])))
  for (k in 1:length(test2))
    test2.names<-c(test2.names,names(xmlChildren(test2[[k]])))
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to read (a) value(s) in a <criteriaValues>."
      }
      #			if (!("criterionID" %in% test1.names)& !("real" %in% test2.names)){
      if (!("criterionID" %in% test1.names) & ("interval" %in% test2.names)){
        
        criteriaVal <- matrix(nrow=0,ncol=4)
        
        vals <- getNodeSet(criteriaValues[[i]], "criterionValue")
        
        if (length(vals)>0){
          for (j in 1:length(vals)){
            
            tmpErr<-try(
{	
  elements <- getNodeSet(vals[[j]], "criteriaSet/element")
  head1<- getNodeSet(elements[[1]], "criterionID")
  head2<- getNodeSet(elements[[2]], "criterionID")
  #									lb <- getNodeSet(vals[[j]], "value/interval/lowerBound/real")
  #									ub <- getNodeSet(vals[[j]], "value/interval/upperBound/real")
  #									if((length(which(criteriaIDs==xmlValue(head1[[1]])))>0)&(length(which(criteriaIDs==xmlValue(head2[[1]])))>0))
  #										criteriaVal <-rbind(criteriaVal,c(which(criteriaIDs==xmlValue(head1[[1]])),which(criteriaIDs==xmlValue(head2[[1]])),as.numeric(xmlValue(lb[[1]])),as.numeric(xmlValue(ub[[1]]))))
  lb <- getNodeSet(vals[[j]], "value/interval/lowerBound")
  ub <- getNodeSet(vals[[j]], "value/interval/upperBound")
  if((length(which(criteriaIDs==xmlValue(head1[[1]])))>0)&(length(which(criteriaIDs==xmlValue(head2[[1]])))>0))
    criteriaVal <-rbind(criteriaVal,c(which(criteriaIDs==xmlValue(head1[[1]])),which(criteriaIDs==xmlValue(head2[[1]])),getNumericValue(lb),getNumericValue(ub)))
}
            )
            if (inherits(tmpErr, 'try-error')){
              err2<-"Impossible to read (a) value(s) in a <criteriaValues>."
            }
          }
        } #if (length(vals)>0){
        if (dim(criteriaVal)[1] == 0)
          criteriaVal <- NULL
        out<-c(out,list(criteriaVal))
        names(out)[length(out)]<-toString(xmlGetAttr(criteriaValues[[i]],"mcdaConcept"))
      }
    }
  }
  else {#if (length(criteriaValues)>0){
    err1<-"No <criteriaValues> found."
  }
  # In case there are <criteriaValues> for sets of criteria, and none for
  # single criteria, out might be empty, but no error is detected. Therefore
  # we add this supplementary control. 
  if (length(out)==0)
    err1<-"No <criteriaValues> found."
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# Returns a list containing the values of single alternatives. 
# Possibility to specify which mcdaConcept should be searched. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getAlternativesValues <- function(tree, alternativesIDs, mcdaConcept = NULL){
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <performanceTable> from the tree (according to the mcdaConcept if necessary)
  alternativesValues <- getNodeSet(tree, paste("//alternativesValues",specification,sep=""))
  
  # create the empty output list and the errors
  out<-list()
  err1<-NULL
  err2<-NULL
  
  
  if (length(alternativesValues)>0){
    for (i in 1:length(alternativesValues)){
      
      # check whether we only have <alternativeID> under <alternativeValue>
      # and <real> under <value> (and not an interval)
      test1<-getNodeSet(alternativesValues[[i]], "alternativeValue")
      test1.names<-NULL
      test2<-getNodeSet(alternativesValues[[i]], "alternativeValue/value")
      test2.names<-NULL
      tmpErr<-try(
{
  for (k in 1:length(test1))
    test1.names<-c(test1.names,names(xmlChildren(test1[[k]])))
  for (k in 1:length(test2))
    test2.names<-c(test2.names,names(xmlChildren(test2[[k]])))
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to read (a) value(s) in a <alternativesValues>."
      }
      if (!("alternativesSet" %in% test1.names)& !("interval" %in% test2.names)){			
        altVal <- matrix(nrow=0,ncol=2)
        
        vals <- getNodeSet(alternativesValues[[i]], "alternativeValue")
        
        if (length(vals)>0){
          for (j in 1:length(vals)){
            tmpErr<-try(
{
  alternativeID <- getNodeSet(vals[[j]], "alternativeID")
  #									val <- getNodeSet(vals[[j]], "value/real")
  #									if (length(which(alternativesIDs==xmlValue(alternativeID[[1]])))>0)
  #										altVal <-rbind(altVal,c(which(alternativesIDs==xmlValue(alternativeID[[1]])),as.numeric(xmlValue(val[[1]]))))
  val <- getNodeSet(vals[[j]], "value")
  if (length(which(alternativesIDs==xmlValue(alternativeID[[1]])))>0)
    altVal <-rbind(altVal,c(which(alternativesIDs==xmlValue(alternativeID[[1]])),getNumericValue(val)))
}
            )
            if (inherits(tmpErr, 'try-error')){
              err2<-"Impossible to read (a) value(s) in a <alternativesValues>."
            }
          }
        } #if (length(vals)>0){
        if (dim(altVal)[1] == 0)
          altVal <- NULL
        out<-c(out,list(altVal))
        names(out)[length(out)]<-toString(xmlGetAttr(alternativesValues[[i]],"mcdaConcept"))
      }
    }
  }
  else {#if (length(alternativesValues)>0){
    err1<-"No <alternativesValues> found."
  }
  # In case there are <alternativesValues> for sets of alternatives, and none for
  # single alternatives, out might be empty, but no error is detected. Therefore
  # we add this supplementary control. 
  if (length(out)==0)
    err1<-"No <alternativesValues> found."
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# Returns a list containing the intervals around the values of single alternatives. 
# Possibility to specify which mcdaConcept should be searched. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getAlternativesIntervalValues <- function(tree, alternativesIDs, mcdaConcept = NULL){
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <alternativesValues> from the tree (according to the mcdaConcept if necessary)
  alternativesValues <- getNodeSet(tree, paste("//alternativesValues",specification,sep=""))
  
  # create the empty output list and the errors
  out<-list()
  err1<-NULL
  err2<-NULL
  
  if (length(alternativesValues)>0){
    for (i in 1:length(alternativesValues)){
      
      # check whether we only have <alternativeID> under <alternativeValue>
      # and <interval> under <value> (and not a real)
      test1<-getNodeSet(alternativesValues[[i]], "alternativeValue")
      test1.names<-NULL
      test2<-getNodeSet(alternativesValues[[i]], "alternativeValue/value")
      test2.names<-NULL
      tmpErr<-try(
{
  for (k in 1:length(test1))
    test1.names<-c(test1.names,names(xmlChildren(test1[[k]])))
  for (k in 1:length(test2))
    test2.names<-c(test2.names,names(xmlChildren(test2[[k]])))
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to read (a) value(s) in a <alternativesValues>."
      }
      #			if (!("alternativesSet" %in% test1.names)& !("real" %in% test2.names)){
      if (!("alternativesSet" %in% test1.names)& ("interval" %in% test2.names)){
        criteriaVal <- matrix(nrow=0,ncol=3)
        
        vals <- getNodeSet(alternativesValues[[i]], "alternativeValue")
        
        if (length(vals)>0){
          for (j in 1:length(vals)){
            tmpErr<-try(
{
  alternativeID <- getNodeSet(vals[[j]], "alternativeID")
  #									lb <- getNodeSet(vals[[j]], "value/interval/lowerBound/real")
  #									ub <- getNodeSet(vals[[j]], "value/interval/upperBound/real")
  #									if (length(which(alternativesIDs==xmlValue(alternativeID[[1]])))>0)
  #										criteriaVal <-rbind(criteriaVal,c(which(alternativesIDs==xmlValue(alternativeID[[1]])),as.numeric(xmlValue(lb[[1]])),as.numeric(xmlValue(ub[[1]]))))
  lb <- getNodeSet(vals[[j]], "value/interval/lowerBound")
  ub <- getNodeSet(vals[[j]], "value/interval/upperBound")
  if (length(which(alternativesIDs==xmlValue(alternativeID[[1]])))>0)
    criteriaVal <-rbind(criteriaVal,c(which(alternativesIDs==xmlValue(alternativeID[[1]])),getNumericValue(lb),getNumericValue(ub)))
}
            )
            if (inherits(tmpErr, 'try-error')){
              err2<-"Impossible to read (a) value(s) in a <alternativesValues>."
            }
          }
        } #if (length(vals)>0){
        if (dim(criteriaVal)[1] == 0)
          criteriaVal <- NULL
        out<-c(out,list(criteriaVal))
        names(out)[length(out)]<-toString(xmlGetAttr(alternativesValues[[i]],"mcdaConcept"))
      }
    }
  }
  else {#if (length(alternativesValues)>0){
    err1<-"No <alternativesValues> found."
  }
  # In case there are <alternativesValues> for sets of alternatives, and none for
  # single alternatives, out might be empty, but no error is detected. Therefore
  # we add this supplementary control. 
  if (length(out)==0)
    err1<-"No <alternativesValues> found."
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

getAlternativesAffectations <- function (tree, alternativesIDs, categoriesIDs,
                                         mcdaConcept = NULL)  {
  specification = ""
  if (!is.null(mcdaConcept))
    specification <- paste("[@mcdaConcept='", mcdaConcept, "']", sep = "")
  alternativesAffectations <- getNodeSet(tree, paste("//alternativesAffectations", 
                                                     specification, sep = ""))
  out <- list()
  err1 <- NULL
  err2 <- NULL
  
  if (length(alternativesAffectations) > 0) {
    for (i in 1:length(alternativesAffectations)) {
      nodes <- getNodeSet(alternativesAffectations[[i]], "alternativeAffectation")
      altAff <- matrix(data = FALSE, nrow = length(alternativesIDs), ncol = length(categoriesIDs))
      
      if (length(nodes) > 0) {
        for (j in 1:length(nodes)) {
          tmpErr <- try({
            alternativeID <- getNodeSet(nodes[[j]], "alternativeID")
            categoryID <- getNodeSet(nodes[[j]], "categoryID")
            categoriesInterval <- getNodeSet(nodes[[j]], "categoriesInterval")
            categoriesSet <- getNodeSet(nodes[[j]], "categoriesSet")
            
            if (length(alternativeID) == 0) {
              err1 <- "Missing <alternativeID> in <alternativeAffectation>."
            }
            else {
              if (length(categoryID) == 1) {
                categoryIndex = which(categoriesIDs == xmlValue(categoryID[[1]]))
                if (length(categoryIndex) > 0) {
                  altAff[which(alternativesIDs == xmlValue(alternativeID[[1]])), categoryIndex] = TRUE
                }
              }
              else if (length(categoriesSet) == 1) {
                elements <- getNodeSet(categoriesSet[[1]], "element/categoryID")
                for (element in elements) {
                  categoryIndex = which(categoriesIDs == xmlValue(element[[1]]))
                  if (length (categoryIndex) > 0) {
                    altAff[which(alternativesIDs == xmlValue(alternativeID[[1]])), categoryIndex] = TRUE
                  }
                }
              }
              else if (length(categoriesInterval) == 1) {
                lowerBound <- getNodeSet(categoriesInterval[[1]], "lowerBound/categoryID")
                upperBound <- getNodeSet(categoriesInterval[[1]], "upperBound/categoryID")
                
                lowerBoundIndex <- 1
                upperBoundIndex <- ncol(altAff)
                
                if (length(lowerBound) > 0) {
                  lowerBoundIndex <- which(categoriesIDs == xmlValue(lowerBound[[1]]))
                  if (length(lowerBoundIndex) == 0) lowerBoundIndex <- 1
                }
                
                if (length(upperBound) > 0) {
                  upperBoundIndex <- which(categoriesIDs == xmlValue(upperBound[[1]]))
                  if (length(upperBoundIndex) == 0) upperBoundIndex <- ncol(altAff)
                }
                
                for (k in 1:ncol(altAff)) {
                  if (k >= lowerBoundIndex && k <= upperBoundIndex) {
                    altAff[which(alternativesIDs == xmlValue(alternativeID[[1]])), k] = TRUE
                  }
                }
              }
            }        
          })
          if (inherits(tmpErr, "try-error")) {
            err2 <- "Impossible to read (a) value(s) in a <alternativesAffectations>."
          }
        }
        
        out <- c(out, list(altAff))
        names(out)[length(out)] <- toString(xmlGetAttr(alternativesAffectations[[i]], 
                                                       "mcdaConcept"))
      }
    }
  }
  else {
    err1 <- "No <alternativesAffectations> found."
  }
  if (length(out) == 0) 
    err1 <- "No <alternativesAffectations> found."
  if (!is.null(err1) | (!is.null(err2))) {
    out <- c(out, list(status = c(err1, err2)))
  }
  else {
    out <- c(out, list(status = "OK"))
  }
  return(out)
}

getCategoriesValues <- function(tree, categoriesIDs, mcdaConcept = NULL) {
  specification = ""
  if (!is.null(mcdaConcept))
    specification <- paste("[@mcdaConcept='", mcdaConcept, "']", sep = "")	
  
  categoriesValues <- getNodeSet(tree, paste("//categoriesValues",
                                             specification, sep = ""))
  
  out<-list()
  err1<-NULL
  err2<-NULL
  
  if (length(categoriesValues)>0){
    for (i in 1:length(categoriesValues)){
      # check whether we only have <categoryID> under <categoryValue>
      # and <real> under <value> (and not an interval)
      test1<-getNodeSet(categoriesValues[[i]], "categoryValue")
      test1.names<-NULL
      test2<-getNodeSet(categoriesValues[[i]], "categoryValue/value")
      test2.names<-NULL
      tmpErr<-try(
{
  for (k in seq_len(length(test1)))
    test1.names<-c(test1.names, names(xmlChildren(test1[[k]])))
  for (k in seq_len(length(test2)))
    test2.names<-c(test2.names, names(xmlChildren(test2[[k]])))
})
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to read (a) value(s) in a <categoryValue>."
      }
      if (!("categoriesSet" %in% test1.names)& !("interval" %in% test2.names)){			
        catVal <- matrix(nrow=0,ncol=2)
        
        vals <- getNodeSet(categoriesValues[[i]], "categoryValue")
        
        if (length(vals)>0){
          for (j in 1:length(vals)){
            tmpErr<-try(
{
  categoryID <- getNodeSet(vals[[j]], "categoryID")
  val <- getNodeSet(vals[[j]], "value")
  if (length(which(categoriesIDs == xmlValue(categoryID[[1]]))) > 0)
    catVal <- rbind(catVal, c(which(categoriesIDs == xmlValue(categoryID[[1]])),
                              getNumericValue(val)))
})
            if (inherits(tmpErr, 'try-error')){
              err2<-"Impossible to read (a) value(s) in a <categoriesValues>."
            }
          }
        }
        
        if (dim(catVal)[1] == 0)
          catVal <- NULL
        out<-c(out,list(catVal))
        names(out)[length(out)] <- toString(xmlGetAttr(categoriesValues[[i]],
                                                       "mcdaConcept"))
      }
    }
  }
  else {
    err1<-"No <categoriesValues> found."
  }
  
  if (length(out) == 0)
    err1<-"No <categoriesValues> found."
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

getCategoriesIntervalValues <- function (tree, categoriesIDs, mcdaConcept = NULL) {
  specification = ""
  if (!is.null(mcdaConcept)) 
    specification <- paste("[@mcdaConcept='", mcdaConcept, "']", sep = "")
  categoriesValues <- getNodeSet(tree, paste("//categoriesValues", 
                                             specification, sep = ""))
  out <- list()
  err1 <- NULL
  err2 <- NULL
  if (length(categoriesValues) > 0) {
    for (i in 1:length(categoriesValues)) {
      test1 <- getNodeSet(categoriesValues[[i]], "categoryValue")
      test1.names <- NULL
      test2 <- getNodeSet(categoriesValues[[i]], "categoryValue/value")
      test2.names <- NULL
      tmpErr <- try({
        for (k in seq_len(length(test1))) test1.names <- c(test1.names, 
                                                           names(xmlChildren(test1[[k]])))
        for (k in seq_len(length(test2))) test2.names <- c(test2.names, 
                                                           names(xmlChildren(test2[[k]])))
      })
      if (inherits(tmpErr, "try-error")) {
        err2 <- "Impossible to read (a) value(s) in a <categoriesValues>."
      }
      
      if (!("categoriesSet" %in% test1.names) & ("interval" %in% test2.names)) {
        catVal <- matrix(nrow = 0, ncol = 3)
        vals <- getNodeSet(categoriesValues[[i]], "categoryValue")
        
        for (j in seq_len(length(vals))) {
          tmpErr <- try({
            categoryID <- getNodeSet(vals[[j]], "categoryID")
            lowerBoundNode <- getNodeSet(vals[[j]], "value/interval/lowerBound")
            upperBoundNode <- getNodeSet(vals[[j]], "value/interval/upperBound")
            
            lowerBound <- NA
            upperBound <- NA
            
            if (length(lowerBoundNode) == 1)
              lowerBound <- getNumericValue(lowerBoundNode)
            
            if (length(upperBoundNode) == 1)
              upperBound <- getNumericValue(upperBoundNode)
            
            if (length(which(categoriesIDs == xmlValue(categoryID[[1]]))) > 0 &&
                  (!is.na(lowerBound) || !is.na(upperBound))) 
              catVal <- rbind(catVal, c(which(categoriesIDs == xmlValue(categoryID[[1]])),
                                        lowerBound,
                                        upperBound))
          })
          if (inherits(tmpErr, "try-error")) {
            err2 <- "Impossible to read (a) value(s) in a <categoriesValues>."
          }
        }
        
        if (dim(catVal)[1] == 0) 
          catVal <- NULL
        out <- c(out, list(catVal))
        names(out)[length(out)] <- toString(xmlGetAttr(categoriesValues[[i]], 
                                                       "mcdaConcept"))
      }
    }
  }
  else {
    err1 <- "No <categoriesValues> found."
  }
  if (length(out) == 0) 
    err1 <- "No <categoriesValues> found."
  if (!is.null(err1) | (!is.null(err2))) {
    out <- c(out, list(status = c(err1, err2)))
  }
  else {
    out <- c(out, list(status = "OK"))
  }
  return(out)
}

# Puts criteriaValues in the XML tree.
# Returns an error if something went wrong.  
# Possibility to specify which mcdaConcept should be written. 

putCriteriaValues <- function(tree, criteriaValues, criteriaIDs, mcdaConcept = NULL){
  
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No root tag found."
  }
  
  if (length(racine)!=0){
    
    
    if (!is.null(mcdaConcept)){
      critVals<-newXMLNode("criteriaValues", attrs = c(mcdaConcept=mcdaConcept), parent=racine, namespace=c())
      
    }
    else{
      
      critVals<-newXMLNode("criteriaValues", parent=racine, namespace=c())
      
    }
    for (i in 1:dim(criteriaValues)[1]){
      tmpErr<-try(
{
  critVal<-newXMLNode("criterionValue", parent=critVals, namespace=c())
  newXMLNode("criterionID", criteriaIDs[criteriaValues[i,1]], parent = critVal, namespace=c())
  val<-newXMLNode("value", parent = critVal, namespace=c())
  newXMLNode("real",criteriaValues[i,2], parent=val, namespace=c())
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to put (a) value(s) in a <criteriaValues>."
      }
    }	
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# Puts criteriaValues in the XML tree.
# Returns an error if something went wrong.  
# Possibility to specify which mcdaConcept should be written. 

putCriteriaMatrix <- function(tree, criteriaMatrix, mcdaConcept = NULL){
  
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No root tag found."
  }
  
  if (length(racine)!=0){
    
    
    if (!is.null(mcdaConcept)){
      criteriaMat<-newXMLNode("criteriaMatrix", attrs = c(mcdaConcept=mcdaConcept), parent=racine, namespace=c())
      
    }
    else{
      
      criteriaMat<-newXMLNode("criteriaMatrix", parent=racine, namespace=c())
      
    }
    
    tmpErr<-try(
{
  for (i in 1:dim(criteriaMatrix)[1]){
    row <-newXMLNode("row", parent=criteriaMat, namespace=c())
    newXMLNode("criterionID", rownames(criteriaMatrix)[i], parent = row, namespace=c())
    for (j in 1:dim(criteriaMatrix)[2]){
      col <-newXMLNode("column", parent=row, namespace=c())
      newXMLNode("criterionID", colnames(criteriaMatrix)[j], parent = col, namespace=c())
      val<-newXMLNode("value", parent = col, namespace=c())
      newXMLNode("real",criteriaMatrix[i,j], parent=val, namespace=c())
    }
  }
}
    )
    if (inherits(tmpErr, 'try-error')){
      err2<-"Impossible to put (a) value(s) in a <criteriaMatrix>."
    }
    
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

putCriterionValue <- function(tree, criterionValue, criteriaIDs = NULL, mcdaConcept = NULL){
  
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    
    if (!is.null(mcdaConcept)){
      critVals<-newXMLNode("criterionValue", attrs = c(mcdaConcept=mcdaConcept), parent=racine, namespace=c())
      
    }
    else{
      
      critVals<-newXMLNode("criterionValue", parent=racine, namespace=c())
    }
    tmpErr<-try(
{
  if (!is.null(criteriaIDs)){
    if (length(criteriaIDs) == 1){
      newXMLNode("criterionID", criteriaIDs, parent = critVals, namespace=c())
    }
    else{
      critSet<-newXMLNode("criteriaSet", parent=critVals, namespace=c())
      for (i in 1:length(criteriaIDs)){
        element<-newXMLNode("element", parent=critSet, namespace=c())
        newXMLNode("criterionID", criteriaIDs[i], parent = element, namespace=c())
      }
    }
  }
  val<-newXMLNode("value", parent = critVals, namespace=c())
  newXMLNode("real",criterionValue, parent=val, namespace=c())
}
    )
    if (inherits(tmpErr, 'try-error')){
      err2<-"Impossible to put (a) value(s) in a <criteriaValues>."
    }	
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

putCriteriaPlot <- function(tree, base64plot, criteriaIDs, mcdaConcept=NULL, name=NULL){
  
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    if ((is.null(mcdaConcept))&(is.null(name)))
    {
      critVal<-newXMLNode("criterionValue", parent=racine, namespace=c())
    }
    else if (is.null(mcdaConcept))
    {
      critVal<-newXMLNode("criterionValue", attrs = c(name=name), parent=racine, namespace=c())
    }
    else if (is.null(name))
    {
      critVal<-newXMLNode("criterionValue", attrs = c(mcdaConcept=mcdaConcept), parent=racine, namespace=c())
    }
    else
    {
      critVal<-newXMLNode("criterionValue", attrs = c(mcdaConcept=mcdaConcept, name=name), parent=racine, namespace=c())
    }
    
    #		if (is.null(mcdaConcept))
    #			critVal<-newXMLNode("criterionValue", attrs = c(name=name), parent=racine, namespace=c())
    #		else
    #			critVal<-newXMLNode("criterionValue", attrs = c(mcdaConcept=mcdaConcept, name=name), parent=racine, namespace=c())
    
    tmpErr<-try(
{
  critSet<-newXMLNode("criteriaSet", parent=critVal, namespace=c())
  for (i in 1:length(criteriaIDs)){
    element<-newXMLNode("element", parent=critSet, namespace=c())
    newXMLNode("criterionID", criteriaIDs[i], parent = element, namespace=c())
  }
  val<-newXMLNode("value", parent = critVal, namespace=c())
  newXMLNode("image",base64plot, parent=val, namespace=c())
}
    )
    if (inherits(tmpErr, 'try-error')){
      err2<-"Impossible."
    }
    
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

putAlternativesPlot <- function(tree, base64plot, alternativesIDs, mcdaConcept=NULL, name=NULL){
  
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    if ((is.null(mcdaConcept))&(is.null(name)))
    {
      critVal<-newXMLNode("alternativeValue", parent=racine, namespace=c())
    }
    else if (is.null(mcdaConcept))
    {
      critVal<-newXMLNode("alternativeValue", attrs = c(name=name), parent=racine, namespace=c())
    }
    else if (is.null(name))
    {
      critVal<-newXMLNode("alternativeValue", attrs = c(mcdaConcept=mcdaConcept), parent=racine, namespace=c())
    }
    else
    {
      critVal<-newXMLNode("alternativeValue", attrs = c(mcdaConcept=mcdaConcept, name=name), parent=racine, namespace=c())
    }
    
    
    #		if (is.null(mcdaConcept)){
    #			critVal<-newXMLNode("alternativeValue", attrs = c(name=name), parent=racine, namespace=c())
    #		}
    #		else
    #		{
    #			critVal<-newXMLNode("alternativeValue", attrs = c(mcdaConcept=mcdaConcept, name=name), parent=racine, namespace=c())	
    #		}
    
    
    tmpErr<-try(
{
  critSet<-newXMLNode("alternativesSet", parent=critVal, namespace=c())
  for (i in 1:length(alternativesIDs)){
    element<-newXMLNode("element", parent=critSet, namespace=c())
    newXMLNode("alternativeID", alternativesIDs[i], parent = element, namespace=c())
  }
  val<-newXMLNode("value", parent = critVal, namespace=c())
  newXMLNode("image",base64plot, parent=val, namespace=c())
}
    )
    if (inherits(tmpErr, 'try-error')){
      err2<-"Impossible."
    }
    
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

putCriteriaPairsValues <- function(tree, criteriaPairsValues, criteriaIDs, mcdaConcept = NULL){
  
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    
    if (!is.null(mcdaConcept)){
      critVals<-newXMLNode("criteriaValues", attrs = c(mcdaConcept=mcdaConcept), parent=racine, namespace=c())
      
    }
    else{
      
      critVals<-newXMLNode("criteriaValues", parent=racine, namespace=c())
      
    }
    for (i in 1:dim(criteriaPairsValues)[1]){
      tmpErr<-try(
{
  critVal<-newXMLNode("criterionValue", parent=critVals, namespace=c())
  critSet<-newXMLNode("criteriaSet", parent=critVal, namespace=c())
  element<-newXMLNode("element", parent=critSet, namespace=c())
  newXMLNode("criterionID", criteriaIDs[criteriaPairsValues[i,1]], parent = element, namespace=c())
  element<-newXMLNode("element", parent=critSet, namespace=c())
  newXMLNode("criterionID", criteriaIDs[criteriaPairsValues[i,2]], parent = element, namespace=c())
  val<-newXMLNode("value", parent = critVal, namespace=c())
  newXMLNode("real",criteriaPairsValues[i,3], parent=val, namespace=c())
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to put (a) value(s) in a <criteriaValues>."
      }
    }	
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
  
  
}

putAlternativesValues <- function(tree, alternativesValues, alternativesIDs, mcdaConcept = NULL){
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    
    if (!is.null(mcdaConcept)){
      altVals<-newXMLNode("alternativesValues", attrs = c(mcdaConcept=mcdaConcept), parent=racine, namespace=c())
      
    }
    else{
      
      altVals<-newXMLNode("alternativesValues", parent=racine, namespace=c())
      
    }
    for (i in 1:dim(alternativesValues)[1]){
      tmpErr<-try(
{
  altVal<-newXMLNode("alternativeValue", parent=altVals, namespace=c())
  newXMLNode("alternativeID", alternativesIDs[alternativesValues[i,1]], parent = altVal, namespace=c())
  val<-newXMLNode("value", parent = altVal, namespace=c())
  newXMLNode("real",alternativesValues[i,2], parent=val, namespace=c())
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to put (a) value(s) in a <alternativesValues>."
      }
    }	
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

putAlternativesIDs <- function(tree, alternativesIDs, mcdaConcept = NULL){
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    
    if (!is.null(mcdaConcept)){
      altIDs<-newXMLNode("alternatives", attrs = c(mcdaConcept=mcdaConcept), parent=racine, namespace=c())
      
    }
    else{
      
      altIDs<-newXMLNode("alternatives", parent=racine, namespace=c())
      
    }
    for (i in 1:length(alternativesIDs)){
      tmpErr<-try(
{
  altID<-newXMLNode("alternative", attrs=c(id=alternativesIDs[i]),parent=altIDs, namespace=c())
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to put an id in <alternative>."
      }
    }	
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}


putAlternativeValue <- function(tree, alternativeValue, alternativesIDs=NULL, mcdaConcept = NULL){
  
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    
    if (!is.null(mcdaConcept)){
      critVals<-newXMLNode("alternativeValue", attrs = c(mcdaConcept=mcdaConcept), parent=racine, namespace=c())
      
    }
    else{
      
      critVals<-newXMLNode("alternativeValue", parent=racine, namespace=c())
    }
    tmpErr<-try(
{
  if (!is.null(alternativesIDs)){
    if (length(alternativesIDs) == 1){
      newXMLNode("alternativeID", alternativesIDs, parent = critVals, namespace=c())
    }
    else{
      critSet<-newXMLNode("alternativesSet", parent=critVals, namespace=c())
      for (i in 1:length(alternativesIDs)){
        element<-newXMLNode("element", parent=critSet, namespace=c())
        newXMLNode("alternativeID", alternativesIDs[i], parent = element, namespace=c())
      }
    }
  }
  val<-newXMLNode("value", parent = critVals, namespace=c())
  newXMLNode("real",alternativeValue, parent=val, namespace=c())
}
    )
    if (inherits(tmpErr, 'try-error')){
      err2<-"Impossible to put (a) value(s) in a <alternativeValue>."
    }	
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

putAlternativesComparisonsLabels <-function(tree, alternativesComparisons, mcdaConcept = NULL){
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    
    if (!is.null(mcdaConcept)){
      altVals<-newXMLNode("alternativesComparisons", attrs = c(mcdaConcept=mcdaConcept), parent=racine, namespace=c())
      
    }
    else{
      
      altVals<-newXMLNode("alternativesComparisons", parent=racine, namespace=c())
      
    }
    
    pairs<-newXMLNode("pairs", parent=altVals, namespace=c())
    
    for (i in 1:dim(alternativesComparisons)[1]){
      tmpErr<-try(
{
  pair<-newXMLNode("pair", parent=pairs, namespace=c())
  initial<-newXMLNode("initial", parent=pair, namespace=c())
  newXMLNode("alternativeID", alternativesComparisons[i,1], parent = initial, namespace=c())
  terminal<-newXMLNode("terminal", parent=pair, namespace=c())
  newXMLNode("alternativeID", alternativesComparisons[i,2], parent = terminal, namespace=c())
  if (dim(alternativesComparisons)[2] > 2)
  {
    val<-newXMLNode("value", parent = pair, namespace=c())
    if (is.na(alternativesComparisons[i,3])){
      newXMLNode("NA",alternativesComparisons[i,3], parent=val, namespace=c())
    }
    else
      newXMLNode("real",alternativesComparisons[i,3], parent=val, namespace=c())
  }
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to put (a) value(s) in a <alternativesComparisons>."
      }
    }	
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

putAlternativesAffectations <- function (tree, alternativesAffectations, 
                                         alternativesIDs, categoriesIDs,
                                         asIntervalsIfPossible = FALSE,
                                         mcdaConcept = NULL) {
  out <- list()
  err1 <- NULL
  err2 <- NULL
  root <- NULL
  tmpErr <- try({
    root <- xmlRoot(tree)
  })
  if (inherits(tmpErr, "try-error")) {
    err1 <- "No <xmcda:XMCDA> found."
  }
  if (length(root) != 0) {
    if (!is.null(mcdaConcept)) {
      alternativesAffectationsNode <- newXMLNode("alternativesAffectations",
                                                 attrs = c(mcdaConcept = mcdaConcept), 
                                                 parent = root,
                                                 namespace = c())
    }
    else {
      alternativesAffectationsNode <- newXMLNode("alternativesAffectations",
                                                 parent = root,
                                                 namespace = c())
    }
    
    for (i in 1:nrow(alternativesAffectations)) {
      tmpErr <- try({
        altAffectation <- newXMLNode("alternativeAffectation",
                                     parent = alternativesAffectationsNode,
                                     namespace = c())
        newXMLNode("alternativeID", alternativesIDs[i],
                   parent = altAffectation, namespace = c())
        
        trueIndices = which (alternativesAffectations[i, ] == TRUE)
        
        if (length(trueIndices) != 0) {
          firstTrue = min(trueIndices)
          lastTrue = max(trueIndices)
          if (asIntervalsIfPossible == TRUE &&
                (lastTrue - firstTrue + 1) == length(trueIndices)) {
            interval <- newXMLNode("categoriesInterval", parent = altAffectation, namespace = c())
            lowerBound <- newXMLNode("lowerBound", parent = interval, namespace = c())
            upperBound <- newXMLNode("upperBound", parent = interval, namespace = c())
            newXMLNode("categoryID", categoriesIDs[firstTrue], parent = lowerBound, namespace = c())
            newXMLNode("categoryID", categoriesIDs[lastTrue], parent = upperBound, namespace = c())
          } else {
            categoriesSet <- newXMLNode("categoriesSet", parent = altAffectation, namespace = c())
            for (j in trueIndices) {
              newXMLNode("categoryID", categoriesIDs[j], parent = categoriesSet, namespace = c())
            }
          }
        } 
      })
      if (inherits(tmpErr, "try-error")) {
        err2 <- "Impossible to put (a) value(s) in a <alternativesAffectations>."
        break
      }
    }
  }
  if (!is.null(err1) | (!is.null(err2))) {
    out <- c(out, list(status = c(err1, err2)))
  }
  else {
    out <- c(out, list(status = "OK"))
  }
  return(out)
}

putAlternativesAffectationsWithValues <- function (tree, alternativesAffectations, 
                                         alternativesIDs, categoriesIDs,
                                         mcdaConcept = NULL) {
  out <- list()
  err1 <- NULL
  err2 <- NULL
  root <- NULL
  tmpErr <- try({
    root <- xmlRoot(tree)
  })
  if (inherits(tmpErr, "try-error")) {
    err1 <- "No <xmcda:XMCDA> found."
  }
  if (length(root) != 0) {
    if (!is.null(mcdaConcept)) {
      alternativesAffectationsNode <- newXMLNode("alternativesAffectations",
                                                 attrs = c(mcdaConcept = mcdaConcept), 
                                                 parent = root,
                                                 namespace = c())
    }
    else {
      alternativesAffectationsNode <- newXMLNode("alternativesAffectations",
                                                 parent = root,
                                                 namespace = c())
    }
    
    for (i in 1:nrow(alternativesAffectations)) {
      alternativeIndex <- alternativesAffectations[i, 1]
      categoryIndex <- alternativesAffectations[i, 2]
      value <- alternativesAffectations[i, 3]
      
      tmpErr <- try({
        altAffectation <- newXMLNode("alternativeAffectation",
                                     parent = alternativesAffectationsNode,
                                     namespace = c())
        newXMLNode("alternativeID", alternativesIDs[alternativeIndex],
                   parent = altAffectation, namespace = c())
        newXMLNode("categoryID", categoriesIDs[categoryIndex],
                   parent = altAffectation, namespace = c())
        valueNode <- newXMLNode("value", parent = altAffectation, namespace = c())
        newXMLNode("real", value, parent = valueNode, namespace = c())
      })
      if (inherits(tmpErr, "try-error")) {
        err2 <- "Impossible to put (a) value(s) in a <alternativesAffectations>."
        break
      }
    }
  }
  if (!is.null(err1) | (!is.null(err2))) {
    out <- c(out, list(status = c(err1, err2)))
  }
  else {
    out <- c(out, list(status = "OK"))
  }
  return(out)
}

putCategoriesValues <- function (tree, categoriesValues, categoriesIDs,
                                 mcdaConcept = NULL) {
  out <- list()
  err1 <- NULL
  err2 <- NULL
  root <- NULL
  tmpErr <- try({
    root <- xmlRoot(tree)
  })
  if (inherits(tmpErr, "try-error")) {
    err1 <- "No <xmcda:XMCDA> found."
  }
  if (length(root) != 0) {
    if (!is.null(mcdaConcept)) {
      categoriesValuesParent <- newXMLNode("categoriesValues",
                                           attrs = c(mcdaConcept = mcdaConcept), 
                                           parent = root, namespace = c())
    }
    else {
      categoriesValuesParent <- newXMLNode("categoriesValues", parent = root, 
                                           namespace = c())
    }
    if (nrow(categoriesValues) > 0) {
      for (i in 1:dim(categoriesValues)[1]) {
        tmpErr <- try({
          categoryValue <- newXMLNode("categoryValue",
                                      parent = categoriesValuesParent,
                                      namespace = c())
          newXMLNode("categoryID", categoriesIDs[categoriesValues[i, 1]],
                     parent = categoryValue, namespace = c())
          if (ncol(categoriesValues) == 2) {
            value <- newXMLNode("value", parent = categoryValue, namespace = c())
            newXMLNode("real", categoriesValues[i, 2], parent = value, namespace = c())
          } else if (ncol(categoriesValues) > 2) {
            values <- newXMLNode("values", parent = categoryValue, namespace = c())

            for (j in seq(2, ncol(categoriesValues))) {
              value <- newXMLNode("value", parent = values, namespace = c())
              newXMLNode("real", categoriesValues[i, j], parent = value, namespace = c())
            }
          }
        })
        if (inherits(tmpErr, "try-error")) {
          err2 <- "Impossible to put (a) value(s) in a <categoriesValues>."
        }
      }
    }
  }
  if (!is.null(err1) | (!is.null(err2))) {
    out <- c(out, list(status = c(err1, err2)))
  }
  else {
    out <- c(out, list(status = "OK"))
  }
  return(out)
}

putCategoriesIntervalValues <- function (tree, categoriesValues, categoriesIDs,
                                         mcdaConcept = NULL) {
  out <- list()
  err1 <- NULL
  err2 <- NULL
  root <- NULL
  tmpErr <- try({
    root <- xmlRoot(tree)
  })
  if (inherits(tmpErr, "try-error")) {
    err1 <- "No <xmcda:XMCDA> found."
  }
  if (length(root) != 0) {
    if (!is.null(mcdaConcept)) {
      categoriesValuesParent <- newXMLNode("categoriesValues",
                                           attrs = c(mcdaConcept = mcdaConcept), 
                                           parent = root, namespace = c())
    }
    else {
      categoriesValuesParent <- newXMLNode("categoriesValues", parent = root, 
                                           namespace = c())
    }
    if (nrow(categoriesValues) > 0) {
      for (i in 1:dim(categoriesValues)[1]) {
        tmpErr <- try({
          if (!is.na(categoriesValues[i, 2]) || !is.na(categoriesValues[i, 3])) {
            categoryValue <- newXMLNode("categoryValue",
                                        parent = categoriesValuesParent,
                                        namespace = c())
            newXMLNode("categoryID", categoriesIDs[categoriesValues[i, 1]],
                       parent = categoryValue, namespace = c())
            value <- newXMLNode("value", parent = categoryValue, namespace = c())
            interval <- newXMLNode("interval", parent = value, namespace = c())
            if (!is.na(categoriesValues[i, 2])) {
              lowerBound <- newXMLNode("lowerBound", parent = interval)
              newXMLNode("real", categoriesValues[i, 2], parent = lowerBound)
            }
            if (!is.na(categoriesValues[i, 3])) {
              upperBound <- newXMLNode("upperBound", parent = interval)
              newXMLNode("real", categoriesValues[i, 3], parent = upperBound)
            }
          }
        })
        if (inherits(tmpErr, "try-error")) {
          err2 <- "Impossible to put (a) value(s) in a <alternativesValues>."
        }
      }
    }
  }
  if (!is.null(err1) | (!is.null(err2))) {
    out <- c(out, list(status = c(err1, err2)))
  }
  else {
    out <- c(out, list(status = "OK"))
  }
  return(out)
}

putErrorMessage <- function(tree, errorMessage, name = NULL){
  
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    # we now check if <methodMessages> exists
    # if not, we have to add it
    
    methMessages<-getNodeSet(tree, "//methodMessages")
    
    if (length(methMessages)==0){
      # if no <methodMessages> can be found, create it
      methMessages<-newXMLNode("methodMessages", parent=racine, namespace=c())
      if (!is.null(name)){
        methMessage<-newXMLNode("errorMessage", attrs = c(name=name), parent=methMessages, namespace=c())
      }
      else{
        methMessage<-newXMLNode("errorMessage", parent=methMessages, namespace=c())
      }
    }
    else
    {
      if (!is.null(name)){
        methMessage<-newXMLNode("errorMessage", attrs = c(name=name), parent=methMessages[[1]], namespace=c())
      }
      else{
        methMessage<-newXMLNode("errorMessage", parent=methMessages[[1]], namespace=c())
      }
    }
    newXMLNode("text", errorMessage, parent = methMessage, namespace=c())		
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

putLogMessage <- function(tree, logMessage, name = NULL){
  
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    # we now check if <methodMessages> exists
    # if not, we have to add it
    
    methMessages<-getNodeSet(tree, "//methodMessages")
    
    if (length(methMessages)==0){
      # if no <methodMessages> can be found, create it
      methMessages<-newXMLNode("methodMessages", parent=racine, namespace=c())
      if (!is.null(name)){
        lMessage<-newXMLNode("logMessage", attrs = c(name=name), parent=methMessages, namespace=c())
      }
      else{
        lMessage<-newXMLNode("logMessage", parent=methMessages, namespace=c())
      }
    }
    else
    {
      if (!is.null(name)){
        lMessage<-newXMLNode("logMessage", attrs = c(name=name), parent=methMessages[[1]], namespace=c())
      }
      else{
        lMessage<-newXMLNode("logMessage", parent=methMessages[[1]], namespace=c())
      }
    }
    newXMLNode("text", logMessage, parent = lMessage, namespace=c())		
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}


putMessage <- function(tree, message, name = NULL){
  
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    # we now check if <methodMessages> exists
    # if not, we have to add it
    
    methMessages<-getNodeSet(tree, "//methodMessages")
    
    if (length(methMessages)==0){
      # if no <methodMessages> can be found, create it
      methMessages<-newXMLNode("methodMessages", parent=racine, namespace=c())
      if (!is.null(name)){
        lMessage<-newXMLNode("message", attrs = c(name=name), parent=methMessages, namespace=c())
      }
      else{
        lMessage<-newXMLNode("message", parent=methMessages, namespace=c())
      }
    }
    else
    {
      if (!is.null(name)){
        lMessage<-newXMLNode("message", attrs = c(name=name), parent=methMessages[[1]], namespace=c())
      }
      else{
        lMessage<-newXMLNode("message", parent=methMessages[[1]], namespace=c())
      }
    }
    newXMLNode("text", message, parent = lMessage, namespace=c())		
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

# Returns a list containing the capacities on the criteria. 
# Possibility to specify which mcdaConcept should be searched. 
# The elements of the list are named according to the mcdaConcept attribute, if it has been defined.

getMobiusCapacities<-function(tree, criteriaIDs, numberOfCriteria, kadditivity, mcdaConcept = NULL){
  
  # if an mcdaConcept has been specified, search according to this attribute
  specification = ""
  if (!is.null(mcdaConcept)) specification <- paste("[@mcdaConcept='",mcdaConcept,"']",sep="")	
  
  # extract the <criteriaValues> from the tree (according to the mcdaConcept if necessary)
  criteriaValues <- getNodeSet(tree, paste("//criteriaValues",specification,sep=""))
  
  # create the empty output list and the errors
  out<-list()
  err1<-NULL
  err2<-NULL
  err3<-NULL
  err4<-NULL
  m<-NULL
  
  
  if (length(criteriaValues)>0){
    for (i in 1:length(criteriaValues)){
      
      # check whether we only have <criteriaSet> under <criterionValue>
      # and <real> under <value> (and not an interval)
      test1<-getNodeSet(criteriaValues[[i]], "criterionValue")
      test1.names<-NULL
      test2<-getNodeSet(criteriaValues[[i]], "criterionValue/value")
      test2.names<-NULL
      tmpErr<-try(
{
  for (k in 1:length(test1))
    test1.names<-c(test1.names,names(xmlChildren(test1[[k]])))
  for (k in 1:length(test2))
    test2.names<-c(test2.names,names(xmlChildren(test2[[k]])))
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to read (a) value(s) in a <criteriaValues>."
      }
      if (!("criterionID" %in% test1.names)& !("interval" %in% test2.names)){			
        criteriaVal <- NULL
        
        vals <- getNodeSet(criteriaValues[[i]], "criterionValue")
        
        if (length(vals)>0){
          
          # calculation of the number of elements of the capacity
          num<-0
          for (k in 1:kadditivity){
            num<-num+choose(numberOfCriteria,k)
          }
          
          # if two few criteria are active, don't read the capacity
          
          if (num==length(vals)){
            
            cap<-c()
            for (j in 1:length(vals)){
              # read the values and the subsets from the tree
              # the subsets are transformed into their code
              # which corresponds to the data in m@subsets
              tmpErr<-try(
{	
  elements <- getNodeSet(vals[[j]], "criteriaSet/element")
  tmp<-c()
  code<-0
  for (k in 1:length(elements)){
    critID<-getNodeSet(elements[[k]], "criterionID")
    code<-code+2^(which(criteriaIDs==xmlValue(critID[[1]]))-1)
  }
  
  #										val <- getNodeSet(vals[[j]], "value/real")
  val <- getNodeSet(vals[[j]], "value")
  if (length(val)==0)
    tmp<-c(code,0)
  else{
    #											tmp<-c(code,as.numeric(xmlValue(val[[1]])))
    tmp<-c(code,getNumericValue(val))
  }
  cap<-rbind(cap,tmp)
}
              )
              if (inherits(tmpErr, 'try-error')){
                err2<-"Impossible to read (a) value(s) in a <criteriaValues>."
              }
            }
            
            # creation of the Mobius capacity
            tmpErr<-try(
{	
  m<-Mobius.capacity(rep(0,num+1),numberOfCriteria,kadditivity)
}
            )	
            if (inherits(tmpErr, 'try-error')){
              err3<-"Impossible to create a Mobius capacity."
            }
            # fill in the values read from the tree
            tmpErr<-try(
{for (k in 1:dim(cap)[1]){	
  m@data[which(m@subsets==cap[k,1])]<-cap[k,2]
}
}
            )	
            if (inherits(tmpErr, 'try-error')){
              err4<-"Errors in filling the Mobius capacity."
            }
          }# if (num==length(vals)){
          else {
            err4 <- "Wrong number of criteria active for capacity to be read."
          }
        } #if (length(vals)>0){
        out<-c(out,list(m))
        names(out)[length(out)]<-toString(xmlGetAttr(criteriaValues[[i]],"mcdaConcept"))
      }
    }
  }
  else {#if (length(criteriaValues)>0){
    err1<-"No <criteriaValues> found."
  }
  # In case there are <criteriaValues> for sets of criteria, and none for
  # single criteria, out might be empty, but no error is detected. Therefore
  # we add this supplementary control. 
  if (length(out)==0)
    err1<-"No <criteriaValues> found."
  if (!is.null(err1)|(!is.null(err2))|(!is.null(err3))|(!is.null(err4))){
    out<-c(out,list(status=c(err1,err2,err3,err4)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}


putCapacity<- function(tree, capacity, criteriaIDs, mcdaConcept = NULL){
  
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    
    if (!is.null(mcdaConcept)){
      critVals<-newXMLNode("criteriaValues", attrs = c(mcdaConcept=mcdaConcept), parent=racine, namespace=c())
      
    }
    else{
      
      critVals<-newXMLNode("criteriaValues", parent=racine, namespace=c())
      
    }
    
    c<-to.data.frame(capacity)
    
    # generate a list containing all the possible subsets of criteria in a "natural" order
    
    subsets<-list(c(),c(1))
    
    compteur <- 2
    
    for (i in 2:length(criteriaIDs)){
      for (j in 1:compteur){
        compteur<-compteur + 1
        subsets[[compteur]]<-c(subsets[[j]],i)
      }
    }
    
    # then create the xml file
    
    for (i in 2:length(capacity@subsets)){
      tmpErr<-try(
{
  critVal = newXMLNode("criterionValue", parent = critVals, namespace=c())
  critSet = newXMLNode("criteriaSet", parent = critVal, namespace=c())
  tmp<-subsets[[capacity@subsets[i]+1]]
  for (j in 1:length(tmp)){
    elt<-newXMLNode("element",parent=critSet, namespace=c())
    newXMLNode("criterionID", criteriaIDs[tmp[j]],parent=elt, namespace=c())
  }
  val = newXMLNode("value", parent = critVal, namespace=c())
  newXMLNode("real",c[i+2,1], parent=val, namespace=c())
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to put (a) value(s) in a <criteriaValues>."
      }
    }
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

putPerformanceTable <- function(tree, performanceTable, mcdaConcept = NULL){
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    
    if (!is.null(mcdaConcept)){
      perfTab<-newXMLNode("performanceTable", attrs = c(mcdaConcept=mcdaConcept), parent=racine, namespace=c())
      
    }
    else{
      
      perfTab<-newXMLNode("performanceTable", parent=racine, namespace=c())
      
    }
    for (i in 1:dim(performanceTable)[1]){
      tmpErr<-try(
{
  altPerf<-newXMLNode("alternativePerformances", parent=perfTab, namespace=c())
  newXMLNode("alternativeID", rownames(performanceTable)[i], parent=altPerf, namespace=c())
  for (j in 1:dim(performanceTable)[2]){
    perf<-newXMLNode("performance", parent=altPerf, namespace=c())
    newXMLNode("criterionID", colnames(performanceTable)[j], parent=perf, namespace=c())
    val<-newXMLNode("value", parent=perf, namespace=c())
    newXMLNode("real", performanceTable[i,j], parent=val, namespace=c())
  }
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to put (a) value(s) in a <performanceTable>."
      }
    }	
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}

putPointsCriterionFunction <- function(tree, points, mcdaConcept = NULL){
  out<-list()
  err1<-NULL
  err2<-NULL
  racine<-NULL
  
  # we first check if there is an <xmcda:XMCDA>
  
  tmpErr<-try(
{
  racine<-xmlRoot(tree)
}
  )
  if (inherits(tmpErr, 'try-error')){
    err1<-"No <xmcda:XMCDA> found."
  }
  
  if (length(racine)!=0){
    
    
    if (!is.null(mcdaConcept)){
      criteria<-newXMLNode("criteria", attrs = c(mcdaConcept=mcdaConcept), parent=racine, namespace=c())
      
    }
    else{
      
      criteria<-newXMLNode("criteria", parent=racine, namespace=c())
      
    }
    for (i in 1:length(points)){
      tmpErr<-try(
{
  criterion <- newXMLNode("criterion", attrs=c(id=names(points)[i]), parent=criteria, namespace=c())
  criterionFunction <- newXMLNode("criterionFunction", parent=criterion, namespace=c())
  pts <- newXMLNode("points", parent=criterionFunction, namespace=c())
  for (j in 1:dim(points[[i]])[1]){
    pt <- newXMLNode("point", parent=pts, namespace=c())
    abs <- newXMLNode("abscissa", parent=pt, namespace=c())
    newXMLNode("real", points[[i]][j,1], parent=abs, namespace=c())
    ord <- newXMLNode("ordinate", parent=pt, namespace=c())
    newXMLNode("real", points[[i]][j,2], parent=ord, namespace=c())
  }
}
      )
      if (inherits(tmpErr, 'try-error')){
        err2<-"Impossible to put (a) value(s) in a <criterionFunction>."
      }
    }	
  }
  if (!is.null(err1)|(!is.null(err2))){
    out<-c(out,list(status=c(err1,err2)))
  }
  else{
    out<-c(out,list(status="OK"))
  }
  return(out)
}



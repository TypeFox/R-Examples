#' Tools for interacting with 'jMetrik'
#' 
#' \code{jmetrik} provides tools for using R and the external program 'jMetrik' togther. In particular, it
#' provides tools for reading and writing files in *.jmetrik format. This format is required by 'jMetrik'
#' and it is a flat file with a header and comma separated values. However, a *.jmetrik file is not simply
#' a CSV file. The header includes much more information about the data than is typically found in a CSV file.
#' 'jMetrik' is an open source Java application for psychometric analysis. It may be downloaded from
#' \url{http://www.ItemAnalysis.com}.
#' 
#' See \code{\link{jmetrikWrite}} and \url{http://www.ItemAnalysis.com} for more information about the *.jmetrik
#' file format.
#'  
#' 
#' @docType package
#' @name jmetrik
NULL


#' Writes a file in *.jmetrik format. 
#' 
#' 'jMetrik' is a stand alone program written in Java. It defines a file format that is just a plain text
#' file with a header and comma delimited values. The header contains information about the variables
#' in the file. There is one row in the header for each variable in the file. The header also includes
#' meta information about the data such as the number of rows. This function will create a *.jmetrik
#' file from a data frame. The jMetrik program and other informaiton is available at
#' \url{http://www.Itemanalysis.com}
#' 
#' @param x A data frame
#' @param fileName The complete path and name of the file to be written. The file siffix must be .jmetrik.
#' @param scoring An optional character vector of item scoring. Each element in this vector has two sets 
#'  of parentheses. The first set contains the response option codes. The second set contains the scores 
#'  assigned to each option. There is a correspondence between each set of parentheses such that the first 
#'  element in the code list corresponds to the first element in the score list.
#'  
#' @param codes An optional character vector of special codes. Each element in this vector has two sets 
#'  of parentheses. The first set contains the missing data, omitted, and not reached codes. 
#'  The second set contains the scores assigned to each code. There is a correspondence between 
#'  each set of parentheses such that the first element in the code list corresponds to the first 
#'  element in the score list.
#'  
#' @param labels An optional character vector of variable labels
#' @author J. Patrick Meyer \email{support@@itemanalysis.com}
#' @examples
#' \dontrun{
#'
#' #Create some data
#' id<-100+seq(1:10)
#' x<-sample(c("A", "B", "C", "D"), 10, replace=TRUE)
#' y<-sample(c("A", "B", "C", "D"), 10, replace=TRUE)
#' z<-sample(c(0,1,2,3), 10, replace=TRUE)
#' sc<-rnorm(10)
#' exdata<-as.data.frame(cbind(id, x, y, z, sc))
#' names(exdata)<-c("id", "item1", "item2", "item3", "score")
#'
#' #A is the correct answer
#' aOK<- "(A,B,C,D)(1,0,0,0)"
#'
#' #B is the correct answer
#' bOK<-"(A,B,C,D)(0,1,0,0)"
#' 
#' #polytomous item scoring
#' poly<-"(0,1,2,3)(0,1,2,3)"
#' 
#' #Special data codes e.g. missing and not reached responses
#' #These can be unique to each item or the same. Here they
#' #are the same.
#' datCodes<-"(NA,OM,NR)(0,0,0)"
#'
#' #Create scoring, special data codes, and labels
#' scoring<-c("", aOK, bOK, poly, "")
#' codes<-c("", rep(datCodes, 3), "")
#' labels<-c("ID variable", "Test item 1", "Test item 2", "Test item 3", "Test score")
#' 
#' #write the file
#' jmetrikWrite(x=exdata, fileName="./test-write.jmetrik",
#'             codes=codes, scoring=scoring, labels=labels)
#'             
#' }
#' 
#' @export    
jmetrikWrite<-function(x, fileName, scoring=NULL, codes=NULL, labels=NULL){
    if(!is.data.frame(x)) stop("x must be a data frame.")
    file.create(fileName)
    conn<-file(fileName, "w")
    
    vNames <- names(x)
    
    #Write start of header
    cat("# VERSION", "jmetrik1", "# METADATA", nrow(x), "# ATTRIBUTES", file=conn, sep="\n")
    
    #Write variable attributes to the header
    for(i in 1:length(vNames)){
        myname<-vNames[i]
        
        #Variable type
        if(length(dim(x))==2){
            if(is.numeric(x[,i])){
                mytype<-"DOUBLE"
            }else{
                mytype<-"STRING"
            } 
        }else{
            if(is.numeric(x)){
                mytype<-"DOUBLE"
            }else{
                mytype<-"STRING"
            } 
        }        
               
        #Item scoring information
        if(!is.null(scoring) & length(scoring)>1){
            myscoring<-paste("\"", scoring[i], "\"", sep="") 
        }else{
            myscoring<-""
        }
        
        #Special code information
        if(!is.null(codes) & length(codes)>1){
            mycodes<-paste("\"", codes[i], "\"", sep="")
        }else{
            mycodes<-""
        }
        
        #Variable labels
        if(!is.null(labels) & length(labels)>1){
            mylabel<-paste("\"", labels[i], "\"", sep="")
        }else{
            mylabel<-""
        }
        
        #write line
        cat(myname, mytype, myscoring, mycodes, mylabel, file=conn, sep=",")
        cat("", file=conn, sep="\n")
    }
    
    cat("# DATA", file=conn, sep="\n")
    write.table(x=x, file=fileName, append=TRUE, row.names=FALSE, col.names=FALSE, qmethod="double", sep=",")

    close(conn)
    unlink(conn)
   
}

#' Reads a *.jmetrik file into a data frame. 
#' 
#' A *.jmetrik file can be created with \code{\link{jmetrikWrite}} or by the 'jMetrik' program. 
#' See \url{http://www.ItemAnalysis.com}.
#' 
#' @param fileName The complete path and file name of the *.jmetrik file that is being read.
#' @return a data frame
#' @examples
#' \dontrun{
#' x<-jmetrikRead(fileName="./myfile.jmetrik")
#' }
#' 
#' @export
jmetrikRead<-function(fileName){
    conn<-file(fileName)
    firstRow<-0
    variableAttributes<-FALSE
    scanning<-TRUE
    vname<-vector(mode = "character", length = 0)
    i<-1
    
    while(length(oneLine <- readLines(conn, n = 1, warn = FALSE)) > 0 & scanning){
        if("# ATTRIBUTES"==oneLine){
            variableAttributes<-TRUE            
        }else if("# DATA"==oneLine){
            scanning<-FALSE
            variableAttributes<-FALSE
        }
        
        if(variableAttributes==TRUE){
            myVector<-(strsplit(oneLine, ","))
            vname[i]<-myVector[1]
            #TODO get other attributes here
        }
    
        firstRow<-firstRow+1
    }#end while
    close(conn)
    
    
    (x<-read.csv(file=fileName, col.names=vname, skip=firstRow))
}


#' WoS exported citation report for search AD=(brigham same anes*) OR AD=(brigham same anaes*)
#' @name BWHCitationReport
#' @docType data
#' @keywords data
NULL
#' H-Index Calculator using Data from a Web of Science (WoS) Citation Report
#' @param readcsv - The read.csv function that will locate the citation report of interest
#' @param hx - The h-index of interest. Can be any positive integer or "c" (without quotations) for cumulative h-index
#' @param date - The last year of interest plus 1. For instance, if the citation report includes the year 2015 but you want h(10) not incuding 2015, therefore 2005-2014, you would enter 2015 for this argument
#' @examples
#'  #calculate the h(10) of the Brigham and Women's Hospital - Department of Anesthesia 
#'  data(BWHCitationReport)
#'  readcsv <- BWHCitationReport
#'  hx <- 10
#'  date <- 2015
#'
#'  hindexcalculator(readcsv, hx, date)
#' @export
#BWHCitationReport
# 

hindexcalculator <- function (readcsv, hx, date) {
  
  data <-readcsv #read.csv("Ref 1.csv")
  
  date <- date #2015
  
  hx <- hx #h(hx) where hx is the number of years being analyzed; c = cumulative
  
  data1 <-data[28:nrow(data),]
  colnames(data1)[1:21] <- c("Title", "Authors", "Corporate Authors", "Editors", "Book Editors", "Source Title", "Publication Date", "Publication Year", "Volume", "Issue","Part Number", "Supplement", "Special Issue", "Begining Page", "Ending Page", "Article Number", "Doi", "Conference Title", "Conference Date", "Total Citations", "Average per year")
  colnames(data1)[22:ncol(data)] <- c(1945:2015)
  pubyear <- data1[8] 
  
  data2 <- cbind(pubyear,data1[,22:ncol(data)])
  c <- ncol(data2)-1
  
  if (hx==c){
    exclude <- 1945
  } else { exclude <- date-hx}
  
  data2[,1] <- as.numeric(as.character(data2[,1]))
  data2[,1]
  data3 <- data2[,1] == date |  data2[,1] < exclude #publications prior to exclude date and in current year
  #data2[data3,] #publications prior to exclude date and in current year
  data4<-data2[!data3,] #excludes publications prior to exclude date and in current year
  
  nyears <- ncol(data4)-1
  
  
  if (hx == c){
    upperlimit<- nyears  #to set the upper limit as the latest year the spreadsheet has: upperlimit <- nyears+1 #2015
    lowerlimit <- (upperlimit-(hx-2)) 
  } else {
    upperlimit<- nyears
    lowerlimit <- (upperlimit-(hx-1))   
  }
  
  firstpublication <- min(data4[,1])
  
  data5<- data4[,lowerlimit:upperlimit]
  data5
  sum <- apply(data5,1,sum)
  
  hindex<-cbind( data5, sum)
  
  pubnum<-nrow(hindex)
  pubnum
  
  totalcitations<-sum(sum)
  totalcitations
  
  hindex <- hindex[order(hindex[,"sum"], decreasing = TRUE),] #sorts the data from highest sum of citations to lowest
  rank <- 1:nrow(hindex) #ranks the sorted data from 1 (highest citation)  
  hindex<-cbind(hindex,rank) #combines the rank list to the table
  
  pubnum<-nrow(hindex) #the number of publications 
  
  
  for(i in 1:pubnum) {
    if (pubnum==0){
      output <- 0
    }
  }
  
  for(i in 1:pubnum) { #goes through each publication
    
    if(hindex[i,"rank"] > hindex[i,"sum"]) {
      #when the program reaches the first rank that is higher than the sum of citations, it sets the rank of the publication directly above it as the h-index of that reference
      y <- hindex[i-1,"rank"] # y is the proposed h-index 
      if (length(y)==0){ #ADDING THIS
        y<-0 
      }
      break
    }
    
  }
  
  for(i in 1:pubnum) { #to account for cases where the rank is equal to the sum of citations 
    if(hindex[i,"rank"]==hindex[i,"sum"]) {
      #the first publication that has a rank equal to the sum of citations
      z <- hindex[i,"rank"] #propose the rank of the first publication that has a rank equal to the sum of citation as the hindex 
      break
    } else { #if no such a thing exist assign a value of zero to this loop
      z<- 0   
    }
  }
  
  if (z>0) { #in the case that there is a case where rank=sum (the first instance)
    output <- z #assign that rank as the h-index, anything larger won't work, anything smaller is accounted for
  } else { #if no equal cases, 
    output <- y
  }
  
  #if (z < 1 && y <1 ){
  #  output <- 0
  #}
  
  if (pubnum==1){ #if there is one publication
    if (hindex[1,"sum"] >=1){ #and the sum of citations for that one publication is at least one
      output<- 1 #the h-index is one
    }
    if (hindex[1,"sum"] == 0) { #if there are no publications
      output <- 0 #the h-index is 0
    }
  }
  
  matrix<-matrix(,1,4)
  if (hx==c){
    colnames(matrix) <- c("Year of First Publication", "H-c", "Total # of Publications", "Sum of the Times Cited")
  }else{
    colnames(matrix) <- c("Year of First Publication", paste("H-", hx,"", sep=""), "Total # of Publications", "Sum of the Times Cited")
  }
  matrix[1,1]<- firstpublication
  matrix[1,2]<- output
  matrix[1,3]<- pubnum
  matrix[1,4]<- totalcitations
  class(matrix)<-"numeric"
  
  #library(xlsx)
  #write.xlsx(matrix, "hindexcalculator-output.xlsx")
  
  return(matrix)
  
}

# Calculates the h(x) - the h-index for the past x years - of a scientist/department/etc. using the exported excel file from a Web of Science (WoS) citation report of a search. Also calculated is the year of first publication, total number of publications, and sum of times cited for the specified period (e.g. for h-10, date of first publication, total number of publications, and sum of times cited in the past 10 years). The exported WoS citation report excel file has to first be saved in a .csv format.
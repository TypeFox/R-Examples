#### Read the text from the web site for the book and create 'TsayFiles'
#### list and subdirectory.
####
####

##
## 1.  Web site
##
FinTS.url <- "http://faculty.chicagogsb.edu/ruey.tsay/teaching/fts2"

localDir <- "~/TsayFiles/"

Tsay.webPage <- scan(FinTS.url, what=character(0), sep="\n")
n.webLines <- length(Tsay.webPage)

##
## 2.  Chapters
##
(Ch <- paste("Chapter ", 1:12, ":", sep=""))

(twelve <- c(paste("0", 1:9, sep=""), 10:12))
(ch <- paste("ch", twelve, sep=""))
names(Ch) <- ch

ch. <- rep(NA, 12)
names(ch.) <- ch
for(i in 1:12)
  ch.[i] <- grep(Ch[i], Tsay.webPage)

Tsay.webPage[ch.]
# ch. = line numbers of chapter breaks

##
## 3.  Exercises 
##
(ex <- grep("Exercise", Tsay.webPage))

chEnd <- c(ch.[-1]-1, n.webLines)
names(chEnd) <- ch
ex. <- chEnd
i <- 1
for(ich in 1:12){
  if(ex[i] < chEnd[ich]){
    ex.[ich] <- ex[i]
    i <- i+1
  }
}

(ChEx <- cbind(ch=ch., ex=ex., end=chEnd))

##
## 4.  Files
##

# All files to be downloaded must have a web address 
Tsay.http <- grep("http://", Tsay.webPage)
Tsay.http. <- Tsay.webPage[Tsay.http]
# but not all web addresses are files

# First delete leading junk 
http.start <- regexpr("http://", Tsay.http.)
Tsay.ht1 <- substring(Tsay.http., http.start)
http.stop <- regexpr("\"", Tsay.ht1)
Tsay.ht2 <- substring(Tsay.ht1, 1, http.stop-1)

Tsay.ht2. <- strsplit(Tsay.ht2[-1], "/")
Tsay.ht3 <- sapply(Tsay.ht2., function(x)x[length(x)])

Tsay.ht4 <- Tsay.ht2[-(1:10)]
names(Tsay.ht4) <- Tsay.ht3[-(1:9)]
Tsay.ht4. <- Tsay.http[-(1:10)]

last4 <- substring(Tsay.ht4, nchar(Tsay.ht4)-3)
table(last4)
#.dat .sor .txt at.f ou.f rats 
#  44    1  120    1    1   14 

fileTypes <- c("dat", "txt", "sor", "rats", "f") 

kTypes <- length(fileTypes)
filesByType <- vector("list", length=kTypes)
names(filesByType) <- fileTypes

for(i in 1:kTypes){
  fti <- paste(".", fileTypes[i], "$", sep="")
  filesByType[[i]] <- grep(fti, Tsay.ht4)
}
(fileLineNums <- sort(unlist(filesByType)))

fLNm <- Tsay.http[-(1:10)][fileLineNums]

##
## 5.  TsayFiles
##

TsayFiles <- vector("list", length=12)
names(TsayFiles) <- ch

for(ich in 1:12){
  sel <- ((ChEx[ich, 1]<=fLNm) & (fLNm<=ChEx[ich, 2]))
#  
  TsayFiles[[ich]]$text <- url2dat(Tsay.ht4[sel]) 
#
  sel2 <- ((ChEx[ich, 2]<=fLNm) & (fLNm<=ChEx[ich, 3]))
#
  TsayFiles[[ich]]$exercises <- url2dat(Tsay.ht4[sel2])
}

# check
readSum <- array(0, dim=c(12, 2, 2), dimnames=list(
       ch, c("text", "exercises"), c("FALSE", "TRUE") ) )
for(ich in 1:12)for(j in 1:2){
  TFij <- TsayFiles[[ich]][[j]]
  {
    if(dim(TFij)[1]>0){
      tij <- table(TFij[, 4]) 
      readSum[ich, j, names(tij)] <- tij
    }
    else
      readSum[ich, j, ] <- 0 
  }
}
# Numbers of URLs found without corresponding files:  
#     text exercises
#ch01    3         2
#ch02    1         2
#ch03    0         4
#ch04    5         1
#ch05    0         1
#ch06    1         0
#ch07    2         1
#ch08    1         0
#ch09    0         0
#ch10    0         1
#ch11    0         0
#ch12    4         1

# Numbers of URLs downloaded 
#     text exercises
#ch01   18         6
#ch02   10        11
#ch03    8         5
#ch04    9         4
#ch05   14         7
#ch06    3         0
#ch07    6         4
#ch08   10         3
#ch09    6         2
#ch10   11         3
#ch11    1         2
#ch12    6         2

sum(readSum[,,"FALSE"])
#  30 
sum(readSum[,,"TRUE"])
# 151

# Look for duplicate names

TsayFiles$ch01$exercises[, 2]
TsayFileNames <- sapply(TsayFiles, function(x)
                        sapply(x, function(x2)x2[, 2]))
uniqueFileNames <- table(unlist(TsayFileNames))
duplicatedFiles <- uniqueFileNames[uniqueFileNames>1]
sort(duplicatedFiles)
# 14 duplicated, 2 triplicated
# = 18 duplicated files 

# = 151 dowloads
#  - 18 duplicates
# = 133 total files created.  

# Spot checks suggest that the
# URLs that generated files with 0 KB 
# did not seem to be accessible by a viewer.  

##
## 6.  Save 'TsayFiles' object 
##
save(TsayFiles, file="TsayFiles.rda")


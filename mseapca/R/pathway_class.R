pathway_class <-
function (filepath) {

if (file.exists(filepath)==FALSE){
stop("File not found")
}

a1 <- dirname(filepath)
setwd(a1)
untar(filepath)

# split of filepath
a <- strsplit(filepath,"/")
b <- a[[1]]

filename <- b[length(b)]

# delete extension
r <- substr(filename,1,nchar(filename)-7)

# directory for extraction
paste(a1,r,sep="/")

# set current directory
setwd(paste(a1,r,sep="/"))

# file list
files <- list.files()

# check (xml file only)
a2 <- length(files) # file list (xml)

# get current directory
cur_dir <- getwd()

A3 <- NaN;
Z1 <- NaN; Z2 <- NaN
Name <- NaN

k1 <- 1
for (i in 1:a2){
A3[i] <- substr(files[i],nchar(files[i])-3,nchar(files[i]))==".xml"
if (A3[i]==1){
w <- paste(cur_dir,files[i],sep="/")
z <- kgml_i(w)
# Z1[k1] <- list(z)

Z1[k1] <- z # metabolite ID
Name[k1] <- names(z)

k1 <- k1+1
}
}

names(Z1) <- Name
return(Z1)
}

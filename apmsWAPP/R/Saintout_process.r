#######  Saint output correction for Tab-misplacements:  function "new.saintoutput"
# Input: file-name of original saint interaction-output = unique-interaction table  -> define path before!
# Output: corrected saint-output of interactions


            ### Correction SAINT output  ###

saintout <- function(test) {
#require(seqinr)
test <- test[2:length(test)]                 
test <- strsplit(test,split="\t")
length(test)                                
#table(unlist(lapply(test, function(x){length(x)} ) ))     #check number of columns: expected 12

lost12 <- which(is.na(unlist(lapply(test, function(x){x[12]} ) ) ))  #column 12 missing

for ( i in lost12) {
a <- strsplit(test[[i]][10], split="|", fixed=TRUE)
ll <- length(a[[1]])
if (ll==1)                          # no ctrl counts, pure shifting of columns 11+12 to the previous ones
  {test[[i]][12]<-test[[i]][11] ;
   test[[i]][11]<-test[[i]][10] ;
   test[[i]][10] <- NA
  }
else                                # column 11 displaced in column 10 
  {
  if(is.integer(a[[1]][ll])==FALSE)
  {test[[i]][12]<-test[[i]][11] ;
  test[[i]][11]<- a[[1]][ll] ;
  g <- c2s(paste(a[[1]][-ll],sep="|",collaps=""))
  test[[i]][10] <- substr(g,1,nchar(g)-1) }
  else print("different case in:",i )
  }
}
return(test)
}
#######
                                                                                 
            ########## SAINT OUTPUT in Data Frame ##############

new.saintoutput <- function(filename){
#require(seqinr)
saint <- readLines(filename)                      # read original saint output
saint.out <- saintout(saint)                      # function to correct tab-mistakes from the original Saint output
saint.out2 <- as.data.frame(do.call("rbind",saint.out))
colnames(saint.out2) <- unlist(strsplit(readLines(filename)[1], split="\t", fixed=TRUE))
return(saint.out2)
}

# nochars.R
# remove letters and characters from a string

nochars<-function(text){
   alphabet<-c('a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',
               'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
               "\\(","\\)","\\.")
   n<-length(alphabet)
   compressed<-gsub(alphabet[1],"",text);
   for (i in 2:n){
      compressed<-gsub(alphabet[i],"",compressed)
   }
   return(compressed)
}

# example
nochars(text='adrian.(0)')
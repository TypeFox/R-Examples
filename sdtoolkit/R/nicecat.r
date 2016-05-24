
#Function for printing out lots of text prettily -
#specifically, will take an unbroken character string and format it
#to print out within the specifed buffer limits [ options('width') ]
#and (hopefully, check to break for spaces)

nicecat <- function(charstring, nobreak=TRUE, width=getOption("width")){

  if(nobreak){  #split words up by spaces, and also remove extra whitespace

    nolbreaks <- strsplit(charstring,"\n")
    nolbreaksvec <- paste(nolbreaks[[1]],collapse=" ")
    aswords <- strsplit(nolbreaksvec," ")
    cat(aswords[[1]][aswords[[1]]!=""],fill=TRUE)


  } else{   #naive way (doesn't care about keeping words together):
    nc <- nchar(charstring)

    nlines <- ((nc-1)%/%width)+1

    startpts <- seq(1,(nlines-1)*width+1,width)
    endpts   <- c(1:nlines)*width
    endpts[length(endpts)] <- nc #replace last line

    for (i in 1:nlines){

      cat(substr(charstring,startpts[i],endpts[i]),"\n")

    }

  }

}
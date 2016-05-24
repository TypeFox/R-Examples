`DUMPLOC`<-function(zloc, dig=12)
  {
	if(missing(dig)) dig =12

    nam = deparse(substitute(zloc))
    cat(file="",paste(sep="", nam,"=list()"), fill=TRUE)
if(length(zloc$x)>0)
{
    cat(file="",paste(sep="", nam, "$x=c(", paste(format(zloc$x, digits=dig), collapse=","), ")") , fill=TRUE)
}
if(length(zloc$y)>0)
{
    cat(file="",paste(sep="", nam,"$y=c(", paste(format(zloc$y, digits=dig), collapse=","), ")") , fill=TRUE)
}

if(length(zloc$lat)>0)
{
    cat(file="",paste(sep="", nam,"$lat=c(", paste(format(zloc$lat, digits=dig), collapse=","), ")") , fill=TRUE)
}

if(length(zloc$lon)>0)
{
    cat(file="",paste(sep="", nam,"$lon=c(", paste(format(zloc$lon, digits=dig+1), collapse=","), ")") , fill=TRUE)
}

if(length(zloc$LON)>0)
{
    cat(file="",paste(sep="", nam,"$LON=c(", paste(format(zloc$LON, digits=dig+1), collapse=","), ")") , fill=TRUE)
}


if(length(zloc$LAT)>0)
{
    cat(file="",paste(sep="", nam,"$LAT=c(", paste(format(zloc$LAT, digits=dig), collapse=","), ")") , fill=TRUE)
}

}



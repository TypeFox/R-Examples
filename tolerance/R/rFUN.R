#Internal Function

rFUN <- function(FUN, r1 = "1", r2 = "2")
{
FUN <- match.fun(FUN)
FUN.char <- deparse(body(FUN))
temp.names <- names(formals(FUN))
name.length <- which.max(nchar(temp.names))
if(name.length==1){
	FUN.char <- gsub(temp.names[1], r1, FUN.char)
	FUN.char <- gsub(temp.names[2], r2, FUN.char)
} else{
	FUN.char <- gsub(temp.names[2], r2, FUN.char)
	FUN.char <- gsub(temp.names[1], r1, FUN.char)
}
FUN.char
}


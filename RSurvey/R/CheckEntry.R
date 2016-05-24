# Controls the character string content within a Tk entry widget.

CheckEntry <- function(obj.class, ent.str="") {

  if (ent.str == "")
    return(ent.str)

  if ("numeric" %in% obj.class) {
    accept.vals <- c(as.character(0:9), "-", "e", "E", ".", "N", "A")
  } else if ("integer" %in% obj.class) {
    accept.vals <- c(as.character(0:9), "-", "e", "E", "N", "A")
  } else if ("logical" %in% obj.class) {
    accept.vals <- c("T", "R", "U", "E", "F", "A", "L", "S", "N")
  } else {
    return(ent.str)
  }

  chr <- unlist(strsplit(ent.str, split=""))
  if (all(chr %in% accept.vals))
    rtn <- ent.str
  else
    rtn <- paste(chr[chr %in% accept.vals], collapse="")

  return(rtn)
}

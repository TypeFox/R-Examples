nippon.palette <- function(){
  palette(c("#2A2A2A", # Kuro (black)
            "#BE0032", # Aka (red)
            "#00B66E", # Midori (green)
            "#006AB6", # Ao (blue)
            "#008E94", # Aomidori (cyan)
            "#DA508F", # Akamurasaki (magenta)
            "#E3C700", # Kiiro (yellow)
            "#767676"  # Haiiro (gray)
            ))
##  cat("use 'palette(\"default\")' to restore original color palette\n")
}

JapaneseColors <- function(names){
  findcolor <- function(n) {
    x <- c(jiscolors[grep(n,jiscolors$Kanji),"RGB"],
           jiscolors[grep(n,jiscolors$Kana),"RGB"],
           jiscolors[grep(n,jiscolors$Romaji),"RGB"])
    if(length(x)>1) warning(paste("two or more colors were resulted as",n))
    return(x)
  }
  res <- sapply(names, findcolor)
  return(res)
}

  

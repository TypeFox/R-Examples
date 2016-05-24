lilyinput <- function(X, file = "Rsong.ly", 
    Major = TRUE, key = "c", clef = c("treble", "bass", "alto", "tenor"), 
    time = "4/4", endbar = TRUE, midi = TRUE, tempo = "2 = 60", 
    textheight = 220, linewidth = 150, indent = 0, fontsize = 14)
{
  clef <- match.arg(clef)
  
  # notes, 97 entries in the pot (a,,, - a'''''):
  if(Major){
  pot <- 
    switch(key,
        d = c("c", "cis", "d", "dis", "e", "f", "fis", "g", "gis", "a", "bes", "b"),
        e = c("c", "cis", "d", "dis", "e", "f", "fis", "g", "gis", "a", "ais", "b"),
        f = c("c", "cis", "d", "es", "e", "f", "fis", "g", "as", "a", "bes", "b"),
        g = c("c", "cis", "d", "es", "e", "f", "fis", "g", "gis", "a", "bes", "b"),
        a = c("c", "cis", "d", "dis", "e", "f", "fis", "g", "gis", "a", "bes", "b"),
        b = c("c", "des", "d", "es", "e", "f", "fis", "g", "as", "a", "bes", "b"),
        es = c("c", "des", "d", "es", "e", "f", "ges", "g", "as", "a", "bes", "b"),
        c("c", "cis", "d", "es", "e", "f", "fis", "g", "gis", "a", "bes", "b")
    )
  }else{
  pot <- 
    switch(key,
        h = c("c", "cis", "d", "dis", "e", "f", "fis", "g", "gis", "a", "bes", "b"),
        cis = c("c", "cis", "d", "dis", "e", "f", "fis", "g", "gis", "a", "ais", "b"),
        d = c("c", "cis", "d", "es", "e", "f", "fis", "g", "as", "a", "bes", "b"),
        e = c("c", "cis", "d", "es", "e", "f", "fis", "g", "gis", "a", "bes", "b"),
        fis = c("c", "cis", "d", "dis", "e", "f", "fis", "g", "gis", "a", "bes", "b"),
        g = c("c", "des", "d", "es", "e", "f", "fis", "g", "as", "a", "bes", "b"),
        c = c("c", "des", "d", "es", "e", "f", "ges", "g", "as", "a", "bes", "b"),
        c("c", "cis", "d", "es", "e", "f", "fis", "g", "gis", "a", "bes", "b")
    )
  }  

  pot <- unlist(lapply(
    c(",,,", ",,", ",", "", "'", "''", "'''", "''''", "'''''"), 
        function(x) paste(pot, x, sep="")))[-c(1:9, 107:108)]
  # Initializing
  slur <- toene <- character(length(X$note)) 
  
  # note: pitch, length, punctuation
  note <- ifelse(is.na(X$note), "r", pot[X$note + 49])
  duration <- ifelse(X$duration %in% 2^(0:8), X$duration, "")
  punctuation <- ifelse(X$punctuation, ".", "")

  # start/end of slurs:
  if(sum(X$slur) %% 2) 
    stop("More starting than ending slurs")
  slur[which(X$slur)] <- 
    rep(c("(", ")"), sum(X$slur) %/% 2)
  # join notes:
  toene <- ifelse(slur == ")", 
    paste(note, duration, punctuation, slur, sep = ""),
    paste(note, duration, punctuation, slur, sep = ""))

  # mode
  mode <- if(Major) "\\major" else "\\minor"        

 # key
  if(Major){
    pot <- c("fis" , "h" , "e" , "a" , "d" , "g" , "c" , "f" , 
        "b" , "es" , "as" , "des" , "ges") 
    if(!(key %in% pot))
        stop("Wrong key, possible major keys are:\n", 
             paste(pot, collapse = " "), "\n")
  }else{
    pot <- c("cis" , "gis" , "dis" , "fis" , "h" , "e" , "a" , 
        "d" , "g" , "c" , "f" , "b" , "es")
    if(!(key %in% pot))
        stop("Wrong key, possible minor keys are:\n", 
             paste(pot, collapse = " "), "\n")
  }
  if(key == "b"){key <- "bes"
  }else {if(key == "h") key <- "b"}

  # generate LilyPond file:
  write(file = file,
        c(paste("#(set-global-staff-size", fontsize, ")"),
        "\\header{tagline = \"\"}",
        "\\melody = {",
        paste("    \\time", time),
        paste("    \\key", key, mode),
        paste("    \\clef", clef),
        paste("   ", toene),
    if(endbar) 
        "   \\bar \"|.\"",
        "  }",
        "  \\paper{",
        paste("    textheight = ", textheight, ".\\mm", sep = ""),
        paste("    linewidth = ", linewidth, ".\\mm", sep = ""),
        paste("    indent = ", indent, ".\\mm", sep = ""),
        "  }",  
     "   \\score{",
     "   \\melody",
     "   \\layout{ }",
     if(midi){ c("    \\midi{",  paste("      \\tempo", tempo), "  }")},
    "   }"))
}

## test:
#  X <- data.frame(note = c(3, 0, 1, 3, -4, -2, 0, 1, 3, 1, 0, -2, NA),  
#    duration = c(2, 4, 8, 2, 2, 8, 8, 8, 8, 4, 4, 1, 1), 
#    punctuation = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), 
#    slur = c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
#
#  lilyinput(X, file = "c:/test2.ly")
#  
# Y <- data.frame(note = c(3, 5, 7, 8, 10, 10, 12, 12, 12, 12, 10, 12, 12, 12, 12, 10, 8, 8, 8, 8, 7, 7, 5, 5, 5, 5, 3, NA),  
#    duration = c(4, 4, 4, 4, 2, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 2, 2, 4, ,4, 4, 4, 2, 1), 
#    punctuation = FALSE, 
#    slur = FALSE)
#    
#lilyinput(Y, file="c:/Bsp.ly")

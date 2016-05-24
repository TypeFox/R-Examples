"odfTranslate" <-
  function(x, toR = TRUE)
{
  longDash <- "\342\200\223"
  Encoding(longDash) <- "UTF-8"
  badQuote1 <- "\342\200\235"
  Encoding(badQuote1) <- "UTF-8"
  badQuote2 <- "\342\200\234"
  Encoding(badQuote2) <- "UTF-8"  
  if(toR)
    {
      x <- gsub("&gt;", ">", x)
      x <- gsub("&lt;", "<", x)
      x <- gsub("&quot;", "\"", x)
      x <- gsub("&apos;", "\'", x)
      x <- gsub("&amp;", "&", x)
      x <- gsub(longDash, "-", x)
      x <- gsub(badQuote1, "\"", x)
      x <- gsub(badQuote2, "\"", x)
    } else {
      ##TODO:  Is this code ever used? (20060808, nc)
      ## & must be the first here
      x <- gsub("&","&amp;",  x)
      x <- gsub(">", "&gt;", x)
      x <- gsub("<", "&lt;",  x)
      x <- gsub("\"", "&quot;",  x)
      x <- gsub( "\'", "&apos;", x)

      ## ODF puts in a code for areas where there are 2+ characters
      ## of white space, probably because XML standard dictates that whitespace
      ## in element content must be normalized. First, we can figure out all of the lengths of
      ## white space characters and convert this to
      ## '<text:s text:c="#"/>'  where # is the length
      spaceCount <- unique(unlist(lapply(strsplit(x, "[^ ]"), function(x) as.numeric(names(table(nchar(x)))))))
      spaceCount <- spaceCount[spaceCount > 1]

      for(i in sort(spaceCount, decreasing = TRUE))
        {
          from <- paste(rep(" ", i), collapse = "")
          to <- paste("<text:s text:c=\"", i, "\"/>", sep = "")
          x <- gsub(from, to, x)
        }
    }
  x
}


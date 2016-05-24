"as.anova" <-
  function(df, heading)
{
  if(!inherits(df, "data.frame"))
    stop("df must be a data frame")
  attr(df, "heading") <- heading
                                        #if the "class" attribute of df already starts with "anova" return(df)
  if(inherits(df, "anova")) {
    dfClasses <- attr(df, "class")
    if(dfClasses[1] == "anova")
      return(df)
  }
  class(df) <- unique(c("anova", class(df)))
  df
}

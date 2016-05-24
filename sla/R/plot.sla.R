plot.sla <- function(x, modelType2Plot = 'A', ...) {
  # For backwards compatibility
  slaObj <- x
  ## construct new data frame
  input.df <- slaObj$INPUT.df
  summary.input <- summary(slaObj$INPUT.df)
  ## 
  group <- input.df[,1] 
  x <- input.df[,2]
  y <- input.df[,3]
  new.df <- data.frame(group, x, y)
  summary.new.df <- summary(new.df)
  n1 <- with(new.df, table(group))[1]
  n2 <- with(new.df, table(group))[2]
  new.group <- c(rep('one', n1), rep('two', n2))
  new.df.2 <- data.frame(new.group, x, y)
  ##  
  with(new.df, plot(y ~ x, xlab = names(input.df)[2], ylab = names(input.df)[3], 
                    col = ifelse(new.group == 'one', 'blue', 'red'), lwd = 2))
  ## get names groups from original data frame for title
  blue.names <- as.character(with(new.df, levels(group)[1]))
  red.names <- as.character(with(new.df, levels(group)[2]))
  ## create subtitle string
  sub.title.string <- paste(blue.names, "[blue]", ".......", red.names, "[red]" )
  ##
  switch(modelType2Plot, 
         A = {
           my_title <- 'Model A: 4 Parms Estimated\nInd Intercepts & Ind Slopes'
           coef1 <- slaObj$Mod.A$coef[1:2]
           coef2 <- slaObj$Mod.A$coef[3:4]
         } ,
         B = {
           my_title <- 'Model B: 2 Parms Estimated\nCom Intercept & Com Slope'
           coef1 <- slaObj$Mod.B$coef[1:2]
           coef2 <- NULL
         } ,
         C = {
           my_title <- 'Model C: 3 Parms Estimated\nInd Intercepts & Com Slope'
           coef1 <- slaObj$Mod.C$coef[c(1,3)]
           coef2 <- slaObj$Mod.C$coef[c(2,3)]
         } ,
         D = {
           my_title <- 'Model D: 3 Parms Estimated\nCom Intercept & Ind Slopes'
           coef1 <- slaObj$Mod.D$coef[c(1,2)]
           coef2 <- slaObj$Mod.D$coef[c(1,3)]
         }
  )
  
  title(my_title, sub = sub.title.string)
  
  # TODO: add these colors as function arguments so user can define
  if(modelType2Plot == "B") {
    col1 <- "black"
  } else {
    col1 <- "blue"
  }
  abline(coef1, col = col1, lwd = 3)
  if(!is.null(coef2)) {
    abline(coef2, col = "red", lwd = 3)
  }
  invisible()
}

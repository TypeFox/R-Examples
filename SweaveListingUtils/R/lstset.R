lstset <- function(taglist, LineLength = getOption("width"),
                   startS = "\\lstset{"){
   print(taglist, LineLength = LineLength, offset.start = nchar(startS),
         withFinalLineBreak = FALSE, first.print = startS)
   cat("}%\n")
   return(invisible())
}

lstdefRstyle <- function(Rset = NULL, LineLength = getOption("width"),
                    add = TRUE){

   oN <- paste("RstyleO", .numberofRstyleDefs, sep = "")
   nN <- paste("RstyleO", .numberofRstyleDefs + 1, sep = "")

   if(add) Rset <- c("style"= oN, Rset)
   if(!is(Rset, "taglist")) Rset <- taglist(list=Rset)

   if(is.null(.numberofRstyleDefs)) .numberofRstyleDefs <<- 1
      else .numberofRstyleDefs <<- .numberofRstyleDefs + 1


   cat("\n")
   lstset(Rset, LineLength = LineLength,
          startS = paste("\\lstdefinestyle{", nN,"}{",sep=""))
   cat("\\lstdefinestyle{Rstyle}{style=",nN,"}\n%\n",sep="")

   return(invisible())
}

.preOrappend <- function(SetFromOptions, NewSet, append = TRUE,
                         withRstyle = FALSE){
#       cat("%\n%",SetFromOptions,"\n%\n")
       set0 <- getSweaveListingOption(SetFromOptions)
#       cat("%\n%",paste("%",set0,"\n"),"\n%\n")
       if(length(NewSet)){
          newnms <- names(NewSet)
          oldnms <- names(set0)
          ooldnms <- oldnms[! (oldnms %in% newnms)]
          if(append)  ret <- c(set0[ooldnms], NewSet)
          else        ret <- c(NewSet, set0[ooldnms])
       }else {
          ret <- set0
       }
       if(withRstyle) if(is.null(ret$style)) ret <- c(style="Rstyle",ret)
       return(ret)
}

.lstsetTempl <- function(set = NULL, setName = "",
                         LineLength = getOption("width"),
                         add = TRUE, startS = "", append = TRUE,
                         withRstyle = FALSE, withPre = TRUE,
                         withPre2 = TRUE, withAft = TRUE){
   if(add) set <- .preOrappend(setName, set, append = append,
                               withRstyle = withRstyle)
   if(!is(set, "taglist")) set <- taglist(list=set)
   if(withPre) cat(ifelse(withPre2,"\n",""),
                   "%----------------\n", sep = "")
   lstset(set, LineLength = LineLength, startS = startS)
   if(withAft) cat("%----------------\n")
   return(invisible())
}

lstsetR <- function(Rset = NULL, LineLength = getOption("width"),
                    add = getSweaveListingOption("addRset"),
                    startS = "\\lstset{",
                    append = TRUE, withRstyle = FALSE)
   .lstsetTempl(set = Rset, setName = "Rset", LineLength = LineLength,
                add = add, startS = startS, append = append,
                withRstyle = withRstyle)

lstsetRd <- function(Rdset = NULL, LineLength = getOption("width"),
                    add = getSweaveListingOption("addRdset"),
                    startS = "\\lstset{",append = TRUE)
   .lstsetTempl(set = Rdset, setName = "Rdset", LineLength = LineLength,
                add = add, startS = startS, append = append)


lstsetRin <- function(Rinset = NULL, LineLength = getOption("width"),
                    add = getSweaveListingOption("addRinset"),
                    startS = "\\lstdefinestyle{Rinstyle}{",
                    append = TRUE)
   .lstsetTempl(set = Rinset, setName = "Rin", LineLength = LineLength,
                add = add, startS = startS, append = append,
                withRstyle = TRUE)

lstsetRout <- function(Routset = NULL, LineLength = getOption("width"),
                    add = getSweaveListingOption("addRoutset"),
                    startS = "\\lstdefinestyle{Routstyle}{",
                    append = TRUE)
   .lstsetTempl(set = Routset, setName = "Rout", LineLength = LineLength,
                add = add, startS = startS, append = append)

lstsetRcode <- function(Rcodeset = NULL, LineLength = getOption("width"),
                    add = getSweaveListingOption("addRcodeset"),
                    startS = "\\lstdefinestyle{Rcodestyle}{",
                    append = TRUE)
   .lstsetTempl(set = Rcodeset, setName = "Rcode", LineLength = LineLength,
                add = add, startS = startS, append = append,
                withRstyle = TRUE)

lstsetRall <- function(Rallset = NULL, LineLength = getOption("width"),
                    add = c("in" = getSweaveListingOption("addRinset"),
                            "out" = getSweaveListingOption("addRoutset"),
                            "code" = getSweaveListingOption("addRcodeset")),
                    startS = c("in" = "\\lstdefinestyle{Rinstyle}{",
                               "out" = "\\lstdefinestyle{Routstyle}{",
                               "code" = "\\lstdefinestyle{Rcodestyle}{"),
                    append = c("in" = TRUE, "out" = TRUE, "code" = TRUE),
                    withOptionsDefAppend = TRUE){

   .prep <- function(u){
      u <- rep(u, length.out = 3)
      ioc <- c("in", "out", "code")
      if(all(ioc %in% names(u)))
         u <- u[ioc]
      names(u) <- ioc
      u}

   add    <- .prep(add)
   startS <- .prep(startS)
   append <- .prep(append)

   cat("%----------------\n")
   lstdefRstyle(Rset = Rallset, LineLength = LineLength)
   if(withOptionsDefAppend){
     .lstsetTempl(set = Rallset, setName = "Rin", LineLength = LineLength,
                add = add["in"], startS = startS["in"], append = append["in"],
                withRstyle = TRUE, withPre = FALSE, withAft = FALSE)
     .lstsetTempl(set = Rallset, setName = "Rout", LineLength = LineLength,
                add = add["out"], startS = startS["out"], append = append["out"],
                withRstyle = TRUE, withAft = FALSE, withPre = FALSE)
     .lstsetTempl(set = Rallset, setName = "Rcode", LineLength = LineLength,
                add = add["code"], startS = startS["code"], append = append["code"],
                withRstyle = TRUE,  withPre = FALSE)
   }else{
     cat("\\lstdefinestyle{Rinstyle}{style=Rstyle}\n")
     cat("\\lstdefinestyle{Routstyle}{style=Rstyle}\n")
     cat("\\lstdefinestyle{Rcodestyle}{style=Rstyle}\n")
   }
   return(invisible())
   }

mdf <-
function (df, pedcols = c(1:3), factorcols = NULL, ycols = NULL, sexcode = NULL, keep=F
           ,relmat = NULL) {
#  mdf()  - version 4 of mdf
#         - remove Id == NA
#         - remove duplicate Id's - including first dup
#         - add SId's which dont match Id to base 
#         - add DId's which dont match Id to base
#         - renumber all Id's
#         - make the original Id a row no
#         - if keep == T - retain  unused columns of df
#         - if keep == F - dont retain unused columns  of df
#         - always retain Id's and factors and add matrix of multivariate traits
#         - Sex should be one of the factors
#         - transform sex codes to NA if not in sexcode[]
#         - the first entry in sexcode[] should be the heterogametic sex
#     Note:  pedcols and factorcols must be either both numeric ( ie col nos) 
#                                              or both character ( ie col names)
#            ycols can be either numeric or character
#            sexcode should be of the same type as df$Sex
#
# check sexcode present
   if(is.null(sexcode)) {
     stop("mdf(): sexcode argument must be present.\n")
   }
#
# check sexcode type
   if(is.character(df$Sex) & is.character(sexcode)){
     sexcharacter <- T
   }
   else if(is.numeric(df$Sex) &  is.numeric(sexcode)){
     sexcharacter <- F
   }
   else if(is.factor(df$Sex)){
     if(is.character(sexcode)){
       sexcharacter <- T
     }
     else if(is.numeric(sexcode)){
       sexcharacter <- F
     }
     else {
       stop("mdf(): df$Sex is a factor and sexcode is not either numeric or character.\n")
     }
   }
   else {
     stop("mdf(): sexcode and df$Sex not the same type.\n")
   }
#
   cat("Pedigree Id check:\n")
   cat("No of rows with Id in original dataframe = ",length(df$Id),"\n")
#
#  transform Sex code to NA if not in sexcode[]
#  df$Sex[match(df$Sex,sexcode,nomatch=0)] <- match(df$Sex,sexcode)
   sexmatch <- match(df$Sex,sexcode)
   count <- 0
   for(i in 1:length(df$Sex)) {
     if(is.na(sexmatch[i])) {
       df$Sex[i] <- NA
       count <- count + 1
     }
   }
   cat("No of sex codes not in sexcode[] so changed to NA = ",count,"\n")
#
#  Remove NA's for Sex
   df1 <- data.frame(df[!is.na(df$Sex),])
   nona <- length(df$Sex) - length(df1$Sex)
   cat("No of rows with Sex == NA removed from dataframe = ",nona,"\n")
#
#  Remove NA's for Id
    df2 <- data.frame(df1[!is.na(df1$Id),])
    nona <- length(df1$Id) - length(df2$Id)
    cat("No of rows with Id == NA removed from dataframe = ",nona,"\n")
#
# Remove duplicate Id's
    df3 <- df2[!duplicated(df2$Id) & !duplicated.first(df2$Id), ]
    nodup <- length(df2$Id) - length(df3$Id)
    cat("No of rows with duplicated Id removed from dataframe = ",nodup,"\n")
    cat("No of rows remaining after duplicates and NA's removed = ",nrow(df3),"\n")
#
# Make a df of unmatched SId's
    SId.r <- match(df3$SId, df3$Id)  # unmatched are NA
    SId.um <- is.na(SId.r)  # T if NA,  F if OK
    SId.um.id <- df3$SId[SId.um]  # list of all unmatched SId's
    SId.um.id.unique <- unique(SId.um.id)   # duplicates  removed
    SId.um.id.unique <- SId.um.id.unique[!is.na(SId.um.id.unique)]  # NA's removed
    noum <- length(SId.um.id.unique)
    if(noum == 0) {
      noum.s <- 0
    }
    else {
      df.SId.um <- data.frame(df3[1:noum,])
      df.SId.um$Id <- SId.um.id.unique
      df.SId.um[,2:ncol(df3)] <- rep(NA,noum)
      df.SId.um[,"Sex"] <- rep(sexcode[1],noum)  # sex male coded 1
      cat("No of SId's with no matching Id = ", nrow(df.SId.um),"\n")
      noum.s <- noum
    }
#
# Make a df of unmatched DId's
    DId.r <- match(df3$DId, df3$Id)  # unmatched are NA
    DId.um <- is.na(DId.r)  # T if NA,  F if OK
    DId.um.id <- df3$DId[DId.um]  # list of all unmatched DId's
    DId.um.id.unique <- unique(DId.um.id)   # duplicates removed
    DId.um.id.unique <- DId.um.id.unique[!is.na(DId.um.id.unique)]  # NA's removed
    noum <- length(DId.um.id.unique)
    if(noum == 0) {
      noum.d <- 0
    }
    else {
      df.DId.um <- data.frame(df3[1:noum,])
      df.DId.um$Id <- DId.um.id.unique
      df.DId.um[,2:ncol(df3)] <- rep(NA,noum)
      df.DId.um[,"Sex"] <- rep(sexcode[2],noum)  # sex female coded 2
      cat("No of DId's with no matching Id = ", nrow(df.DId.um),"\n")
      noum.d <- noum
    }
#
# Merge the unmatched Sid's and DId's with df3
    if(noum.s > 0 && noum.d > 0) {
      rnames <- c(as.character(df.SId.um$Id),as.character(df.DId.um$Id),as.character(df3$Id))
      df4 <- data.frame(rbind(df.SId.um, df.DId.um, df3),check.rows=T,row.names=rnames)
    }
    else if(noum.s > 0 && noum.d == 0) {
      rnames <- c(as.character(df.SId.um$Id),as.character(df3$Id))
      df4 <- data.frame(rbind(df.SId.um,df3),check.rows=T,row.names=rnames)
    }
    else if(noum.s == 0 && noum.d > 0){
      rnames <- c(as.character(df.DId.um$Id),as.character(df3$Id))
      df4 <- data.frame(rbind(df.DId.um, df3),check.rows=T,row.names=rnames)
    }
    else if(noum.s == 0 && noum.d == 0) {
      df4 <- df3
    }
    cat("Length of dataframe with base Id's added = ",nrow(df4),"\n")
#
# Make sure all factorcols are factors in df4
    for(i in factorcols){
      df4[,i] <- factor(df4[,i])
    }
# Do keep option
    if(!keep) {
      # dont keep unused cols
#      df5 <- data.frame(df4[,pedcols], df4[,factorcols])
      df5 <- df4[,c(pedcols,factorcols)]
    }
    else if (keep){
      # keep all cols including ycols
      df5 <- df4
    }
#
# Renumber Id's
    cat("Renumber pedigree Id's:\n")
    mdf <- pedrenum(df5)
#
# Make multivariate df
    cat("Add matrix of multivariate traits:\n")
    mdf[, "Ymat"] <- as.matrix(df4[, ycols]) 
    # note gets ycols from df4 (not in df5 if keep == F)
#
# Do relmat option
    if(is.null(relmat)) {
      cat("Return mdf as a normal dataframe:\n")
      return(mdf)
    }
    else {
      # setup for nadiv()
      cat("Setup pedigree for nadiv():\n")
      ped <- matrix(0,length(mdf$Id),3)
      ped2 <- matrix(0,length(mdf$Id),4)
      ped[,1] <- mdf[,"Id"]
      ped[,2] <- mdf[,"DId"]
      ped[,3] <- mdf[,"SId"]
      ped2[,1:3] <- ped
      if(sexcharacter) {
        ped2[,4] <- as.character(mdf[,"Sex"])
      }
      else {
        ped2[,4] <- as.numeric(as.character(mdf[,"Sex"]))
      }
      dimnames(ped2)[[2]] <- c("ID","Dam","Sire","Sex")
      dimnames(ped)[[2]] <- c("ID","Dam","Sire")
      cat("Make relationship matrices:\n")
      a <- NULL
      aa <- NULL
      d <- NULL
      ad <- NULL
      dd <- NULL
      s <- NULL
      dsim <- NULL
      e <- NULL
      if(any(relmat == "E")) {
        e <- Diagonal(length(mdf$Id))
      }
      if(any(relmat == "A")){
        a <- makeA(ped)
      }
      if(any(relmat == "AA")){
        aa <- makeAA(ped)$AA
      }
      if(any(relmat == "D")) {
        d <- makeD(ped,invertD=F)$D
      }
      if(any(relmat == "Dsim")) {
        dsim <- makeDsim(ped,N=1000,invertD=F)$D
      }
      if(any(relmat == "AD")){
        ad <- makeDomEpi(ped,output="AD",invertD=F)$AD
      }
      if(any(relmat == "DD")){
        dd <- makeDomEpi(ped,output="DD",invertD=F)$DD
      }
      if(any(relmat == "S")){
        s <- makeS(ped2,heterogametic=sexcode[1],returnS=T,DosageComp="ngdc")$S
      }
      if(any(relmat == "S.hori")){
        s <- makeS(ped2,heterogametic=sexcode[1],returnS=T,DosageComp="hori")$S
      }
      if(any(relmat == "S.hedo")){
        s <- makeS(ped2,heterogametic=sexcode[1],returnS=T,DosageComp="hedo")$S
      }
      if(any(relmat == "S.hoha")){
        s <- makeS(ped2,heterogametic=sexcode[1],returnS=T,DosageComp="hoha")$S
      }
      if(any(relmat == "S.hopi")){
        s <- makeS(ped2,heterogametic=sexcode[1],returnS=T,DosageComp="hopi")$S
      }
      cat("Return mdf as an object of class mdf:\n")
      cat(" containing the dataframe as mdf$df:\n")
      cat(" and the relationship matrices as mdf$rel:\n")
      outlist <- list(df=mdf,rel=list(e=e,a=a,aa=aa,d=d,dsim=dsim,ad=ad,dd=dd,s=s))
      class(outlist) <- "mdf"
      return(outlist)
    }
}

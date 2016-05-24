
fix.names.manually <- function (master.list) {
  #  master.list=roster.collection[[kk]][,c("number","pos","last","first","numfirstlast")]
    
  #one name, two players.
  master.list$first[which(master.list$last=="PICARD" & master.list$first=="ALEXANDRE" & master.list$pos == "D")] <- "ALEXANDRE R."
  master.list$first[which(master.list$last=="GREEN" & master.list$first=="MIKE" & master.list$pos == "C")] <- "MICHAEL G."
  
  #manual fixes.
  master.list$last[which(master.list$last=="ANDERSSON" & master.list$first=="CRAIG")] <- "ANDERSON"
  master.list$first[which(master.list$last=="ANTROPOV" & master.list$first=="NIKOLAI")] <- "NIK"
  master.list$first[which(master.list$last=="AULD" & master.list$first=="ALEXANDER")] <- "ALEX"
  master.list$first[which(master.list$last=="AXELSSON" & master.list$first=="PER JOHAN")] <- "P.J."
  master.list$first[grep("P\\. *J\\. *",master.list$first)] <- "P.J."
  master.list$first[which(master.list$last=="BAILEY" & master.list$first=="JOSHUA")] <- "JOSH"
  master.list$first[which(master.list$last=="BARCH" & master.list$first=="KRYSTOFER")] <- "KRYS"
  master.list$first[which(master.list$last=="BARKER" & master.list$first=="CAMERON")] <- "CAM"
  master.list$first[which(master.list$last=="BERGFORS" & master.list$first=="NICKLAS")] <- "NICLAS"
  master.list$first[which(master.list$last=="BLACKBURN" & master.list$first=="DANIEL")] <- "DAN"
  master.list$first[which(master.list$last=="BLAKE" & master.list$first=="ROBERT")] <- "ROB"
  master.list$first[which(master.list$last=="BLUNDEN" & master.list$first=="MICHAEL")] <- "MIKE"
  master.list$first[which(master.list$last=="BOURQUE" & master.list$first=="CHRISTOPHER")] <- "CHRIS"
  master.list$first[which(master.list$last=="BOYNTON" & master.list$first=="NICHOLAS")] <- "NICK"
  master.list$first[which(master.list$last=="BRIERE" & master.list$first=="DANNY")] <- "DANIEL"
  master.list$first[which(master.list$last=="BRYZGALOV" & master.list$first=="ILJA")] <- "ILYA"
  master.list$first[which(master.list$last=="BURROWS" & master.list$first=="ALEXANDRE")] <- "ALEX"
  master.list$first[which(master.list$last=="CAMMALLERI" & master.list$first=="MICHAEL")] <- "MIKE"
  master.list$first[which(master.list$last=="CARCILLO" & master.list$first=="DANIEL")] <- "DAN"
  master.list$first[which(master.list$last=="CARLE" & master.list$first=="MATTHEW")] <- "MATT"
  master.list$first[which(master.list$last=="CLEARY" & master.list$first=="DAN")] <- "DANIEL"
  master.list$first[which(master.list$last=="CLEARY" & master.list$first=="DANNY")] <- "DANIEL"
  master.list$first[which(master.list$last=="CORVO" & master.list$first=="JOSEPH")] <- "JOE"
  master.list$first[which(master.list$last=="CRABB" & master.list$first=="JOSEPH")] <- "JOEY"
  master.list$first[which(master.list$last=="CROMBEEN" & master.list$first=="BJ")] <- "B.J."
  master.list$first[which(master.list$last=="CROMBEEN" & master.list$first=="BRANDON")] <- "B.J."
  master.list$first[which(master.list$last=="DADONOV" & master.list$first=="EVGENII")] <- "EVGENY"
  
  master.list$first[which(master.list$last=="DOWD" & master.list$first=="JAMES")] <- "JIM"
  master.list$first[which(master.list$last=="DOWELL" & master.list$first=="JACOB")] <- "JAKE"
  master.list$first[which(master.list$last=="DRAZENOVIC" & master.list$first=="NICHOLAS")] <- "NICK"

  master.list$first[which(master.list$last=="DUMBA" & master.list$first=="MATHEW")] <- "MATT"
  
  master.list$first[which(master.list$last=="DUMONT" & master.list$first=="J P")] <- "JEAN-PIERRE"
  master.list$first[which(master.list$last=="DUMONT" & master.list$first=="J-P")] <- "JEAN-PIERRE"
  master.list$first[which(master.list$last=="EARL" & master.list$first=="ROBBIE")] <- "ROBERT"

  master.list$first[which(master.list$last=="ENSTROM" & master.list$first=="TOBY")] <- "TOBIAS"  ## 2014-10-19
  
  master.list$first[which(master.list$last=="FERNANDEZ" & master.list$first=="EMMANUEL")] <- "MANNY"
  master.list$first[which(master.list$last=="FROLOV" & master.list$first=="ALEXANDER")] <- "ALEX"
  master.list$first[which(master.list$first=="TJ")] <- "T.J."
  master.list$first[which(master.list$last=="GAUTHIER" & master.list$first=="DENIS JR.")] <- "DENIS"
  master.list$first[which(master.list$last=="GIGUERE" & master.list$first=="J")] <- "JEAN-SEBASTIEN"
  master.list$first[which(master.list$last=="GIRARDI" & master.list$first=="DAN")] <- "DANIEL"
  master.list$first[which(master.list$last=="GREENE" & master.list$first=="ANDY")] <- "ANDREW"
  master.list$first[which(master.list$last=="GREER" & master.list$first=="MICHAEL")] <- "MIKE"
  master.list$first[which(master.list$last=="GROSSMAN" & master.list$first=="NIKLAS")] <- "NICKLAS"
  master.list$last[which(master.list$last=="GROSSMAN" & master.list$first=="NICKLAS")] <- "GROSSMANN"
  master.list$first[which(master.list$last=="GUENIN" & master.list$first=="NATE")] <- "NATHAN"
  master.list$first[which(master.list$last=="HALKO" & master.list$first=="STEVE")] <- "STEVEN"
  master.list$first[which(master.list$last=="HIGGINS" & master.list$first=="CHRISTOPHER")] <- "CHRIS"

  master.list$first[which(master.list$last=="HAVLAT" & master.list$first=="MARTY")] <- "MARTIN"

  master.list$first[which(master.list$last=="HERR" & master.list$first=="MATTHEW")] <- "MATT"

  
  master.list$last[which(master.list$last=="HILLEN III" & master.list$first=="JOHN")] <- "HILLEN"
  master.list$first[which(master.list$last=="HILLEN" & master.list$first=="JOHN")] <- "JACK"
  master.list$first[which(master.list$last=="HOLIK" & master.list$first=="ROBERT")] <- "BOBBY"
  master.list$first[which(master.list$last=="HOWARD" & master.list$first=="JAMES")] <- "JIMMY"
  
  master.list$first[which(master.list$last=="IRWIN" & master.list$first=="MATTHEW")] <- "MATT"
  master.list$first[which(master.list$last=="JACKMAN" & master.list$first=="RICHARD")] <- "RIC"
  master.list$first[which(master.list$last=="JACQUES" & master.list$first=="J-F")] <- "JEAN-FRANCOIS"
  master.list$first[which(master.list$last=="JOHANSSON" & master.list$first=="MATTIAS")] <- "MATHIAS"
  master.list$first[which(master.list$last=="KALINSKI" & master.list$first=="JONATHON")] <- "JON"
  
  master.list$last[which(master.list$last=="KASTSITSYN")] <- "KOSTITSYN"
  master.list$first[which(master.list$last=="KOSTITSYN" & master.list$first=="SIARHEI")] <- "SERGEI"

  master.list$first[which(master.list$last=="KILLORN" & master.list$first=="ALEXANDER")] <- "ALEX"
  master.list$first[which(master.list$last=="KING" & master.list$first=="DWAYNE")] <- "D.J."
  
  master.list$first[which(master.list$first=="DJ")] <- "D.J."
  master.list$first[which(master.list$last=="KNUBLE" & master.list$first=="MICHAEL")] <- "MIKE"
  master.list$first[which(master.list$last=="KOLANOS" & master.list$first=="KRYSTOFER")] <- "KRYS"
  master.list$first[which(master.list$last=="KOMISAREK" & master.list$first=="MICHAEL")] <- "MIKE"
  master.list$first[which(master.list$last=="KONDRATIEV" & master.list$first=="MAX")] <- "MAXIM"
  master.list$first[which(master.list$last=="KOVALEV" & master.list$first=="ALEXEI")] <- "ALEX"
  master.list$last[which(master.list$last=="KRONVALL")] <- "KRONWALL"
  master.list$first[which(master.list$last=="LEGACE" & master.list$first=="EMMANUEL")] <- "MANNY"
  master.list$first[which(master.list$last=="LETANG" & master.list$first=="KRISTOPHER")] <- "KRIS"

  master.list$first[which(master.list$last=="MACIAS" & master.list$first=="RAYMOND")] <- "RAY"

  master.list$first[which(master.list$last=="MACLEAN" & master.list$first=="DONALD")] <- "DON"
  master.list$last[which(master.list$last=="MAGNAN-GRENIER")] <- "MAGNAN"
  master.list$first[which(master.list$last=="MAYOROV" & master.list$first=="MAXIM")] <- "MAKSIM"
  master.list$first[which(master.list$last=="MCCOLLUM" & master.list$first=="TOM")] <- "THOMAS"
  master.list$first[which(master.list$last=="MCGILLIS" & master.list$first=="DAN")] <- "DANIEL"
  master.list$last[which(master.list$last=="MEYER IV")] <- "MEYER"
  master.list$first[which(master.list$last=="MEYER" & master.list$first=="FREDDY")] <- "FREDERICK"
  master.list$first[which(master.list$last=="MILLER" & master.list$first=="ANDREW")] <- "DREW"

  master.list$first[which(master.list$last=="MILLS" & master.list$first=="BRADLEY")] <- "BRAD"
  
  master.list$first[which(master.list$last=="MODANO" & master.list$first=="MICHAEL")] <- "MIKE"
  master.list$first[which(master.list$last=="MODIN" & master.list$first=="FREDDY")] <- "FREDRIK"
  master.list$first[which(master.list$last=="NEIL" & master.list$first=="CHRISTOPHER")] <- "CHRIS"

  master.list$first[which(master.list$last=="NIETO" & master.list$first=="MATTHEW")] <- "MATT"

  master.list$first[which(master.list$last=="ODUYA" & master.list$first=="DAVID JOHNNY")] <- "JOHNNY"
  master.list$first[which(master.list$last=="ODUYA" & master.list$first=="JOHN")] <- "JOHNNY"
  master.list$last[which(master.list$last=="ORTMYER" & master.list$first=="JED")] <- "ORTMEYER"
  master.list$first[which(master.list$last=="OVECHKIN" & master.list$first=="ALEXANDER")] <- "ALEX"
  
  master.list$first[which(master.list$last=="PARENTEAU" & master.list$first=="PIERRE")] <- "P.A."
  master.list$first[which(master.list$last=="PARENTEAU" & master.list$first=="PA")] <- "P.A."
  master.list$first[which(master.list$last=="PELLEY" & master.list$first=="RODNEY")] <- "ROD"
  master.list$first[which(master.list$last=="PEVERLEY" & master.list$first=="JOHN")] <- "RICH"

  master.list$first[which(master.list$last=="POULIOT" & master.list$first=="MARC")] <- "MARC-ANTOINE"

  master.list$first[which(master.list$last=="PROSPAL" & master.list$first=="VINNY")] <- "VACLAV"
  master.list$first[which(master.list$last=="PURCELL" & master.list$first=="EDWARD")] <- "TEDDY"

  master.list$last[which(master.list$last=="PUSHKAREV" & master.list$first=="KONSTANTIN")] <- "PUSHKARYOV"

  master.list$first[which(master.list$last=="REINHART" & master.list$first=="MAXWELL")] <- "MAX"
  
  
  master.list$first[which(master.list$last=="REINPRECHT" & master.list$first=="STEVE")] <- "STEVEN"
  master.list$first[which(master.list$last=="RISSMILLER" & master.list$first=="PAT")] <- "PATRICK"
  master.list$first[which(master.list$last=="RUPP" & master.list$first=="MICHAEL")] <- "MIKE"
  master.list$first[which(master.list$last=="SANTORELLI" & master.list$first=="MICHAEL")] <- "MIKE"
  master.list$first[which(master.list$last=="SCUDERI" & master.list$first=="ROBERT")] <- "ROB"

  master.list$first[which(master.list$last=="SESTITO" & master.list$first=="TOMMY")] <- "TOM"
  master.list$last[which(master.list$last=="SHISKANOV" & master.list$first=="TIMOFEI")] <- "SHISHKANOV"
  master.list$first[which(master.list$last=="SILLINGER" & master.list$first=="MICHAEL")] <- "MIKE"
  
  master.list$first[which(master.list$last=="SIM" & master.list$first=="JON")] <- "JONATHAN"
  master.list$first[which(master.list$last=="SIMON" & master.list$first=="BEN")] <- "BENJAMIN"
  master.list$first[which(master.list$last=="STAJAN" & master.list$first=="MATTHEW")] <- "MATT"
  
  master.list$first[which(master.list$last=="STEEN" & master.list$first=="ALEXANDER")] <- "ALEX"
  master.list$last[which(master.list$last=="ST LOUIS" & master.list$first=="MARTIN")] <- "ST. LOUIS"
  master.list$first[which(master.list$last=="STORTINI" & master.list$first=="ZACHERY")] <- "ZACK"
  master.list$last[which(master.list$last=="ST PIERRE" & master.list$first=="MARTIN")] <- "ST. PIERRE"
  master.list$last[which(master.list$last=="STREBAK" & master.list$first=="MARTIN")] <- "STRBAK"
  master.list$first[which(master.list$first=="PK")] <- "P.K."

  master.list$first[which(master.list$last=="TAYLOR" & master.list$first=="TIMOTHY")] <- "TIM"
  master.list$first[which(master.list$last=="THOMAS" & master.list$first=="TIMOTHY JR.")] <- "TIM"
  master.list$first[which(master.list$last=="THOMAS" & master.list$first=="WILLIAM")] <- "BILL"


  
  master.list$first[which(master.list$first=="RJ")] <- "R.J."
  master.list$first[which(master.list$last=="VALICEVIC" & master.list$first=="ROBERT")] <- "ROB"
  master.list$first[which(master.list$last=="VALIQUETTE" & master.list$first=="STEVE")] <- "STEPHEN"
  master.list$first[which(master.list$last=="VANDERMEER" & master.list$first=="JAMES")] <- "JIM"
  master.list$first[which(master.list$last=="VARLAMOV" & master.list$first=="SIMEON")] <- "SEMYON"
  master.list$last[which(master.list$last=="VANDE VELDE" & master.list$first=="CHRIS")] <- "VANDEVELDE"

  master.list$first[which(master.list$last=="WHITE" & master.list$first=="COLIN (JOHN)")] <- "COLIN"
  
  
  master.list$first[which(master.list$last=="WOZNIEWSKI" & master.list$first=="ANDREW")] <- "ANDY"
  master.list$first[which(master.list$last=="WYMAN" & master.list$first=="JT")] <- "JAMES"

  master.list$first[which(master.list$last=="YORK" & master.list$first=="MICHAEL")] <- "MIKE"
  master.list$first[which(master.list$last=="ZHERDEV" & master.list$first=="NIKOLAY")] <- "NIKOLAI"
  master.list$first[which(master.list$last=="ZOLNIERCZYK" & master.list$first=="HARRISON")] <- "HARRY"
  
  master.list <- master.list[order(master.list$last,
                                   master.list$first),]
  #master.list$numfirstlast <- with(master.list, paste(num, first, last))
  
  
  return(master.list)
}





.simpleCap <- function(x, ch=" ", ch2=ch) {
    s <- strsplit(x, ch)[[1]]
    paste(toupper(substring(s, 1, 1)), tolower(substring(s, 2)),
          sep = "", collapse = ch)
}

manual.patches <- function (roster.unique) {

    roster.unique$first <- sapply(roster.unique$first, .simpleCap)
    roster.unique$first <- sapply(roster.unique$first, .simpleCap, "-")
    #roster.unique$first <- sapply(roster.unique$first, .simpleCap, "\\.")
    
    roster.unique$last <- sapply(roster.unique$last, .simpleCap)
    roster.unique$last <- sapply(roster.unique$last, .simpleCap, "-")
    #roster.unique$last <- sapply(roster.unique$last, .simpleCap, "\\.")
    
    substr(roster.unique$last[grep("^Mac", roster.unique$last)], 4,4) <-
        toupper(substr(roster.unique$last[grep("^Mac", roster.unique$last)], 4,4))
    substr(roster.unique$last[grep("^Mc", roster.unique$last)], 3,3) <-
        toupper(substr(roster.unique$last[grep("^Mc", roster.unique$last)], 3,3))
    roster.unique$last <- gsub("^Van ", "van ", roster.unique$last)
    roster.unique$firstlast <- paste(roster.unique$first, roster.unique$last)
    return(roster.unique)
    
}


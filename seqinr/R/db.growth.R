get.db.growth <- function(where = "ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/relnotes.txt" )
{
  if (!capabilities("http/ftp")) 
    stop("capabilities(\"http/ftp\") is not TRUE")
  ftp.proxy.bck <- Sys.getenv("ftp_proxy")
  if (ftp.proxy.bck != "") {
      warning("I'am trying to neutralize proxies")
      Sys.setenv("no_proxy" = "")
  }

  embl <- where
  tmp <- readLines( embl )
  idx <- grep("Release(.+) Month", tmp)
  tmp <- tmp[ (idx + 2):length(tmp) ]
  tmp <- strsplit( tmp, split = " " )
  not.empty <- function(x)
  {
    x <- x[nchar(x) > 0 ]
  }
  tmp <- sapply( tmp, not.empty )
  tmp <- data.frame( matrix(unlist(tmp), ncol = 4, byrow = TRUE) )
  names(tmp) <- c("Release", "Month", "Entries", "Nucleotides")

  tmp[,1] <- as.double( as.character(tmp[,1]))
  tmp[,3] <- as.double( as.character(tmp[,3]))
  tmp[,4] <- as.double( as.character(tmp[,4]))

  date  <- strsplit(as.character(tmp[,2]), split="/")
  date.to.num <- function(x)
  {
    x <- as.double( x )
    return( (x[1]-1)/12 + x[2] )
  }
  date <- sapply(date, date.to.num)
  tmp <- data.frame( cbind(tmp, date) )
  return(tmp)
}

dia.db.growth <- function( get.db.growth.out = get.db.growth(), 
  Moore = TRUE, ... )
{
  embl <- "ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/relnotes.txt"
  op <- par(no.readonly = TRUE)
  par( bg = "blue" )
  par( fg = "yellow" )
  par( col = "yellow" )
  par( col.axis = "yellow" )
  par( col.lab = "yellow" )
  par( col.main = "yellow" )
  par( col.sub = "yellow" )

  Nucleotides <- get.db.growth.out$Nucleotides
  Month <- get.db.growth.out$Month
  date <- get.db.growth.out$date
    
  plot( date, log10(Nucleotides) , pch = 20,
    main = paste("The exponential growth of the DDBJ/EMBL/Genbank content\n",
           "Last update:", 
            Month[nrow(get.db.growth.out)]),
    xlab = "Year", ylab = "Log10 number of nucleotides",
    sub = paste("Source:", embl),
    ... )
  abline(lm(log10(Nucleotides)~date),col="yellow")
  lm1 <- lm(log(Nucleotides)~date)
  mu <- lm1$coef[2] # slope
  dbt <- log(2)/mu # doubling time
  dbt <- 12*dbt # in months

  if( Moore )
  {
    x <- mean(date)
    y <- mean(log10(Nucleotides))
    a <- log10(2)/1.5
    b <- y - a*x

    for( i in seq(-10,10,by=0.5) )
      if( i != 0 )
        abline( coef=c(b+i, a), col="black" )
    legend( x = 1990, y = 7, legend= c(paste("Observed doubling time:", 
      round(dbt,1),"months"),"Moore's doubling time : 18 months"), 
      lty = c(1,1), col = c("yellow","black"))
  }
  else
  {
    legend( x = 1990, y = 7, legend=paste("Observed doubling time:", 
      round(dbt,1), "months"), lty = 1, col = "yellow")
  }
  par( op )
} 

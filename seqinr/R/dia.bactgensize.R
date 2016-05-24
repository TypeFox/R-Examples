dia.bactgensize <- function(
  fit = 2, p = 0.5, m1 = 2000, sd1 = 600, m2 = 4500, sd2 = 1000, p3 = 0.05,
  m3 = 9000, sd3 = 1000, maxgensize = 20000,
  source = c("ftp://pbil.univ-lyon1.fr/pub/seqinr/data/goldtable15Dec07.txt",
    "http://www.genomesonline.org/DBs/goldtable.txt"))
{
#
# Use local source by default:
#
  source <- source[1]
#
# Build source of data string:
#
  if(source == "ftp://pbil.univ-lyon1.fr/pub/seqinr/data/goldtable15Dec07.txt"){
    sodtxt <- "Source of data: GOLD (Genomes OnLine Database) 15 Dec 2007"
  } else {
    sodtxt <- paste("Source of data: GOLD (Genomes OnLine Database)",date())
  }
#
# Read data from GOLD:
#
  alldata <- read.table(source, header = TRUE, sep = "\t",
                        comment.char = "", quote = "")
  SUPERKINGDOM <- 1 # col number
  kingdom <- alldata[, SUPERKINGDOM]
  prodata <- alldata[ kingdom == "Archaea" | kingdom == "Bacteria", ]
  data <- prodata[, c("GENUS", "SPECIES", "SIZE.kb.")]
  names(data) <- c("genus", "species", "gs")
  data <- data[complete.cases(data), ]
#
# Remove data > maxgensize:
#
  data <- data[data$gs <= maxgensize, ]
#
# Use Kb scale
#
  sizeKb <- data$gs
  n <- length(sizeKb)
#
# Graphics
#
  x <- seq( min(sizeKb), max(sizeKb), le=200)
  mybreaks <- seq(min(sizeKb),max(sizeKb),length=15)
  vscale <- diff(mybreaks)[1]*n

  if(fit == 0)
  {
    hst <- hist(sizeKb, freq = TRUE, 
      breaks = mybreaks,
      main=paste("Genome size distribution for", n, "bacterial genomes"),
      xlab="Genome size [Kb]",
      ylab="Genome count",
      col="lightgrey")

    dst <- density(sizeKb)
    lines(x=dst$x, y=vscale*dst$y)
    legend(x = max(sizeKb)/2, y = 0.75*max(hst$counts), lty=1,
      "Gaussian kernel density estimation")
      
    mtext(sodtxt)
  }
##########################################
#
# Fitting a mixture of two normal distributions
#
###########################################
  if(fit == 2)
  {
    logvraineg <- function(param, obs) 
    {
      p <- param[1]
      m1 <- param[2]
      sd1 <- param[3]
      m2 <- param[4]
      sd2 <- param[5]

      -sum(log(p*dnorm(obs,m1,sd1)+(1-p)*dnorm(obs,m2,sd2)))
    }

    nlmres <- suppressWarnings(nlm(logvraineg, c(p, m1, sd1, m2, sd2), obs=sizeKb))
    estimate <- nlmres$estimate

    y1 <- vscale*estimate[1]*dnorm(x, estimate[2], estimate[3])
    y2 <- vscale*(1-estimate[1])*dnorm(x, estimate[4], estimate[5])
    dst <- density(sizeKb)

    hst <- hist(sizeKb, plot = FALSE, breaks = mybreaks)    
    ymax <- max(y1, y2, hst$counts, vscale*dst$y)
    
    
    hist(sizeKb, freq = TRUE, ylim=c(0,ymax),
      breaks = mybreaks,
      main=paste("Genome size distribution for", n, "bacterial genomes"),
      xlab="Genome size [Kb]",
      ylab="Genome count",
      col="lightgrey")
      
    lines(x, y1, col="red", lwd=2)
    lines(x, y2, col="blue", lwd=2)

    text(x = max(sizeKb)/2, y = ymax, pos=4, "Maximum likelihood estimates:")

    text(x = max(sizeKb)/2, y = 0.95*ymax, col="red", pos = 4, cex=1.2,
    substitute(hat(p)[1] == e1~~hat(mu)[1] == e2~~hat(sigma)[1] == e3, 
    list(e1 = round(estimate[1],3), 
       e2 = round(estimate[2],1),
       e3 = round(estimate[3],1))))

    text(x = max(sizeKb)/2, y = 0.90*ymax, col="blue", pos = 4, cex=1.2,
    substitute(hat(p)[2] == q~~hat(mu)[2] == e4~~hat(sigma)[2] == e5, 
    list(q = round(1 - estimate[1],3), 
       e4 = round(estimate[4],1),
       e5 = round(estimate[5],1))))
       
    lines(x=dst$x, y=vscale*dst$y)
    legend(x = max(sizeKb)/2, y = 0.75*ymax, lty=1,
       "Gaussian kernel density estimation")
       
    mtext(sodtxt)
  }
##########################################
#
# Fitting a mixture of three normal distributions
#
###########################################
  if(fit == 3)
  {
    logvraineg <- function(param, obs) 
    {
      p <- param[1]
      m1 <- param[2]
      sd1 <- param[3]
      m2 <- param[4]
      sd2 <- param[5]
      p3 <- param[6]
      m3 <- param[7]
      sd3 <- param[8]

      -sum(log(p*dnorm(obs,m1,sd1)
             +(1-p-p3)*dnorm(obs,m2,sd2)
             +p3*dnorm(obs,m3,sd3)))
    }

    nlmres <- suppressWarnings(nlm(logvraineg, c(p, m1, sd1, m2, sd2, p3, m3, sd3), obs=sizeKb))
    estimate <- nlmres$estimate

    y1 <- vscale*estimate[1]*dnorm(x, estimate[2], estimate[3])
    y2 <- vscale*(1-estimate[1]-estimate[6])*dnorm(x, estimate[4], estimate[5])
    y3 <- vscale*estimate[6]*dnorm(x, estimate[7], estimate[8])

    hst <- hist(sizeKb, plot = FALSE, breaks = mybreaks)    
    ymax <- max(y1, y2, y3, hst$counts)

    hist(sizeKb, freq = TRUE, ylim=c(0,ymax),
      breaks = mybreaks,
      main=paste("Genome size distribution for", n, "bacterial genomes"),
      xlab="Genome size [Kb]",
      ylab="Genome count",
      col="lightgrey")
          
    lines(x, y1, col="red", lwd=2)
    lines(x, y2, col="blue", lwd=2)
    lines(x, y3, col="green3", lwd=2)
  
    text(x = max(sizeKb)/2, y = ymax, pos=4, "Maximum likelihood estimates:")

    text(x = max(sizeKb)/2, y = 0.95*ymax, col="red", pos = 4, cex=1.2,
    substitute(hat(p)[1] == e1~~hat(mu)[1] == e2~~hat(sigma)[1] == e3, 
    list(e1 = round(estimate[1],3), 
       e2 = round(estimate[2],1),
       e3 = round(estimate[3],1))))

    text(x = max(sizeKb)/2, y = 0.90*ymax, col="blue", pos = 4, cex=1.2,
    substitute(hat(p)[2] == q~~hat(mu)[2] == e4~~hat(sigma)[2] == e5, 
    list(q = round(1 - estimate[1]-estimate[6],3), 
       e4 = round(estimate[4],1),
       e5 = round(estimate[5],1))))
       
    text(x = max(sizeKb)/2, y = 0.85*ymax, col="green3", pos = 4, cex=1.2,
    substitute(hat(p)[3] == p3~~hat(mu)[3] == e7~~hat(sigma)[3] == e8, 
    list(p3 = round(estimate[6],3), 
       e7 = round(estimate[7],1),
       e8 = round(estimate[8],1))))     
       
    dst <- density(sizeKb)
    lines(x=dst$x, y=vscale*dst$y)
    legend(x = max(sizeKb)/2, y = 0.75*ymax, lty=1,
       "Gaussian kernel density estimation")
       
    mtext(sodtxt)
  }
#
# Return invisibly the dataset
#
  invisible(data)
}













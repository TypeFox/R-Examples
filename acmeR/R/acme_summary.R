#' Report Summary of ACME Dataset
#' 
#' Provide various summaries for an event-level dataset, including estaimtes
#' and plots.
#' 
#' This function takes in an event-level dataset of placement and searches,
#' and reports various statistics related to ACME. This includes empirical 
#' information such as search intervals, carcass information summaries, and
#' available scavenger information. This function also reports on estimated
#' values R* and bleed-through rate. All information is printed to the console.
#' 
#' Users have the option of creating plots of empirical data overlayed on estimated 
#' distributions for both persistance and search efficiency.
#' 
#'
#'@param fname Data, either a string for csv files or a data frame name 
#'@param spec Species subset. Default (empty string) includes all species in 
#'data set.
#'@param blind Logical. If TRUE, assumes FT are unaware of carcasses.
#'@param plot_scav Logical. If TRUE, persistence plot is output. Default is FALSE.
#'@param plot_srch Logical. If TRUE, search efficiency plot is output. Default 
#' is FALSE.
#'@param ps Character string, name of Postscript file generated. Default (empty
#' string) displays plots on console and does not generate postscript file.
#'@return \code{acme.summary} returns a list with the following components:
#' \item{params}{5-element vector: alpha and rho parameters for the 
#' Weibull persistence distribution, a and b parameters for the 
#' exponentially decreasing search proficiency, and bt as the bleed-through
#' parameter.}
#' \item{Rstar}{Reduction factor (inverse of the inflation factor for mortality
#' counts).}
#' \item{T}{First component of R* - carcass fraction found on first search.}
#' \item{I}{Average interval length, in days.}
#' 
#'@import graphics grDevices
#'@export 
#'@examples
#'\dontrun{
#'#If altamont is a file in the working directory
#'acme.summary('altamont.csv', spec = "BHCO")
#'
#'#To include plots
#'acme.summary('altamont.csv',spec = "BHCO", plot_scav = TRUE, plot_srch = TRUE)
#'}

##################################################################
# acme.summary():  Report summary of ACME data (-> UI)
#
acme.summary <- function(fname, spec="", blind=TRUE, ps="", plot_scav = FALSE,
                         plot_srch = FALSE)
{
  if(missing(fname)){
    make.csv(altamont_internal,fname = 'altamont.csv')
    stop("Missing csv file 'fname'.  Example 'altamont.csv' created in
     working directory.  Try acme.summary('altamont.csv').");
  }
  rd  <- read.data(fname, spec=spec, blind=blind);
  
  name <- deparse(substitute(fname))
  cat(paste("Data set \'", name, "\' includes records of ",
              ncarc <- dim(rd$scav)[1], " placed carcasses\nof ",
              nsp <- length(unique(rd$scav$Species)),
              " species:\n", sep=""));
  tab <- table(rd$scav$Species);
  num <- sort(tab, decreasing=T);
  spc <- names(num);
  nsn <- nam <- character(nsp);
  for(i in 1:nsp) {
    nam[i] <- str2str(spec2name(spc[i]),22);
    nsn[i] <- paste(nam[i], "(", spc[i],":", int2str(num[i],3),")",sep="");
  }
  SpecPerRow <- 2;
  for(i in seq(1,nsp,SpecPerRow)) {
    cat("  ");
    cat(paste(paste(nsn[i:min(i+(SpecPerRow-1),nsp)],
                    sep="",collapse=",  "),"\n"));
  }
  NP.num <- length(rd$NP.ID);
  if(NP.num == 1) {
    cat(paste("In addition, 1 carcass of species \'",
              spec2name(rd$NP.Spec), "\' (", rd$NP.Spec,
              ") was\ndiscovered but not placed (and ",
              "is not used in this analysis).\n",sep=""));
  } else if(NP.num > 1) {
    cat(paste("\nIn addition, ", NP.num, " carcasses ",
              "of ", NP.Nsp<-length(rd$NP.Spec), " species were ",
              "discovered\nbut not placed (and are not used in",
              " this analysis):\n",sep=""));
    for(i in seq(1, NP.Nsp, 6)) {
      cat("    ");
      cat(paste(rd$NP.Spec[i:min(i+5,NP.Nsp)],
                sep="",collapse=", "),"\n");
    }
  }
  Day   <- 86400;            # Seconds per day: 60*60*24
  Range <- format(range(rd$scav[,"Placed"]), '%a, %b %d, %Y');
  fnd.pct <- sse(rd);
  fnd.cnt <- round(fnd.pct * ncarc);
  cat(paste("\nField technicians discovered ", fnd.cnt, " of these ",
            ncarc, " placed carcasses\n(", round(100*fnd.pct,2),sep="")); 
  cat(paste("%) on at least one occasion.  Carcass placements ranged\n",
            "from ", Range[1], " to ", Range[2], ".\n\n", sep=""));
  cat(paste("FT  Search intervals averaged ",
            round(rd$Ik["mu"],2), " days, with sd = ",
            round(rd$Ik["sd"],2), " days.\n", sep=""));
  cat(paste("PFM Check  intervals averaged ",
            round(rd$Sk["mu"],2), " days, with sd = ",
            round(rd$Sk["sd"],2), " days.\n", sep=""));
  if(!blind) cat(paste("Unblinded FT searches (those after the first",
            " discovery) were\ntreated as PFM checks.\n", sep=""));
  scav <- mle.wei(rd, spec=spec, v=FALSE);
  cat(paste("\nThe estimated scavenger removal survival function is:\n\n",
            "    P[ remain > t days ] = exp(-(", round(scav$rho,4),"*t)^",
            round(scav$alp,4),
            "),\n\n(Weibull distribution) with mean removal time of ",
            round(scav$tij,2)," days.\n",sep=""));
  se <- abs(scav$alp-1)/scav$alp.se;
  if     (se>-qnorm(0.001)) { exp.assump <- "completely implausible"; }
  else if(se>-qnorm(0.010)) { exp.assump <- "rather implausible"; }
  else if(se>-qnorm(0.050)) { exp.assump <- "implausible"; }
  else if(se>-qnorm(0.100)) { exp.assump <- "dubious but plausible"; }
  else          { exp.assump <- "plausible"; }
  cat(paste("Weibull shape parameter alpha=",round(scav$alp,4)," is ",
            round(se <- abs(scav$alp-1)/scav$alp.se,1),
            " standard errors away\nfrom 1, ",
            "making the Exponential Distribution ",
            exp.assump, "\n(P-value < ", signif(2*pnorm(se,lower.tail=FALSE),2),
            ").\n", sep=""));
  if(nchar(ps)>0) {
    #    postscript(file=ps, paper="special", height=8, width=12, hori=FALSE);
    postscript(file=ps, horizontal=TRUE);
  }
  #Plotting, if switches are on;
  if(plot_scav | plot_srch){
    opar <- par(no.readonly=TRUE);
    par(mar=opar$mar+c(1,1,0,0), mfrow=c(2,1));
    if(plot_scav) scav.plot(rd, spec=spec, add=TRUE);
    if(plot_srch) srch.plot(rd, spec=spec, add=TRUE, verb=FALSE);
    par(opar)
    
    if(nchar(ps)>0) {
      dev.off();
      print(paste("Plot stored as postscript file",ps));
    }
  }
  #Obtaining parameters
  nv    <- naive.srch(rd, spec);
  srch  <- mle.srch(rd, spec=spec, v=FALSE);
  arabt <- c(alp=scav$alp, rho=scav$rho,
             a=srch$a.hat, b=srch$b.hat, bt=srch$bt.hat);  
  nbleed <- sum(bleed(rd));
  if(nbleed>0) bleed <- "occurred"
  else bleed <- "was not observed";
  cat(paste("\nField Technicians tried ", nv$nsrch, " times to find ",
            "carcasses known to be present,\nand succeeded ", nv$succ[3],
            " times (", round(100*nv$succ[3]/nv$nsrch,2),
            "% of the time).  Of the total number ",  nv$tot, 
            " of\ncarcasses placed, ", nv$tot-nv$ncarc,
            " were removed before the first FT search.  Of the\n",
            "remaining ", nv$ncarc, " carcasses ever present during ",
            "a search, ", nv$succ[1], " were ever\ndiscovered (",
            round(100*nv$succ[1]/nv$ncarc,2), "%), ", nv$succ[2],
            " of those (", round(100*nv$succ[2]/nv$nsrch,2),
            "%) on the first try.\n",
            "The estimated FT search proficiency is\n\n    ",
            "P[ discovery | carcass age = t days ] = exp(-",
            round(srch$a.hat,4),"-",round(srch$b.hat,4), "*t),\n\n",
            "starting at ", round(100*exp(-srch$a.hat),2), 
            "% for fresh carcasses",sep=""));
  #  Does S(t) fall over 5% in 100 days? -log1p(-0.05)/100=0.000513
  if(srch$b.hat > 0.0005) cat(paste(", falling to ",
            round(100*exp(-srch$a.hat-srch$b.hat*rd$Ik["mu"]),2), 
            "% after one\n", round(rd$Ik["mu"],1), 
            "-day FT search interval, ",
            round(100*exp(-srch$a.hat-2*srch$b.hat*rd$Ik["mu"]),2), 
            "% after two, and ",
            round(100*exp(-srch$a.hat-3*srch$b.hat*rd$Ik["mu"]),2), 
            "% after three.",sep=""))
  else cat(paste(", and remaining nearly constant\nafter that.", sep=""));
  I     <- round(rd$Ik["mu"]);
  Rest <- Rst(I, arabt);
  Rstar <- Rest$R;
  T <- Rest$T0;
  cat(paste("\n\nCarcasses were later discovered after an initial miss ",
            nv$succ[1]-nv$succ[2], " times, and\nafter some earlier",
            " miss ", nbleed," times, suggesting bleed-through ",
            bleed, ".\nEstimated bleed-through rate is ",
            round(100*srch$bt.hat,2),"%.\n",
            "\nACME mortality estimates for ", I, "-day search intervals ",
            "that discover C\ncarcasses are:\n",
            "          M* = ", round(1/Rstar,2)," * C = C/R*, for ",
            "R* = ", round(Rstar,4),".\n", sep=""));
  return(list(params = arabt, Rstar = Rstar, T = T, I = I))
}
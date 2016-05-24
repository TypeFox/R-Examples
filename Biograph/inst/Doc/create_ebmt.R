

rm(list=ls())
library (Biograph)
library (mstate)
data (ebmt4)
attach(ebmt4)
#  2279 patients transplanted at the EBMT between 1985 and 1998;15 variables
 # id Patient identification number
 # rec Time in days from transplantation to recovery or last follow-up
 # rec.s Recovery status; 1 = recovery, 0 = censored
 # ae Time in days from transplantation to adverse event (AE) or last follow-up
 # ae.s Adverse event status; 1 = adverse event, 0 = censored
 # recae Time in days from transplantation to both recovery and AE or last follow-up
 # plag.s Recovery and AE status; 1 = both recovery and AE, 0 = no recovery or no AE or censored
 # rel Time in days from transplantation to relapse or last follow-up
 # rel.s Relapse status; 1 = relapse, 0 = censored
 # srv Time in days from transplantation to death or last follow-up
 # srv.s Relapse status; 1 = dead, 0 = censored
 # year Year of transplantation; factor with levels "1985-1989", "1990-1994", "1995-1998"
#  agecl Patient age at transplant; factor with levels "<=20", "20-40", ">40"
#  proph Prophylaxis; factor with levels "no", "yes"
#  match Donor-recipient gender match; factor with levels "no gender mismatch", "gender mismatch"

namstates <- c("T","P","A","Z","R","D")
nsample <- nrow(ebmt4)
# --------   Get dates at ENTRIES into states (remove censoring)   --------
rec2 <- ifelse (rec.s==1,rec,NA)       # recovery
ae2 <- ifelse (ae.s==1,ae,NA)          # adverse event (AE)
ae2 <- ifelse (!is.na(rec2)&!is.na(ae2)&ae2>rec2,NA,ae2)
rec2 <- ifelse (!is.na(rec2)&!is.na(ae2)&rec2>ae2,NA,rec2) 
recae2 <- ifelse (recae.s==1,recae,NA) # determine recovery after an adverse event
rel2 <- ifelse (rel.s==1,rel,NA)        # relapse
srv2 <- ifelse (srv.s==1,srv,NA)        # death
# 305 experience relapse before death, but transition Rel-> Death is not feasible (see transition matrix)
# Hence number of deaths is 533
srv2 <- ifelse (!is.na(rel2)&!is.na(srv2)&rel2<srv2,NA,srv2)
# --------   status variable: event or censoring   -------------
censored <- 1-srv.s   # Censored: not death before end of study
day.cens <- ifelse(censored==1,srv,NA) # censoring: duration since transplant
# relapse = end of observation (censoring)
day.cens <- ifelse (rel.s==1,NA,day.cens)
#    end of observation (follow-up): death (censored=0) or censoring (censored=1)
#end <-srv 
end <- ifelse (rel.s==1,rel2,srv)  # censored at relapse

# -------   construct data frame with dates at transition   ---------
days <- data.frame(P=rec2,A=ae2,Z=recae2,R=rel2,D=srv2)

# ------------ Arrange transitions in chronological order  ----------
f<- Sequences.ind.0 (d=days,namstates=namstates,absorb=c("R","D"))
path <- as.character(f$path)

# ------------  covariates  -------------------
cov <- c("year","agecl","proph","match")
attach (ebmt4)
options(stringsAsFactors=FALSE)
# ------------  Create Biograph object "EBMT"   -------------
namtrans <- paste("Tr",1:3,sep="")
EBMT <- data.frame (ID=ebmt4$id,born=rep(0,nsample),start=rep(0,nsample),end=end,year=year,agecl=agecl,proph=proph,match=match,path=as.character(path),f$d[,1:(max(nchar(path))-1)])
namcov <- c("year","agecl","proph","match")
colnames(EBMT) <- c("ID","born","start","end",namcov,"path",namtrans)
locpat <- locpath(EBMT)
attr(EBMT,"format.date") <- "day"
attr(EBMT,"format.born") <- "day"
param <- Parameters(EBMT)
attr(EBMT,"param") <- param
 zz8 <- "/Users/frans/Documents/R/0 0 MAC/Package/Biograph.TEST/Chapters/AnnexA/EBMT/"
setwd(zz8)
save (EBMT,file="EBMT.RData")

 
mkSEERold<-function(df,seerHome="~/data/SEER",dataset=c("00","73","92"), SQL = TRUE, mkDFs=FALSE){
  #	if (!require(LaF)) print("R Package LaF must be installed")
  print("mkSEERold is deprecated. It only works with the current release, i.e. of April 2016")
  colTypes=c("integer","string",rep("integer",5),"double")    # double, integer, categorical and string
  colWidths=c(4,7,2,1,1,1,2, 10) 
  colNames = c('popyear','X0','popseer','poprace','origin','popsex', 'popage','population')
  seerHome=path.expand(seerHome)
  dataset=match.arg(dataset)
  dir.create(file.path(seerHome, dataset),showWarnings = FALSE)
  cat("\n\n\nProcessing dataset",dataset," \n")
  ptm <- proc.time()
  for (i in c("19agegroups","singleages")) {
    f=switch(dataset,
             "00" = file.path(seerHome,paste0("populations/expanded.race.by.hispanic/yr2000_2013.ca_ky_lo_nj_ga/",i,".txt")),
             "73" = file.path(seerHome,paste0("populations/white_black_other/yr1973_2013.seer9/",i,".txt")),
             "92" = file.path(seerHome,paste0("populations/expanded.race.by.hispanic/yr1992_2013.seer9.plus.sj_la_rg_ak/",i,".txt")))
    laf<-laf_open_fwf(f,column_types=colTypes,column_widths=colWidths,column_names = colNames)
    if (i=="19agegroups") {
      pops=laf[,colNames[-c(2,5)]] 
      if (mkDFs) save(pops,file=fp<-file.path(seerHome, dataset,"pops.RData"))  }
    else {
      popsa=laf[,colNames[-c(2,5)]]
      if (mkDFs) save(popsa,file=fp<-file.path(seerHome, dataset,"popsa.RData"))  
    }
  }
  
  #   s=switch(dataset,
  #            "00" = file.path(seerHome,"populations/expanded.race.by.hispanic/yr2000_2011.ca_ky_lo_nj_ga/singleages.txt"),
  #            "73" = file.path(seerHome,"populations/white_black_other/yr1973_2011.seer9/singleages.txt"),
  #            "92" = file.path(seerHome,"populations/expanded.race.by.hispanic/yr1992_2011.seer9.plus.sj_la_rg_ak/singleages.txt"))
  #   ptm <- proc.time()
  #   laf<-laf_open_fwf(s,column_types=colTypes,column_widths=colWidths,column_names = colNames)
  #   popsa<-laf[,colNames[-c(2,5)]]
  #   save(popsa,file=fp<-file.path(seerHome, dataset,"popsa.RData"))  
  #   delT=proc.time() - ptm  
  #   cat("The single ages population file of SEER dataset ",dataset," was successfully written to:\n ",fp," in ",delT[3]," seconds.\n",sep="")
  
  if (SQL) {
    cat("Deleting old and creating new SQLite database all.db and adding popultation tables to it.\n")
    dbf<-file.path(seerHome, dataset,"all.db")
    unlink(dbf)
    m <- dbDriver("SQLite")
    con <- dbConnect(m, dbname = dbf)
    dbWriteTable(con, "pops", pops,overwrite=TRUE)
    dbWriteTable(con, "popsa", popsa,overwrite=TRUE)
  }	
  delT=proc.time() - ptm  
  cat("The population files of SEER dataset ",dataset," were processed in ",delT[3]," seconds.\n",sep="")
  
  
  
  ptm <- proc.time()
  cancers=c('breast','digothr','malegen','femgen','other','respir','colrect','lymyleuk','urinary','test') 
  p=switch(dataset,
           "00" = file.path(seerHome,"incidence/yr2000_2013.ca_ky_lo_nj_ga/"),
           "73" = file.path(seerHome,"incidence/yr1973_2013.seer9/"),
           "92" = file.path(seerHome,"incidence/yr1992_2013.sj_la_rg_ak/"))
  y=df[which(df$names!=" "),"names"]; cat("Cancer Data:\nThe following fields will be written:\n");	print(y)
  for (k in 1:9)	{
    laf<-laf_open_fwf(file.path(p,paste0(toupper(cancers[k]),'.TXT')), 
                      column_types=df$type,    # double, integer, categorical and string
                      column_widths=df$width)
    DF=laf[,which(df$names!=" ")]
    colnames(DF)<-y
    DF=mapCancs(DF)
    if (SQL) dbWriteTable(con, cancers[k], DF,overwrite=TRUE)
    if (mkDFs) save(DF,file=fp<-file.path(seerHome,dataset,paste0(cancers[k],'.RData')))
    cat(k,") ",cancers[k],"\n",sep="")
  }
  delT=proc.time() - ptm  
  cat("Cancer files of SEER ",dataset," were written to ",file.path(seerHome, dataset)," in ",delT[3]," seconds.\n")
  if (SQL) {
    cat("\nThe following tables are in ",file.path(seerHome, dataset,"all.db"),":\n")
    print(dbListTables(con)) 	
    dbDisconnect(con)
  }
#   seerHome="~/data/SEER";dataset="92"
#if(.Platform$OS.type=="unix") {
  cat("\nFiles in ",file.path(seerHome, dataset),"are:\n")
  s=dir(file.path(seerHome, dataset),full.names=T)
  d=file.info(s)[,c("size","mtime")] 
  rws=rownames(d)
  dbs=sum(d$size[grep("db",rws)])
  dfs=sum(d$size[grep("RData",rws)])
  d$size=paste0(round(d$size/1e6,1)," MB")
  print(d) 
  cat("\nThe RData files sum to ",round(dfs/1e6,1)," MB and the db files sum to ",
      round(dbs/1e6,1)," MB, i.e. the fraction is ",round(dfs/dbs,2),".\n",sep="")
#} 
}

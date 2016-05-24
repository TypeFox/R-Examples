readHMD <- function(username, password,
  country=c(AUS='Australia', AUT='Austria', BLR='Belarus',
            BEL='Belgium', BGR='Bulgaria', CAN='Canada',
            CHL='Chile', CZE='Czech Republic', DNK='Denmark',
            EST='Estonia', FIN='Finland', FRA='France',
            DEU='Germany', HUN='Hungary', ISL='Iceland',
            IRL='Ireland', ISR='Israel', ITA='Italy',
            JPN='Japan', LVA='Latvia', LTU='Lithuania',
            LUX='Luxemburg', NDL='Netherlands', NZL='New Zealand',
            NOR='Norway', POL='Poland', PRT='Portugal',
            RUS='Russia', SVK='Slovakia', SVN='Slovenia',
            ESP='Spain', SWE='Sweden', CHE='Switzerland',
            TWN='Taiwan', GBR='U.K.', USA='U.S.A.',
            UKR='Ukraine'),
  sex=c('m', 'f', 'b'), HMDurl='http://www.mortality.org/hmd',
  dataType = 'lt',
  ltCol=c('m', 'q', 'a', 'l', 'd', 'L', 'T', 'e'),
  cohper = c(coh='cohort', per='periodic'), ageInterval=c(1, 5),
  yearInterval=c(1, 5, 10), url, ... )
{
##
## 1.  Construct the desired url
##
  require(RCurl)
  missURL <- missing(url)
  if(missURL){
      if(length(country)!=1)
          stop('country must be specified')
      if(nchar(country)==3){
          Ctry <- country
      } else {
          Ctries=c(AUS='Australia', AUT='Austria', BLR='Belarus',
            BEL='Belgium', BGR='Bulgaria', CAN='Canada',
            CHL='Chile', CZE='Czech Republic', DNK='Denmark',
            EST='Estonia', FIN='Finland', FRA='France',
            DEU='Germany', HUN='Hungary', ISL='Iceland',
            IRL='Ireland', ISR='Israel', ITA='Italy',
            JPN='Japan', LVA='Latvia', LTU='Lithuania',
            LUX='Luxemburg', NDL='Netherlands', NZL='New Zealand',
            NOR='Norway', POL='Poland', PRT='Portugal',
            RUS='Russia', SVK='Slovakia', SVN='Slovenia',
            ESP='Spain', SWE='Sweden', CHE='Switzerland',
            TWN='Taiwan', GBR='U.K.', USA='U.S.A.',
            UKR='Ukraine')
          ct <- match.arg(country)
          sel <- (Ctries==ct)
          if(sum(sel)!=1)
              stop('Error in the country name.')
          Ctry <- names(Ctries)[sel]
      }
      sx <- match.arg(sex)
      if(length(cohper)>1){
          pc <- 'coh'
      } else {
          if(nchar(cohper)==3){
              pc <- cohper
          } else {
              pc. <- match.arg(cohper)
              sel <- (c('cohort', 'periodic') == pc.)
              pc <- c('coh', 'per')[sel]
          }
      }
      if(length(ageInterval)>1){
          AI <- 1
      } else {
          sel <- (ageInterval %in% c(1, 5))
          if(sum(sel)!=1)
              stop('ageInterval must be 1 or 5; is ',
                   ageInterval)
          AI <- ageInterval[sel]
      }
      if(length(yearInterval)>1){
          YI <- 1
      } else {
          sel <- (yearInterval %in% c(1, 5, 10))
          if(sum(sel)!=1)
              stop('yearInterval must be 1, 5 or 10;  is ',
                   yearInterval)
          YI <- yearInterval[sel]
      }
      url <- paste(HMDurl, '/', Ctry, '/STATS/', sx, dataType, pc, '_',
                   AI, 'x', YI, '.txt',  sep='')
  }
##
## 2.  connection
##
  userpwd <- paste(username,":",password,sep="")
  txt <- RCurl::getURL(url,userpwd=userpwd, ...)
  conChk <- textConnection(txt)
##
## 3.  Read
##
  HMD <- readLines(conChk)
  if(length(HMD)<11){
      print(HMD)
      close(conChk)
      stop('No data table found at ', url)
  }
  close(conChk)
#
  con <- textConnection(txt)
  hmd <- try(read.table(con,skip=2,header=T,na.strings=".",
                    stringsAsFactors=FALSE))
  close(con)
  if((class(hmd)=='try-error') || !missURL ||
     (dataType != 'lt') ){
      if(class(hmd) == 'try-error')
          warning('Error in read.table trapped by "try"; ',
                  ' returning what was read')
      return(list(URL=txt, getURL=txt, readLines=HMD,
                  read.table=hmd))
  }
  j <- hmd[,"Year"]==hmd[1,"Year"]
  x <- as.numeric(gsub("\\+","",as.character(hmd[j,"Age"])))
##
## 4.  Extract life table data
##
#  if(missURL && (dataType=='lt')){
#  col <- match(tolower(sx),tolower(colnames(hmd)))
  col <- match.arg(ltCol)
#      ltCols=c('m', 'q', 'a', 'l', 'd', 'L', 'T', 'e')
  Col <- pmatch(col, names(hmd))
#  y <- matrix(as.numeric(hmd[,col]),nrow=length(x))
  y <- matrix(hmd[, Col], nrow=length(x))
  colnames(y) <- unique(hmd[,"Year"])
  rownames(y) <- as.character(hmd[j,"Age"])
#  if(substr(file,1,2)=="Mx")
  yname <- c(m='Mortality rate', q='Mortality probability',
             a='Survival time for mortalities',
             l='Number of survivors',
             d='Number of deaths',
             L='Person-years in interval',
             T='Person-years remaining',
             e='Life expectancy')[col]
#
  return(structure(list(x=x,y=y,time=sort(unique(hmd[,"Year"])),
                        xname="Age",yname=yname),class=c("fts","fds")))
#  } else {
#      return(hmd)
#  }
}

ci5zip = 'http://ci5.iarc.fr/CI5plus/old/CI5plus.zip'
canZip = 'http://www20.statcan.gc.ca/tables-tableaux/cansim/csv/01030550-eng.zip'

cancerRates = function(area = "canada",
   year=2000, sex=c("M", "F"), site="Lung") {

  sexes = c("M"=1, "F"=2)[sex]



  areaCodes =
  c("Canada"=1240099,
  "Newfoundland"=1240899,
"Prince Edward Island"=1240799,
"Nova Scotia"=1240699,
"Ontario"=1241199,
"Manitoba"=1240399,
"Saskatchewan"=1241399,
"Alberta"=1240199,
"British Columbia"=1240299,
"New Zeland"=5540099,
"Sweden"=7520099,
"Slovenia"=7050099,
"Slovakia"=7030099,
Norway=5780099,
Latvia=4280099,
Lithuania=4400099,
Iceland=3520099,
Finland=2460099,
Estonia=2330099,
Denmark=2080099,
"Czech Republic"=2030099,
"Costa Rica"=1880099,
USA=8400199,
Iowa=8400899,
"New Mexico"=8401399
)
if(is.character(area)) {
area = areaCodes[grep(paste("^", area[1], sep=""), 
				names(areaCodes), ignore.case=TRUE)]
} 


result = list()

rates=NULL
for(Dsex in names(sexes)) {
fs<-paste("http://ci5.iarc.fr/CI5plus/old/Table4r.asp?registry=",area,(paste("&period=",year,sep="",collapse="")),"&sex=", sexes[Dsex],"&window=1&text=1&stat=0&submit=Execute",sep="")
tempn = scan(fs, what="a", quiet=TRUE)
theurl=(paste("http://ci5.iarc.fr/", 
					gsub("^HREF=", "", grep("href=/data", 
									iconv(tempn,to="UTF-8"), value=TRUE, 
									ignore.case=TRUE)), sep=""))
theurl = url(theurl)
forAttribute = scan(theurl, what="a", sep="\t", nmax=1, quiet=TRUE)
result[[Dsex]]=read.table(theurl, header=TRUE,skip=1, 
		fill=TRUE, sep="\t", as.is=TRUE)
}

iarcSite=NULL
for (Dsex in names(sexes)) {
 cancerTypes =  result[[Dsex]][,1]
  
 cancerTypes = gsub("[[:space:]]+$", "", cancerTypes)

siterow = grep(site, cancerTypes,ignore.case=TRUE)
if(length(siterow) > 1 )    {
  warning(paste("matched ", paste(cancerTypes[siterow], collapse=","), "\n from cancer types", paste(cancerTypes, collapse=","), 
  "\n using", cancerTypes[siterow[1]]))
  siterow=siterow[1]
  }
if(length(siterow) ) {
  x=result[[Dsex]][siterow,]
iarcSite =gsub("[[:space:]]+$", "", x[1]) 
  x = x[,-c(1,grep("Total|Unk|ICD|CR|ASR", 
						  names(x), ignore.case=TRUE))]
  x = as.vector(x)
  names(x) = gsub("^X|\\.$", "", names(x))
  rates[[Dsex]]=x
  }
}
rates = unlist(rates)/100000
names(rates) = gsub("\\.", "_", names(rates))

attributes(rates)$site = iarcSite
forAttribute = strsplit(forAttribute, "\\(")[[1]]
attributes(rates)$area = gsub("[[:space:]]+$", "", forAttribute[1])
attributes(rates)$year = gsub("\\)$", "", forAttribute[length(forAttribute)])

return(rates)

}

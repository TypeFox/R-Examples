
CompactPhysicalConstantNames<-function(s)
{
# compact the physical constant names to camelCase
# Input: string with the physical constant names'
# Output: string with names in camelCase
R<-gsub('\\(.*?\\)','',s)
R<-gsub('\\[.*?\\]','',R)
R<-gsub(',',' ',R)
R<-gsub('-',' ',R)
R<-gsub('/','',R)
R<-gsub('\\.',' ',R)
R<-gsub(' +',' ',R)
R<-gsub('\u00e5','a',R)
R<-gsub('\u00f6','o',R)
R<-gsub(' ([a-z])','\\U\\1',R,  perl = TRUE)
R<-gsub(' ','',R)
R
}

CompactFactorNames<-function(s)
{
# compact the factor names to camelCase
# Input: string with the factor names'
# Output: string with names in camelCase
R<-gsub('\\(.*?\\)','',s)
R<-gsub('\\[.*?\\]','',R)
R<-gsub(',',' ',R)
R<-gsub('-',' ',R)
R<-gsub('/','',R)
R<-gsub(' +',' ',R)
R<-gsub('\u00e5','a',R)
R<-gsub('\u00f6','o',R)
R<-gsub(' ([a-z])','\\U\\1',R,  perl = TRUE)
R<-gsub(' ','',R)
R
}

ShortenTitles<-function(s){
R<-gsub('\\(.*?\\)','',s)
R<-gsub('\\[.*?\\]','',R)
R<-gsub(' +',' ',R)
R
}

ShortenLongNames<-function(s)
{
R<-gsub('square','sqr',s, perl=TRUE, ignore.case=TRUE)
R<-gsub('second','sec',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('minute','min',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('kilogram','kg',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('foot','ft',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('milliliter','ml',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('millimiter','mm',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('kilometer','km',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('calorie','cal',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('fahrenheit','f',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('degree','deg',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('british','uk',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('celsius','c',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('centimeter','cm',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('acceleration','accel',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('standard','std',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('conventional','convtnl',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('kilopascal','kpascal',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('thermal','th',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('unit','un',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('kelvin','k',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('to ?the ?fourth','4th',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('cubic','cub',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('megapascal','mpascal',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('metric','SI',R, perl=TRUE, ignore.case=TRUE)
R
}

FixFactorNames<-function(s)
{
R<-gsub('\u00b0','',s, perl=TRUE, ignore.case=TRUE)
R<-gsub('\u00e5','a',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('\u00f6','o',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('\u00b7','*',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('\u03a9','ohm',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('\u03bc','u',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('\u03b3','y',R, perl=TRUE, ignore.case=TRUE)
R<-gsub('\u2032',"'",R, perl=TRUE, ignore.case=TRUE)
R
}

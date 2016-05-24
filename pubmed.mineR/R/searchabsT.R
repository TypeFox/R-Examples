setGeneric("searchabsT", function(object,yr, include, restrict, exclude) standardGeneric("searchabsT")); 
setMethod("searchabsT","Abstracts", function(object,yr, include, restrict, exclude){  

if(missing(yr) == TRUE){yr = "NONE"}; if(missing(include) == TRUE){include = "NONE"}; if(missing(restrict) == TRUE){restrict = "NONE"}; if(missing(exclude) == TRUE){exclude = "NONE"};
if (yr[1] != "NONE") {

for (j in seq(along = include )) { if (j == 1) resy1 <- getabsT(object,yr[j],FALSE) else if (j == 2) {resy2 <- getabsT(object,yr[j],FALSE); resy <- combineabs(resy1,resy2)} else if ( j > 2) resy <- combineabs(resy,getabsT(object,yr[j],FALSE));
if (length(yr) == 1) resy <- resy1}}
else resy = object;


if (include[1] != "NONE"){ 
            {for (j in seq(along = include)) {if ( j == 1 ) {reso1 <- getabs(resy,include[j],FALSE)}  
else if ( j == 2 )  {reso2 <- getabs(resy,include[j],FALSE); reso <- combineabs(reso1,reso2)} 
else if ( j > 2) {temp1 <- getabs(resy,include[j],FALSE); reso <- combineabs(reso, temp1)}}}; 
if (length(include) == 1) {reso <- reso1}}
else reso = resy;


if ( restrict[1] != "NONE") {resz = reso;
            for (j in seq(along = restrict)) { resz <- getabs(resz,restrict[j],FALSE)}}
else {resz = reso};



if (exclude[1] != "NONE") {resn <- resz;   
            
for (j in seq(along = exclude)) { resn <- removeabs(resn,exclude[j],FALSE) }}
else {resn = resz}; 
##sendabs(resn, "out.txt");
return(resn)}
)

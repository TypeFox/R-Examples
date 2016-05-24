
# 070605: MS modified to insert NA when valid OS reference could not be found.
 
".selectstring" <- function (text, first, last = 1e+06){
    storage.mode(text) <- "character"
    n <- max(lt <- length(text), length(first), length(last))
    if (lt < n) 
        text <- rep(text, length = n)
    substr(text, first, last)
}

"epi.convgrid" <- function(os.refs){

name <- c("SW","SX","SY","SZ","TV","SV","SQ","SR","SS","ST","SU","TQ","TR","SM","SN","SO","SP","TL","TM","NQ","NL","NF","NA","HV","HQ","HL","SG","SB","NW",
"NR","NM","NG","NB","HW","HR","HM","SH","SC","NX","NS","NN","NH","NC","HX","HS","HN","SJ","SD","NY","NT","NO","NJ","ND","HY","HT","HO","SK","SE",
"NZ","NU","NP","NK","NE","HZ","HU","HP","TF","TA","OV","OQ","OL","OF","OA","JV","JQ","JL","TG","TB","OW")

easting <- c(100,200,300,400,500,0,0,100,200,300,400,500,600,100,200,300,400,500,600,0,0,0,
0,0,0,0,100,100,100,100,100,100,100,100,100,100,200,200,200,200,200,200,200,200,
200,200,300,300,300,300,300,300,300,300,300,300,400,400,400,400,400,400,400,400,400,400,
500,500,500,500,500,500,500,500,500,500,600,600,600)

northing <- c(0,0,0,0,0,0,100,100,100,100,100,100,100,200,200,200,200,
200,200,600,700,800,900,1000,1100,1200,300,400,500,600,700,800,900,1000,
1100,1200,300,400,500,600,700,800,900,1000,1100,1200,300,400,500,600,700,
800,900,1000,1100,1200,300,400,500,600,700,800,900,1000,1100,1200,300,400,
500,600,700,800,900,1000,1100,1200,300,400,500)

cells <- as.data.frame(cbind(easting, northing))
cells <- as.data.frame(cbind(name, cells))

x.coord = 0
y.coord = 0
res <- as.matrix(cbind(x.coord, y.coord))

for(i in 1:length(os.refs)){

grid <- os.refs[i]
grid <- .selectstring(grid, 1, 2)

coords <- os.refs[i]
easting <- as.numeric(.selectstring(coords, 3, 5)) * 100
northing <- as.numeric(.selectstring(coords, 6, 8)) * 100
id <- cells$name == grid
tmp <- cells[id, 1:3]
tmp <- cbind(((tmp$easting * 1000) + easting), ((tmp$northing * 1000) + northing))

   if(dim(tmp)[1] == 0) tmp <- matrix(c(NA, NA), nrow = 1)

rval <- rbind(res, tmp)
   }
rval <- rval[-1,]
return(rval)

}

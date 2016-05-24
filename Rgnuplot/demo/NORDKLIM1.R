
data(NordklimData)
# nordklim <- read.table(system.file('extdata/NordklimData.tab', package = 'Rgnuplot'), stringsAsFactors=FALSE, header=TRUE)

NKmonths <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
# choose Helsinki (code 304) and country (code 'FIN') Precipitation (code 601)
nordklimHelsinkiPrecipitation <- NordklimData[which((NordklimData$NordklimNumber == 304) & (NordklimData$CountryCode == "FIN") & (NordklimData$ClimateElement == 601)), c("FirstYear", 
    NKmonths)]
nordklimHelsinkiPrecipitation <- as.matrix(nordklimHelsinkiPrecipitation)
# save to a data file
GpMatrixr2gnu(nordklimHelsinkiPrecipitation, "NORDKLIM-Helsinki-prec.dat")
# 2D plot with yearly average of data OK!
Gprun("#set terminal png;set output \"NORDKLIM1-1.png\"\nset title \"Helsinki Precipitation\"\nplot \"NORDKLIM-Helsinki-prec.dat\" using 1:(($2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13)/12.0) with lines tit \"yearly average\"", 
    TRUE)

# 2D plot with 4 months of data OK!
Gprun("#set terminal png;set output \"NORDKLIM1-2.png\"\nset title \"Helsinki Precipitation\"\nplot \"NORDKLIM-Helsinki-prec.dat\" using 1:2 with lines tit \"January\",  \"NORDKLIM-Helsinki-prec.dat\" using 1:3 with lines tit \"February\",\\\n\"NORDKLIM-Helsinki-prec.dat\" using 1:4 with lines tit \"March\",  \"NORDKLIM-Helsinki-prec.dat\" using 1:5 with lines tit \"April\"", 
    TRUE)

# convert the data file from 12 columns for the monthly data to 2 rows
GpCols2rows("NORDKLIM-Helsinki-prec.dat", "NORDKLIM-Helsinki-prec-columns.dat")
# 3D all months/years lines
Gprun("unset key\n#unset ztics\nset grid xtics\nset grid ytics\nset ztics scale 5\nset xmtics\nset pm3d\nset xrange[1:12]\nsplot \"NORDKLIM-Helsinki-prec-columns.dat\" using 1:2:3  w l  linewidth 2 linecolor 3", 
    TRUE)

# 2D heatmap
Gprun("#set terminal png;set output \"NORDKLIM1-4.png\"\nset xrange [ 1891 : 2002 ]\nset yrange [ 1 : 12 ] # hide the year column\nset ymtics\nset pm3d map\nset palette defined (0 \"black\", 2000 \"white\")\nset title \"Helsinki Precipitation\"\nsplot \"NORDKLIM-Helsinki-prec.dat\" u ($2+1891):1:3 matrix w pm3d notit", 
    TRUE)


# 3D all months/years
Gprun("#set terminal png;set output \"NORDKLIM1-4.png\"\nset xrange [ 1891 : 2002 ]\nset yrange [ 1 : 12 ] # hide the year column\nset zrange [ 0:2000 ]\nset ymtics\nset palette defined (0 \"black\", 2000 \"white\")\nset title \"Helsinki Precipitation\"\nsplot \"NORDKLIM-Helsinki-prec.dat\" u ($2+1891):1:3 matrix w pm3d notit", 
    TRUE)

GpShowPaletteColornames()
 

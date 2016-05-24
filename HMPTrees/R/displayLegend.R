displayLegend <-
function(colors, divisions, title="Confidence Value"){
if(missing(divisions)) 
divisions <- c(0, .1, 1, 10, 100, 1000, 10000)

divisions <- sort(divisions)
if(missing(colors)) 
colors <- c("red", "orange", "yellow", "green" , "cyan", "blue")

if(length(divisions) > (length(colors)+1)) #need more colors, dont care if more colors than divisons
colors <- c(colors, rep(colors[length(colors)], (length(divisions) - length(colors)-1)))

lgd <- NULL
for(num in length(divisions):2)
lgd <- c(lgd, paste(divisions[num], "-", divisions[num-1], sep=""))

palette(colors)
plot(0, 0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
legend(-.75, .75, legend=lgd, col=rev(palette()), pch=19, title=title)
}

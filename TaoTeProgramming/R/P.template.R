"P.template" <-
function (filename="xxx.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=3, width=5)
                par(mar=c(5,4,0,2) + .1, cex.axis=.7, cex.lab=.7)
        } 

        if(length(filename)) dev.off()
}


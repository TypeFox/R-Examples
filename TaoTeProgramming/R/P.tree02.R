"P.tree02" <-
function (filename="tree02.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=6.5, width=4)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	tree(branches=30, seed=12)

        if(length(filename)) dev.off()
}


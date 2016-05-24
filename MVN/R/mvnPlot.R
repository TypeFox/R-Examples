mvnPlot <-
function(object, type=c("persp","contour"), default=TRUE, ...){
        
        type <- match.arg(type)
        if (!(class(object)[1]  %in% c("mardia", "hz", "royston"))) stop("Object must be in one of the following classes: \"mardia\", \"hz\", \"royston\" ")
        p <- ncol(object@dataframe)
        if (p != 2) stop("Plots are available for bivariate normal distributions. Number of variables exceed 2.")
        
        dataframe <- object@dataframe
        
        data.kde <- kde2d(dataframe[,1], dataframe[,2], n=100)
        
        if (type == "persp"){
            if (default){
                persp(data.kde, theta = 1, phi = 30, border = NA, shade = 0.5,
                      box=T, xlab = colnames(dataframe)[1], ylab = colnames(dataframe)[2],
                      main = "", zlab = "Density")
            } 
            else {
                persp(data.kde, ...)
            }
        }
        
        if (type == "contour"){
            if (default){
                contour(data.kde, nlevels=20, xlab = colnames(dataframe)[1], ylab = colnames(dataframe)[2])
            }
            else {
                contour(data.kde, ...)
            }
        }
    }

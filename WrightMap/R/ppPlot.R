ppPlot <-
function(thetas, thresholds = NULL, est, SE, main.title = "Person Probability Plot", cut.left = 0, cut.right = .94, cut.lab.adj = c(1,.5),...) {
	
	prob.calc <- function( x, theta.est, pr = .5 ){

        (exp( 1 * ( theta.est - x ) ) / ( 1 + exp( 1 * ( theta.est - x ) )) - pr)^2

    }
    
    cut.at <- function(x) {optimize( prob.calc, interval = c(-10,10), theta.est = est, pr = x)$minimum}
    
    cutpoints <- sapply(c(.2,.4,.5,.6,.8),cut.at)
    
    
    wrightMap(thetas, thresholds, main.title = main.title, person.points = est, person.range = c(est - SE, est + SE), cutpoints = cutpoints, cut.lab.text = c("20%","40%","50%","60%","80%"), cut.left = cut.left, cut.right = cut.right, cut.lab.adj = cut.lab.adj,...)



    
    
}

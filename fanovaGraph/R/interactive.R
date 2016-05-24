plotTk <- function(graphlist, delta.layout = 0.01) {
    # requireNamespace(tcltk) # included in namespace
    d <- graphlist$d
    dall <- graphlist$V
    totalInt <- graphlist$tii[,1]
    tii.layout <- threshold(graphlist, delta = delta.layout, scaled = TRUE)$tii[,1]
    E.layout <- t(combn(d,2)[,tii.layout>0])
    g.layout <- graph(as.vector(t(E.layout)), n = d, directed = FALSE)
    layout <- layout.fruchterman.reingold(g.layout)
    max.delta <- max(totalInt/dall + 0.001)
    
    plotting <- function(delta) {
        graph <- threshold(graphlist, delta=delta, scaled = TRUE)
        plot.graphlist(graph, layout = layout)
        n.CL <- length(graph$cliques)
        title(main = paste("delta =", round(delta, 5)))
        title(sub = paste("number of cliques =", n.CL))
    }
    
    variable <- tclVar(delta.layout)  # start value
    refresh <- function(...) {
        # function for every change of input
        delta <- as.numeric(tclvalue(variable))
        plotting(delta)
    }
    m <- tktoplevel()  # Erstellen des Eingabefensters
    tkwm.title(m, "input window")
    cutFrame <- tkframe(m)
    cutSlider <- tkscale(cutFrame, command = refresh, from = 0, to = max.delta, 
        resolution = 5e-04, orient = "horiz", variable = variable)
    tkpack(tklabel(cutFrame, text = "Delta:"), side = "left")
    tkpack(cutFrame, cutSlider, side = "bottom")
}

# necessary for checks:
if(getRversion() >= "2.15.1")  globalVariables(c("delta", "manipulate", "slider")) 

plotManipulate <- function(graphlist, delta.layout = 0.01) {
    requireNamespace(manipulate)
    d <- graphlist$d
    dall <- graphlist$V
    totalInt <- graphlist$tii[,1]
    tii.layout <- threshold(graphlist, delta = delta.layout, scaled = TRUE)$tii[,1]
    E.layout <- t(combn(d,2)[,tii.layout>0])
    g.layout <- graph(as.vector(t(E.layout)) , n = d, directed = FALSE)
    layout <- layout.fruchterman.reingold(g.layout)
    max.delta <- max(totalInt/dall + 0.001)
    
    plotting <- function(delta) {
        graph <- threshold(graphlist, delta=delta, scaled = TRUE)
        plot.graphlist(graph, layout = layout)
        n.CL <- length(graph$cliques)
        title(main = paste("delta =", round(delta, 5)))
        title(sub = paste("number of cliques =", n.CL))
    }
    manipulate(plotting(delta), delta = slider(0, max.delta, initial = delta.layout, 
        step = 5e-04))
} 
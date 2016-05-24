## This defines a class of hypothesis trees. It should
## allow us to manipulate the hypotheses, the results of
## tests associated with them, and their hierarchical
## relation. We set a class to make this manpiulation
## easier and to prevent the bugs that arise from complexity.
## Note that the tree is defined by its adjacency matrix.
setClass("hypothesesTree", representation = list(tree = "matrix",
                               p.vals = "data.frame",
                               alpha = "numeric"))

setMethod("initialize", "hypothesesTree", function(.Object, ...) {
          value <- callNextMethod()
          value
      })

setMethod("show", "hypothesesTree", function(object)
      {
          tree <- slot(object, "tree")
          p.vals <- slot(object, "p.vals")
          alpha <- slot(object, "alpha")

          cat("hFDR adjusted p-values:", "\n")
          print(p.vals)
          cat('---', '\n')          
          cat('Signif. codes:  0 \'***\'', alpha / 50, '\'**\'', alpha / 5, '\'*\'', alpha, '\'.\'', 2 * alpha, '\'-\' 1', '\n')          
      })

setMethod("summary", "hypothesesTree", function(object) {
          tree <- slot(object, "tree")
          p.vals <- slot(object, "p.vals")
          alpha <- slot(object, "alpha")

          igraph.tree <- graph.edgelist(tree)
          nHypotheses <- length(V(igraph.tree))
          nEdges <- length(E(igraph.tree))
          cat("Number of hypotheses:", nHypotheses, "\n")
          FDR.control <- EstimatedHFDRControl(object)
          cat("Number of tree discoveries:", FDR.control$n.tree.discoveries, "\n")
          cat("Estimated tree FDR:", FDR.control$tree, "\n")
          cat("Number of tip discoveries:", FDR.control$n.tip.discoveries, "\n")
          cat("Estimated tips FDR:", FDR.control$tip, "\n")
          cat("\n", "hFDR adjusted p-values:", "\n")
          n.to.print <- min(nrow(p.vals), 10)
          print(p.vals[order(p.vals$adjp)[1:n.to.print], ])
          if(n.to.print < nrow(object@p.vals)) {
              cat('[only 10 most significant hypotheses shown]', '\n')
          }
          cat('---', '\n')
          cat('Signif. codes:  0 \'***\'', alpha / 50, '\'**\'',
              alpha / 5, '\'*\'', alpha, '\'.\'',
              2 * alpha, '\'-\' 1', '\n')          
      })

setMethod("plot", "hypothesesTree",
          function(x,..., adjust = TRUE,
                   return_script = FALSE, width = 900,
                   height = 500, base_font_size = 12,
                   output_file_name = paste('hyp_tree', gsub("[^\\d]+", "", Sys.time(), perl=TRUE), '.html', sep = "")) {
              PlotHypTree(x, adjust, return_script, width,
                          height, base_font_size, output_file_name)
          })

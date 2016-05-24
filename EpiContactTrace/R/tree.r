## Copyright 2013-2014 Stefan Widgren and Maria Noremark,
## National Veterinary Institute, Sweden
##
## Licensed under the EUPL, Version 1.1 or - as soon they
## will be approved by the European Commission - subsequent
## versions of the EUPL (the "Licence");
## You may not use this work except in compliance with the
## Licence.
## You may obtain a copy of the Licence at:
##
## http://ec.europa.eu/idabc/eupl
##
## Unless required by applicable law or agreed to in
## writing, software distributed under the Licence is
## distributed on an "AS IS" basis,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
## express or implied.
## See the Licence for the specific language governing
## permissions and limitations under the Licence.

##' Build a graph tree from the NetworkStructure
##'
##' @param network_structure a data.frame from the call
##' \code{NetworkStructure} with a \code{ContactTrace} object
##' @return A \code{list} with the two fields \code{ingoing} and
##' \code{outgoing}. The fields are \code{NULL} or contain a
##' \code{data.frame} with the tree. The fields are \code{NULL} if
##' there are no in- or outgoing contacts.
##' @keywords internal
build_tree <- function(network_structure)
{
    stopifnot(is.data.frame(network_structure))
    root <- unique(network_structure$root)
    stopifnot(identical(length(root), 1L))

    tree.in <- network_structure[network_structure$direction == "in",]
    tree.out <- network_structure[network_structure$direction == "out",]

    result <- list(ingoing=NULL, outgoing=NULL)

    root_node <- data.frame(node = root[1],
                            parent = NA_character_,
                            level = 0,
                            stringsAsFactors = FALSE)

    if(nrow(tree.in)) {
        i <- order(tree.in$distance, tree.in$source)
        tree.in <- tree.in[i, c('source', 'distance')]
        tree.in <- tree.in[!duplicated(tree.in$source),]
        tree.in$parent <- NA_character_
        colnames(tree.in)[1:2] <- c('node', 'level')
        tree.in <- tree.in[, colnames(root_node)]

        for(lev in rev(seq_len(max(tree.in$level)))) {
            for(src in tree.in$node[tree.in$level == lev]) {
                if(lev > 1) {
                    i <- which(network_structure$direction == "in"
                               & network_structure$distance == lev)
                    dst <- network_structure$destination[i]
                    dst <- unique(dst)
                } else {
                    dst <- root
                }

                stopifnot(length(dst)>0)
                tree.in$parent[tree.in$level == lev
                               & tree.in$node == src] <- dst[1]
            }
        }

        tree.in <- rbind(root_node, tree.in)
        rownames(tree.in) <- NULL
        result$ingoing <- tree.in
    }

    if(nrow(tree.out)) {
        i <- order(tree.out$distance, tree.out$destination)
        tree.out <- tree.out[i, c('destination', 'distance')]
        tree.out <- tree.out[!duplicated(tree.out$destination),]
        tree.out$parent <- NA_character_
        colnames(tree.out)[1:2] <- c('node', 'level')
        tree.out <- tree.out[, colnames(root_node)]

        for(lev in rev(seq_len(max(tree.out$level)))) {
            for(dst in tree.out$node[tree.out$level == lev]) {
                if(lev > 1) {
                    i <- which(network_structure$direction == "out"
                               & network_structure$distance == lev)
                    src <- network_structure$source[i]
                    src <- unique(src)
                } else {
                    src <- root
                }

                stopifnot(length(src)>0)
                tree.out$parent[tree.out$level == lev
                                & tree.out$node == dst] <- src[1]
            }
        }

        tree.out <- rbind(root_node, tree.out)
        rownames(tree.out) <- NULL
        result$outgoing <- tree.out
    }

    return(result)
}

##' Position nodes in a tree
##'
##' This function determines the coordinates for each node in a
##' tree.
##' @param tree
##' @param x The x coordinate of the root node
##' @param y The y coordinate of the root node
##' @param orientation The orientation of the tree. \code{North}, the
##' root is at the top. \code{South}, the root is at the
##' bottom. \code{East}, the root is at the left. \code{West}, th root
##' is at the right.
##' @param sibling_separation The minimum distance between adjacent
##' siblings of the tree
##' @param subtree_separation The minimum distance between adjacent
##' subtrees of a tree.
##' @param level_separation The fixed distance between adjacent levels
##' of the tree.
##' @param left_size The left size of a node.
##' @param right_size The right size of a node.
##' @param top_size  The top size of a node.
##' @param bottom_size The bottom size of a node.
##' @keywords internal
##' @references \itemize{
##'   \item John Q. Walker II, A node positioning algorithm for general tress.\cr
##'   \url{http://www.cs.unc.edu/techreports/89-034.pdf}
##'}
position_tree <- function(tree,
                          x = 0,
                          y = 0,
                          orientation = c("North", "South", "East", "West"),
                          sibling_separation = 4,
                          subtree_separation = 4,
                          level_separation = 1,
                          left_size = 1,
                          right_size = 1,
                          top_size = 1,
                          bottom_size = 1)
{
    ## Clean up the positioning of small sibling subtrees
    apportion <- function(node) {
        left_most <- first_child(node)
        neighbor <- left_neighbor(left_most)
        compare_depth <- 1L
        depth_to_stop <- max_depth() - node_level(node)

        while(all(!is.null(left_most),
                  !is.null(neighbor),
                  compare_depth <= depth_to_stop)) {
            ## Compute the location of left_most and where it
            ## should be with respect to neighbor.
            left_modsum <- 0
            right_modsum <- 0
            ancestor_left_most <- left_most
            ancestor_neighbor <- neighbor

            for(i in seq_len(compare_depth)) {
                ancestor_left_most <- parent(ancestor_left_most)
                ancestor_neighbor <- parent(ancestor_neighbor)
                right_modsum <- right_modsum + modifier(ancestor_left_most)
                left_modsum <- left_modsum + modifier(ancestor_neighbor)
            }

            ## Find the move_distance, and apply it to Node's subtree.
            ## Add appropriate portions to smaller interior subtrees.
            move_distance <- ((prelim(neighbor) +
                               left_modsum +
                               subtree_separation +
                               mean_node_size(left_most, neighbor)) -
                              (prelim(left_most) + right_modsum))

            if(move_distance > 0) {
                ## Count interior sibling subtrees in left siblings
                temp_node <- node
                left_siblings <- 0

                while(all(!is.null(temp_node),
                          !identical(temp_node, ancestor_neighbor))) {
                    left_siblings <- left_siblings + 1
                    temp_node <- left_sibling(temp_node)
                }

                if(!is.null(temp_node)) {
                    ## Apply portions to appropriate leftsibling
                    ## subtrees
                    portion <- move_distance / left_siblings
                    temp_node <- node

                    while(!identical(temp_node, ancestor_neighbor)) {
                        set_prelim(temp_node, prelim(temp_node) + move_distance)
                        set_modifier(temp_node, modifier(temp_node) + move_distance)
                        move_distance <- move_distance - portion
                        temp_node <- left_sibling(temp_node)
                    }
                } else {
                    return(NULL)
                }
            }

            compare_depth <- compare_depth + 1L
            if(is_leaf(left_most)) {
                left_most <- get_left_most(node, 0L, compare_depth)
            } else {
                left_most <- first_child(left_most)
            }
            neighbor <- left_neighbor(left_most)
        }

        return(NULL)
    }

    ##
    ## Help functions to work with nodes
    ##
    is_leaf <- function(node) {
        return(!has_child(node))
    }

    first_child <- function(node) {
        children <- tree$node[!is.na(tree$parent) & (tree$parent == node)]
        if(length(children)>0) {
            return(children[1])
        }
        return(NULL)
    }

    get_left_most <- function(node, level, depth) {
        if(level >= depth) {
            return(node)
        } else if(is_leaf(node)) {
            return(NULL)
        } else {
            right_most <- first_child(node)
            left_most <- get_left_most(right_most, level + 1L, depth)

            while(all(is.null(left_most),
                      has_right_sibling(right_most))) {
                right_most <- right_sibling(right_most)
                left_most <- get_left_most(right_most, level + 1L, depth)
            }

            return(left_most)
        }
    }

    has_child <- function(node) {
        return(!is.null(first_child(node)))
    }

    node_index <- function(node) {
        ## Deterimine row index to the node
        i <- which(tree$node == node)
        stopifnot(identical(length(i), 1L))
        return(i[1])
    }

    node_level <- function(node) {
        return(tree$level[node_index(node)])
    }

    left_neighbor <- function(node) {
        if(!is.null(node)) {
            n <- tree$node[tree$level == node_level(node)]
            stopifnot(node %in% n)
            i <- which(node == n)
            stopifnot(identical(length(i), 1L))
            if(i[1] > 1)
                return(n[i[1]-1])
        }
        return(NULL)
    }

    parent <- function(node) {
        p <- tree$parent[node_index(node)]
        stopifnot(identical(is.na(p), FALSE))
        return(p[1])
    }

    root <- function() {
        return(tree$node[1])
    }

    max_depth <- function() {
        return(max(tree$level))
    }

    ##
    ## Help functions to work with the size of nodes
    ##
    mean_node_size <- function(left_node, right_node) {
        node_size <- 0

        if(any(identical(orientation, "North"),
               identical(orientation, "South"))) {
            if(!is.null(left_node))
                node_size <- node_size + get_right_size(left_node)
            if(!is.null(right_node))
                node_size <- node_size + get_left_size(right_node)
        } else if(any(identical(orientation, "East"),
                      identical(orientation, "West"))) {
            if(!is.null(left_node))
                node_size <- node_size + get_top_size(left_node)
            if(!is.null(right_node))
                node_size <- node_size + get_bottom_size(right_node)
        }

        return(node_size)
    }

    get_left_size <- function(node) {
        return(tree$left_size[node_index(node)])
    }

    get_right_size <- function(node) {
        return(tree$right_size[node_index(node)])
    }

    get_top_size <- function(node) {
        return(tree$top_size[node_index(node)])
    }

    get_bottom_size <- function(node) {
        return(tree$bottom_size[node_index(node)])
    }

    ##
    ## Help functions to work with coordinates
    ##
    prelim <- function(node) {
        return(tree$prelim[node_index(node)])
    }

    set_prelim <- function(node, prelim) {
        tree$prelim[node_index(node)] <<- prelim
    }

    modifier <- function(node) {
        return(tree$modifier[node_index(node)])
    }

    set_modifier <- function(node, modifier) {
        tree$modifier[node_index(node)] <<- modifier
    }

    set_x <- function(node, x) {
        tree$x[node_index(node)] <<- x
    }

    set_y <- function(node, y) {
        tree$y[node_index(node)] <<- y
    }

    ##
    ## Help functions to work with siblings
    ##
    siblings <- function(node) {
        i <- node_index(node)
        parent <- tree$parent[i]
        if(is.na(parent)) {
            ## Check that node is root
            stopifnot(identical(node_level(node), 0L))
            siblings <- node
        } else {
            siblings <- tree$node[!is.na(tree$parent) &
                                  (tree$parent == parent)]
        }
        stopifnot(node %in% siblings)
        return(siblings)
    }

    has_left_sibling <- function(node) {
        return(!is.null(left_sibling(node)))
    }

    has_right_sibling <- function(node) {
        return(!is.null(right_sibling(node)))
    }

    left_sibling <- function(node) {
        s <- siblings(node)
        i <- which(node == s)
        stopifnot(identical(length(i), 1L))
        if(i[1] > 1)
            return(s[i[1]-1])
        return(NULL)
    }

    right_sibling <- function(node) {
        s <- siblings(node)
        i <- which(node == s)
        stopifnot(identical(length(i), 1L))
        if(i[1] < length(s))
            return(s[i[1]+1])
        return(NULL)
    }

    ## Every node of the tree is assigned a preliminary x-coordinate
    ## (held in column prelim). In addition, internal nodes are given
    ## modifiers, which will be used to move their offspring to the
    ## right (held in column modifier).
    first_walk <- function(node) {
        ## Set the default modifier value.
        set_modifier(node, 0)

        if(any(is_leaf(node),
               node_level(node) == max_depth())) {
            if(has_left_sibling(node)) {
                ## Determine the preliminary x-coordinate based on:
                ##  - the preliminary x-coordinate of the left sibling,
                ##  - the separation between sibling nodes, and
                ##  - the mean size of left sibling and current node.
                set_prelim(node,
                           prelim(left_sibling(node)) +
                           sibling_separation +
                           mean_node_size(left_sibling(node), node))
            } else {
                ## No sibling on the left to worry about.
                set_prelim(node, 0)
            }
        } else {
            ## This node is not a leaf, so call this procedure
            ## recursively for each of its offspring.
            right_most <- first_child(node)
            left_most <- right_most
            first_walk(left_most)
            while(has_right_sibling(right_most)) {
                right_most <- right_sibling(right_most)
                first_walk(right_most)
            }

            mid_point <- (prelim(left_most) + prelim(right_most)) / 2

            if(has_left_sibling(node)) {
                set_prelim(node,
                           prelim(left_sibling(node)) +
                           sibling_separation +
                           mean_node_size(left_sibling(node), node))

                set_modifier(node, prelim(node) - mid_point)

                apportion(node)
            } else {
                set_prelim(node, mid_point)
            }
        }
    }

    check_extents_range <- function(x_temp, y_temp) {
        return(TRUE)
    }

    ## Each node is given a final x-coordinate by summing its
    ## preliminary x-coordinate and the modifiers of all the node's
    ## ancestors. The y-coordinate depends on the height of the
    ## tree. If the actual position of an interior node is right of
    ## its preliminary place, the subtree rooted at the node must be
    ## moved right to center the sons around the father. Rather than
    ## immediately readjust all the nodes in the subtree, each node
    ## remembers the distance to the provisional place in a modifier
    ## field. In this second pass down the tree, modifiers are
    ## accumulated and applied to every node.
    second_walk <- function(node, modsum) {
        if(node_level(node) <= max_depth()) {
            if(identical(orientation, "North")) {
                x_temp <- x + prelim(node) + modsum
                y_temp <- y - node_level(node) * level_separation
            } else if(identical(orientation, "South")) {
                x_temp <- x + prelim(node) + modsum
                y_temp <- y + node_level(node) * level_separation
            } else if(identical(orientation, "East")) {
                x_temp <- x - node_level(node) * level_separation
                y_temp <- y + prelim(node) + modsum
            } else if(identical(orientation, "West")) {
                x_temp <- x + node_level(node) * level_separation
                y_temp <- y + prelim(node) + modsum
            } else {
                stop('Undefined orientation')
            }

            ## Check that x_temp and y_temp are of the proper size.
            if(check_extents_range(x_temp, y_temp)) {
                set_x(node, x_temp)
                set_y(node, y_temp)

                if(has_child(node)) {
                    ## Apply the modifier value for this node to
                    ## all its offspring.
                    second_walk(first_child(node), modsum + modifier(node))
                }

                if(has_right_sibling(node)) {
                    second_walk(right_sibling(node), modsum)
                }
            } else {
                stop('Tree outside drawable extents range')
            }
        }

        return(NULL)
    }

    orientation <- match.arg(orientation)
    tree$level <- as.integer(tree$level)
    tree$x <- NA_real_
    tree$y <- NA_real_
    tree$prelim <- NA_real_
    tree$modifier <- NA_real_
    tree$left_size <- left_size
    tree$right_size <- right_size
    tree$top_size <- top_size
    tree$bottom_size <- bottom_size

    first_walk(root())

    if(any(identical(orientation, "North"),
           identical(orientation, "South"))) {
        x <- x - prelim(root())
    } else if(any(identical(orientation, "East"),
                  identical(orientation, "West"))) {
        y <- y - prelim(root())
    }

    second_walk(root(), 0)

    return(tree[, c('node', 'parent', 'level', 'x', 'y')])
}

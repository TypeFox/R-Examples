#' Decision Tree Plot
#'
#' This function will draw decision tree based on the tree information list. All the plotting functions 
#' come frome rpart package.
#'
#' @param mytree a list returned by the decTreexxx function contains the information of the decision tree.
#' @param width, etc. parameters for ploting, see rpart package for more details. You can try serveral times
#' for the best display.
#' @param height height
#' @param uniform TRUE or FALSE
#' @param branch branch
#' @param margin margin
#' @param nspace nspace
#' @param minbranch minbranch
#' @param all all
#' @param fancy fancy
#' @param use.n use.n
#' @param fwidth fwidth
#' @param fheight fheight
#' @return the plot will be save as a png file in your working directory.
#' @export
#' @author Li Qinglong <liqinglong0830@@163.com>
#' @seealso plotDecROC
plotDecTree <- function(mytree, 
                        width = 1000, 
                        height = 750, 
                        uniform = TRUE, 
                        branch = 0, 
                        margin = 0.1, 
                        nspace = -1, 
                        minbranch = 0.3, 
                        all = FALSE, 
                        fancy = TRUE, 
                        use.n = TRUE, 
                        fwidth = 0.5, 
                        fheight = 0.5)
{
    # We only need to know mytree[[3]]
    treeNodeInfo = mytree[[3]]

    # 1 Generate frame
    frame = data.frame(matrix(nrow = length(treeNodeInfo), ncol = 7))
    colnames(frame) = c("row.names", "var", "n", "dev", "yval", "Split_var_id", "Split_value")
    for (row in 1:nrow(frame))
    {
        frame[row, 1] = convertNodeId(treeNodeInfo[[row]]$Node_id)
        if (treeNodeInfo[[row]]$isLeaf)
            frame[row, 2] = "<leaf>"
        else
            frame[row, 2] = treeNodeInfo[[row]]$Split_var_name
        frame[row, 3] = treeNodeInfo[[row]]$Size
        frame[row, 5] = treeNodeInfo[[row]]$predicted_rank
        frame[row, 6] = treeNodeInfo[[row]]$Split_var_id
        frame[row, 7] = treeNodeInfo[[row]]$Split_value
    }
    for (row in 1:nrow(frame))
    {
        if (frame[row, 2] == "<leaf>")
            frame[row, 4] = frame[row, 3]
        else
            frame[row, 4] = frame[frame[, 1] == frame[row, 1] * 2 + 1, 3]
    }
    # turn frame in to pre ordr
    frame = frame[visitNodePreorder(frame[, 1]), ]

    #2 generate split and csplits
    str_var_info = mytree$Var_Info
    # sperate by "\t"
    str_var_info = strsplit(str_var_info[6:length(str_var_info)], "\t")
    # var_info contains variable information, ncat means n categories
    var_info = data.frame(matrix(nrow = length(str_var_info), ncol = 3))
    colnames(var_info) = c("name", "type", "ncat")
    for (i in 1:length(str_var_info))
    {
        var_info$name[i] = str_var_info[[i]][1]
        var_info$type[i] = str_var_info[[i]][2]
        var_info$ncat[i] = as.numeric(str_var_info[[i]][3])
    }

    internal_node_rows = which(frame[,2] != "<leaf>")
    splits = data.frame(matrix(nrow = sum(frame[,2] != "<leaf>"), ncol = 4))
    colnames(splits) = c("row.names", "count",    "ncat", "index")
    splits$row.names = frame$var[internal_node_rows]
    splits$count = frame$n[internal_node_rows]
    csplits = matrix(2, ncol = max(var_info$ncat[var_info$type == 'd']),
                     nrow = sum(var_info$type[frame$Split_var_id[internal_node_rows]+1] == "d"))
    index = 0
    for (row in 1:nrow(splits))
    {
        var_ind = frame$Split_var_id[internal_node_rows[row]] + 1
        if (var_info$type[var_ind] == "c")
        {
            # continuous variable, index is the spliting value, ncat = -1
            splits$ncat[row] = -1
            splits$index[row] = frame$Split_value[internal_node_rows[row]]
        } else {
            # discrete variable,index means the index-th row of csplit, ncat is the size of classes
            # csplit, 1 < left, 3 > right, 2 unknown at the last
            # Split_value is a vector, value of 1 means to be classified to left, otherwise 3
            splits$ncat[row] = var_info$ncat[var_ind]
            index = index + 1
            splits$index[row] = index
            csplits[index, 1:splits$ncat[row]] = 3
            csplits[index, frame$Split_value[internal_node_rows[row]]] = 1
        }
    }

    # 3 build xlevels
    var_levels = list()
    for (i in 1:length(str_var_info))
    {
        if (str_var_info[[i]][2] == "d")
        {
            name = str_var_info[[i]][1]
            choice_size = as.numeric(str_var_info[[i]][4])
            lvls = list(str_var_info[[i]][5:(4+choice_size)])
            names(lvls) = name
            var_levels = c(var_levels, lvls)
        }
    }
    # 4 create tree obj
    x1 = frame[, 2:5]
    rownames(x1) = frame$row.names
    x2 = as.matrix(splits[, 2:4])
    rownames(x2) = splits$row.names
    mat = as.matrix(splits[, 2:4])
    x3 = csplits
    colnames(x3) = NULL
    obj = list(frame = x1, splits = x2, csplit = x3)
    attr(obj, "xlevels") = var_levels
    # 5 binomial tree plot
    png("mytree.png", width = 1000, height = 750)
    myplot.rpart(obj, uniform = T, branch = 0, margin = 0.1)
    parms = list(uniform = TRUE, branch = 0, nspace = -1, minbranch = 0.3)
    mytext.rpart(obj, parms = parms, all = F, fancy = T, use.n = T, fwidth = 0.5, fheight = 0.5, minlength = 1)
    dev.off()    
}

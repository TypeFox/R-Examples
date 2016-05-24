rpart.predict.leaves <- function (rp, newdata, type = "where") 
{
#
# With a tip of the hat to Yuji at
# http://stackoverflow.com/questions/5102754/
# search-for-corresponding-node-in-a-regression-tree-using-rpart
#
if (type == "where")
    rp$frame$yval <- 1:nrow(rp$frame)
else if (type == "leaf")
    rp$frame$yval <- rownames(rp$frame)
    else
    stop ("Type must be 'where' or 'leaf'")
return (predict(rp, newdata=newdata, type = "vector"))
}


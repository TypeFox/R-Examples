`check` <-
function (x, v, l1, l2) 
{
    if (v[l1, l2] > 0) {
        print("WARNING from routine ", x, ": element", l1, l2, 
            "already non-zero")
    }
    v[l1, l2] <- x
    v
}


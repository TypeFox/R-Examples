muprob <-
function (cumprob, nobs, ncategoriesm1) 
{
    ans <- matrix(cumprob, nobs, ncategoriesm1, TRUE)
    ans <- rbind(ans[, 1L], diff(t(ans)), 1 - ans[, ncategoriesm1])
    ans <- c(ans)
    ans
}


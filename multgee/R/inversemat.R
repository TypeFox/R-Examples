inversemat <-
function (x, IM) 
{
    switch(IM, cholesky = chol2inv(chol(x)), solve = solve(x, diagmod(rep.int(1, nrow(x)))),
        qr.solve = qr.solve(x, diagmod(rep.int(1, nrow(x)))))
}

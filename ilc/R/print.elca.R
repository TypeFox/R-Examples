print.elca <-
function(x, ...){
    cat('\n', paste(rep('-', 60), collapse=''))
    cat('\n\t Extended Lee-Carter Regression:')
    cat('\n\t Fitted Model:', x$model)
    cat('\n', paste(rep('-', 60), collapse=''))
    cat('\n Call:', deparse(x$call), fill=T)
    cat(' Error Structure:', x$adjust)
    cat('\n Data Source:', x$label, 'over')
    cat('\n   calendar years:', rb(str.range(x$year)), ', ages:',
        rb(str.range(x$age)))
    cat('\n   and groups:', names(x$ag))
    cat('\n Deviance convergence in:', x$conv.iter, 'iterations\n')
    print(structure(list(dev=x$mdev, df=x$df), class='coef'))
}

fix0202 <-
function (mu = 0, sig = 1, eta = 0, kappa = 1, alf = 2, bet = 1, 
    mdl = -1) 
{
    cenprob <- pnorm(mdl, eta, kappa)
    gam11.exp <- fix0203(gamobs.11, mdl, eta + 10 * kappa, mu, 
        sig, eta, kappa, alf, bet) + gamcen.11(mdl, mu, sig, 
        eta, kappa, alf, bet) * cenprob
    gam12.exp <- fix0203(gamobs.12, mdl, eta + 10 * kappa, mu, 
        sig, eta, kappa, alf, bet) + gamcen.12(mdl, mu, sig, 
        eta, kappa, alf, bet) * cenprob
    gam22.exp <- fix0203(gamobs.22, mdl, eta + 10 * kappa, mu, 
        sig, eta, kappa, alf, bet) + gamcen.22(mdl, mu, sig, 
        eta, kappa, alf, bet) * cenprob
    delmat <- matrix(c(gam11.exp, gam12.exp, gam12.exp, gam22.exp), 
        2, 2)
    delmat.inv <- solve(delmat)
    psi11.exp <- fix0203(psiobs.1sq, mdl, eta + 10 * kappa, mu, 
        sig, eta, kappa, alf, bet) + psicen.1sq(mdl, mu, sig, 
        eta, kappa, alf, bet) * cenprob
    psi22.exp <- fix0203(psiobs.2sq, mdl, eta + 10 * kappa, mu, 
        sig, eta, kappa, alf, bet) + psicen.2sq(mdl, mu, sig, 
        eta, kappa, alf, bet) * cenprob
    psi12.exp <- fix0203(psiobs.12, mdl, eta + 10 * kappa, mu, 
        sig, eta, kappa, alf, bet) + psicen.12(mdl, mu, sig, 
        eta, kappa, alf, bet) * cenprob
    covmat <- matrix(c(psi11.exp, psi12.exp, psi12.exp, psi22.exp), 
        2, 2)
    asympt.var <- delmat.inv %*% covmat %*% solve(t(delmat))
    params <- c(mu = mu, sig = sig, eta = eta, kappa = kappa, 
        alf = alf, bet = bet, mdl = mdl)
    list(params = params, asympt.var = asympt.var, psi.var = covmat, 
        gam.matinv = delmat.inv, gam.mat = delmat)
}

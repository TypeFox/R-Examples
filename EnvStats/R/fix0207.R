fix0207 <-
function (mu = 0, sig = 1, eta = 0, kappa = 1, alf = 2, bet = 1, 
    mdl = -1) 
{
    r <- (mdl - mu)/sig
    w <- (mdl - eta)/kappa
    s.r <- s.0.m.singly.censored(r, alf, bet)
    phi.w <- dnorm(w)
    int.psi <- fix0203(etaobs.1, mdl, eta + 10 * kappa, mu, sig, 
        eta, kappa, alf, bet)
    int.psi2 <- fix0203(etaobs.2, mdl, eta + 10 * kappa, mu, 
        sig, eta, kappa, alf, bet)
    int.chi <- fix0203(chiobs.1, mdl, eta + 10 * kappa, mu, sig, 
        eta, kappa, alf, bet)
    int.chi2 <- fix0203(chiobs.2, mdl, eta + 10 * kappa, mu, 
        sig, eta, kappa, alf, bet)
    m.11 <- int.psi - phi.w * s.r
    m.12 <- int.psi2 - w * phi.w * s.r
    m.21 <- (1 - r * s.r) * phi.w + int.chi
    m.22 <- w * phi.w * (1 - r * s.r) + int.chi2
    (-1/(kappa * sig)) * matrix(c(m.11, m.12, m.21, m.22), 2, 
        2, byrow = TRUE)
}

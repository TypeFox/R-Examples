s <- summary(orings.model)
sqrt(diag(s$cov.unscaled)) -> st.err; st.err
coef(orings.model)[2] + c(-1,1) * st.err[2] * qnorm(0.975)
exp(coef(orings.model)[2] + c(-1,1) * st.err[2] * qnorm(0.975))

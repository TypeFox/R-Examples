`Q_function_addendo2` <-
function(C0, m_0, P_0_n, y_0_n) {
        Q_addendo2 <- log(det(C0)) + sum(diag( solve(C0) %*% (P_0_n + (y_0_n - m_0) %*% t((y_0_n - m_0))) ))
        return(Q_addendo2)
}


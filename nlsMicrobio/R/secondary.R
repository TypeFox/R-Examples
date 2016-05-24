"cpm_T" <- as.formula(sqrtmumax ~ sqrt(((T >= Tmin) & (T <= Tmax)) * muopt * (T - Tmax) * (T - Tmin)^2 / ((Topt - Tmin) * ((Topt - Tmin) * (T - Topt) - (Topt - Tmax) * (Topt + Tmin - 2.0 * T)))))
"cpm_pH_4p" <- as.formula(sqrtmumax ~ sqrt(((pH >= pHmin) & (pH <= pHmax)) * muopt * (pH - pHmin) * (pH - pHmax) / ((pH - pHmin) * (pH - pHmax) - (pH - pHopt)^2)))
"cpm_pH_3p" <- as.formula(sqrtmumax ~ sqrt(((pH >= pHmin) & (pH <= (2 * pHopt - pHmin))) * muopt * (pH - pHmin) * (pH - (2 * pHopt - pHmin)) / ((pH - pHmin) * (pH - (2 * pHopt - pHmin)) - (pH - pHopt)^2)))
"cpm_aw_3p" <- as.formula(sqrtmumax ~ sqrt((aw >= awmin) * muopt * (aw - 1) * (aw - awmin)^2 / ((awopt - awmin) * ((awopt - awmin) * (aw - awopt) - (awopt - 1) * (awopt + awmin - 2.0 * aw)))))
"cpm_aw_2p" <- as.formula(sqrtmumax ~ sqrt((aw >= awmin) * muopt * (aw - awmin)^2 / (1 - awmin)^2))
"cpm_T_pH_aw" <- as.formula(sqrtmumax ~ sqrt(((T >= Tmin) & (T <= Tmax) & (pH >= pHmin) & (pH <= (pHmax)) & (aw >= awmin)) * muopt * (T-Tmax) * (T - Tmin)^2 / ((Topt - Tmin) * ((Topt - Tmin) * (T - Topt) - (Topt - Tmax) * (Topt + Tmin - 2.0 * T))) * (pH - pHmin) * (pH - pHmax) / ((pH - pHmin) * (pH - pHmax) - (pH - pHopt)^2) * (aw - 1) * (aw - awmin)^2 / ((awopt - awmin) * ((awopt - awmin) * (aw - awopt) - (awopt - 1) * (awopt + awmin - 2.0 * aw)))))



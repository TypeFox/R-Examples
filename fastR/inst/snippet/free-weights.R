( 4.96 - 5.0 ) / (0.05 / sqrt(10) ) -> z; z    # test statistic
2 * (pnorm(z))                                 # 2-sided p-value
2 * (1 - pnorm(abs(z)))                        # 2-sided p-value again

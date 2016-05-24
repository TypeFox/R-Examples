-0.0946 + c(-1,1) * 0.01713 * qt(0.975,df=14)  # CI by hand
confint(treb.model,"projectileWt")            # CI using confint()

\dontrun{#REX
library(psd)

##
## Confidence intervals from taper numbers
##

sp <- spectral_properties(as.tapers(1:50), p=0.95, db.ci=TRUE)

# standard errors as a function of tapers
par(las=1)
plot(stderr.chi.upper ~ taper, sp, type="s",
       ylim=c(-10,20), yaxs="i", xaxs="i",
       xlab=expression("number of tapers ("* nu/2 *")"), ylab="dB",
       main="Spectral uncertainties")
mtext("(additive factor)", line=.3)
lines(stderr.chi.lower ~ taper, sp, type="s")
lines(stderr.chi.median ~ taper, sp, type="s", lwd=2)
lines(stderr.chi.approx ~ taper, sp, type="s", col="red",lwd=2)
# indicate K needed to reach 3 dB wide confidence interval (p=.95)
abline(v=33, lty=3)
legend("topright",
        c(expression("Based on "* chi^2 *"(p,"*nu*") and (1-p,"*nu*")"),
          expression(""* chi^2 *"(p=0.5,"*nu*")"),
          "approximation"),
       lwd=c(1,3,3),
       col=c("black","black","red"),
       bg="grey98")

}#REX

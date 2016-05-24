cohd2delta <-
function (d) {#U1 in Table 2.2.1 in Cohen's Statistical Power Analysis, p. 22
 return (((2*pnorm(d/2))-1)/pnorm(d/2)) }


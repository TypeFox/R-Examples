
# sine animation example Initialize the gnuplot handle
h1 <- Gpinit()
print("*** example of gnuplot control through C ***\n")
# loop from 0.1 to 1.0
for (phase in (1:100)/10) {
    # reset and plot the sine function shifted by the loop variable
    Gpresetplot(h1)
    Gpcmd(h1, "plot sin(x+%g)", phase)
}
# loop from 0.1 to 1.0
for (phase in (100:1)/10) {
    # reset and plot the sine function shifted by the loop variable
    Gpresetplot(h1)
    Gpcmd(h1, "plot sin(x+%g)", phase)
}
# close gnuplot handle
h1 <- Gpclose(h1)
 

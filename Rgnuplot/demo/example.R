
# example of plotting slopes, equations, styles, lists of points and multiple output windows gnuplot_i, by Nicolas Devillard Initialize variables
NPOINTS <- 50
SLEEP_LGTH <- 2
x <- vector("numeric", NPOINTS)
y <- vector("numeric", NPOINTS)
# Initialize the gnuplot handle
print("*** example of gnuplot control through C ***\n")
h1 <- Gpinit()
h2 <- Gpinit()

# Slopes

Gpsetstyle(h1, "lines")

print("*** plotting slopes\n")
print("y = x\n")
GpplotSlope(h1, 1, 0, "unity slope")
Sys.sleep(SLEEP_LGTH)

print("y = 2*x\n")
GpplotSlope(h1, 2, 0, "y=2x")
Sys.sleep(SLEEP_LGTH)

print("y = -x\n")
GpplotSlope(h1, -1, 0, "y=-x")
Sys.sleep(SLEEP_LGTH)

# Equations

Gpresetplot(h1)
print("\n\n*** various equations\n")
print("y = sin(x)\n")
GpplotEquation(h1, "sin(x)", "sine")
Sys.sleep(SLEEP_LGTH)

print("y = log(x)\n")
GpplotEquation(h1, "log(x)", "logarithm")
Sys.sleep(SLEEP_LGTH)

print("y = sin(x)*cos(2*x)\n")
GpplotEquation(h1, "sin(x)*cos(2*x)", "sine product")
Sys.sleep(SLEEP_LGTH)

# Styles

Gpresetplot(h1)
print("\n\n")
print("*** showing styles\n")

print("sine in points\n")
Gpsetstyle(h1, "points")
GpplotEquation(h1, "sin(x)", "sine")
Sys.sleep(SLEEP_LGTH)

print("sine in impulses\n")
Gpsetstyle(h1, "impulses")
GpplotEquation(h1, "sin(x)", "sine")
Sys.sleep(SLEEP_LGTH)

print("sine in steps\n")
Gpsetstyle(h1, "steps")
GpplotEquation(h1, "sin(x)", "sine")
Sys.sleep(SLEEP_LGTH)

# User defined 1d and 2d point sets

Gpresetplot(h1)
Gpsetstyle(h1, "impulses")
print("\n\n")
print("*** user-defined lists of doubles\n")
for (i in 1:NPOINTS) {
    x[i] = i * i
}
GpplotX(h1, x, NPOINTS, "user-defined doubles")
Sys.sleep(SLEEP_LGTH)

print("*** user-defined lists of points\n")
for (i in 1:NPOINTS) {
    x[i] = i
    y[i] = i * i
}
Gpresetplot(h1)
Gpsetstyle(h1, "points")
GpplotXY(h1, x, y, NPOINTS, "user-defined points")
Sys.sleep(SLEEP_LGTH)


# Multiple output screens

print("\n\n")
print("*** multiple output windows\n")
Gpresetplot(h1)
Gpsetstyle(h1, "lines")
h2 = Gpinit()
Gpsetstyle(h2, "lines")
h3 = Gpinit()
Gpsetstyle(h3, "lines")
h4 = Gpinit()
Gpsetstyle(h4, "lines")

print("window 1: sin(x)\n")
GpplotEquation(h1, "sin(x)", "sin(x)")
Sys.sleep(SLEEP_LGTH)
print("window 2: x*sin(x)\n")
GpplotEquation(h2, "x*sin(x)", "x*sin(x)")
Sys.sleep(SLEEP_LGTH)
print("window 3: log(x)/x\n")
GpplotEquation(h3, "log(x)/x", "log(x)/x")
Sys.sleep(SLEEP_LGTH)
print("window 4: sin(x)/x\n")
GpplotEquation(h4, "sin(x)/x", "sin(x)/x")
Sys.sleep(SLEEP_LGTH * 5)

# close gnuplot handles
h1 <- Gpclose(h1)
h2 <- Gpclose(h2)
h3 <- Gpclose(h3)
h4 <- Gpclose(h4)

 

# Ajuste de parametros graficos
require(grDevices)
require(graphics)
opar <- par(ask = dev.interactive(orNone = TRUE))

# Creacion de un objeto loca.p
o <- new("loca.p", x = c(-1, 1, 0), y = c(0, 0, 1))

# Imprimiendo un resumen del objeto
print(o)

# Grafica de los puntos de demanda
plot(o)

# Evaluacion de la funcion en el punto (0, 0.5)
zsum(o, x=0, y=0.5)

# Grafica de las curvas de nivel de la funcion objetivo
contour.loca.p(o)

# Grafica 3D de la funcion objetivo
persp.loca.p(o)

# Grafica 3D
persp.loca.p(o, col=cm.colors(10000), border=FALSE, shade=TRUE, theta=50, phi=5, ltheta=135)

# Otra grafica 3D
persp.loca.p(o, col=cm.colors(10000), border=FALSE, shade=TRUE, theta=50, phi=5, ltheta=135, lphi=90)

# Busqueda del minimo
zsummin(o)

# Nuevo objeto aleatorio loca.p con 10 puntos de demanda
p <- rloca.p(10)

# Busqueda del minimo 
sol <- zsummin(p)

# El minimo
sol

# Evaluacion de la funcion en el minimo
zsum(p, sol[1], sol[2])

# Calculo del tiempo de ejecucion
system.time(zsummin(p))

# Restauracion de los parametros graficos
par(opar)

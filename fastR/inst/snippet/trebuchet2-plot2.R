trebuchet2.bandplot <- xyplot(distance~projectileWt,
    data=trebuchet2, ylim=c(2.5,10.5),
    panel=panel.lmbands,
    conf.lty=1,
    pred.lty=1
    )

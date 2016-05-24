predict(treb.model,newdata=data.frame(projectileWt=44))
predict(treb.model,newdata=data.frame(projectileWt=44),
    interval='confidence')
predict(treb.model,newdata=data.frame(projectileWt=44),
    interval='prediction')

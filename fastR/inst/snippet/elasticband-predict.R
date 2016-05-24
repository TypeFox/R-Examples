predict(eband.model,newdata=data.frame(stretch=30))
predict(eband.model,newdata=data.frame(stretch=30),interval='confidence')
predict(eband.model,newdata=data.frame(stretch=30),interval='prediction')

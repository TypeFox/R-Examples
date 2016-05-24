
Gprun("k=0\nbind Home \"k=k+1; replot\"\nplot sin(x+k)", TRUE)

Gprun("bind \"p\" \"set terminal png;set output '00000.png';replot;set terminal wxt\"\nplot sin(x)", TRUE)

Gprun("bind \"p\" \"n=strftime('%d.%b.%Y-%H.%M.%.3S.png',time(0.0));set terminal png;set output n;replot;set terminal wxt\"\nplot sin(x)", TRUE)
 

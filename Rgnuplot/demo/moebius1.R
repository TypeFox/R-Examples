# world map on a moebius strip
Gprun("bind \"p\" \"n=strftime('moebiusA%d.%b.%Y-%H.%M.%.3S.png',time(0.0));set terminal png;set output n;replot;set terminal wxt\"\nset angles degrees\nset parametric\nset hidden\nf(u,v) = cos(u) + v*cos(u/2.0) * cos(u)\ng(u,v) = sin(u)+v*cos(u/2.0)*sin(u)\nh(u,v) = .5*v*sin(u/2.0)\nunset key;unset tics;unset border;unset colorbox\n#set view equal xyz\nset hidden3d\nk=200.0\nsplot \"NOAACoastline.dat\" u (f(-$1,$2/k)):(g(-$1,$2/k)):(h(-$1,$2/k))  w l lc rgb \"black\" notit\n", 
    TRUE) 

a = a+1
set xrange [0:50]
set yrange [0:50]
unset key
plot "matrix-0.95,L=50.txt" pt 7 ps 2.5
pause 0
if(a<=230000) reread
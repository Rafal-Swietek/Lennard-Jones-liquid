a = a+1
set xrange [0:40]
set yrange [0:40]
unset key
plot "matrix.txt" pt 7 ps 1.5
pause 2
if(a<=30000) reread
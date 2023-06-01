#Script para gnuplot que genera la trayectoria del cohete 
set terminal gif animate delay 20.2
set output 'Cohete.gif'

set xlabel 'X '
set ylabel 'Y '

set xrange [-0.2:1.2]
set yrange [-0.2:1.2]

set style fill solid 1.0 border -1

set object circle at -0.2,-0.2 size scr 0.1 fc rgb "blue"


do for [ii=0:49] {
  plot 'Tierra-luna.txt' u 1:2 every ::ii::ii w p pt 3 ps 2 lc rgb "black" t 'Cohete','' u 3:4 every ::ii::ii  with circles lc rgb "gray0"  t 'Luna'
}

set output

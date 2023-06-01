# Establece el título de la gráfica
set title "Grafica de H' en funcion del tiempo"

# Establece el nombre del eje x y el rango de valores
set xlabel "t"

# Establece el nombre del eje y y el rango de valores
set ylabel "H'"

# Establece el tipo de línea y el título de la leyenda
plot "Hamiltoniano.txt" using 1:2 with lines title "H'"

                                        



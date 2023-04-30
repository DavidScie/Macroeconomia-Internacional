#Modelo Ramsey con shock a la MIT
#Parámetros
beta = 0.98 #factor de descuento intertemporal
delta = 0.2 #depreciación del capital

alpha = 1/3 
r = (1/beta) - 1 #tasa de interés real
ctiny = 6.3829e-4 
crit = 1          # Criterio de convergencia          
tol = 0.0001 #tolerancia para determinar la convergencia        

#Stock de capital estacionario y malla
k0 = ((1/alpha)*(r + delta))^(1/(alpha - 1))    #Es el estado estacionario del Capital
kpoints = 1000                                  # Puntos del krid
dev = 0.9                                       # Desviación en torno al Estado Estacionario  
kmin = k0 * (1 - dev)                          
kmax = k0 * (1 + dev) 
dk = (kmax - kmin) / (kpoints - 1)              # Calcula la tasa de incremento
k = kmin:dk:kmax                                # Distintos valores para el capital   
#parte distinta del codigo de matlab
tv=zeros(size(k))               #vector de ceros de tamaño 1000*1
indice_i = zeros(size(k)) 
indice_i = floor.(Int, indice_i)  
k = collect(k)              # vector columna de los krid de puntos de k con tasa de crecimiento positiva

#Adivinanza funcion de valor
v = zeros(size(k)) #es un vector {Float64} igual que tv, debe ser asi para que se puedan restar en el error_i 

#Posibilidades de consumo
c = (((k.^alpha) + (1 - delta)*k)*ones(size(k'))) - ones(size(k))*k' 
# es una matriz de 1000*1000, el cual genera todas las posibilidades de consumo.
c[c.==0].=ctiny;                     # Saca los consumos negativos
c[c.<=0].=ctiny/1000000000000000000000000   # Saca los consumos pequeños, por la función logaritmica
u = log.(c[:,:])                            
#es una matriz que genera la Funcion de utilidad con todas las posibilidades

#aplicación del CMT
j = 1 #para ver en cuantas iteraciones converge
while crit > tol
    for i=1:kpoints
        tv[i],indice_i[i] = findmax(u[i,:] + beta*v) #findmax entrega el maximo valor y su indice (posición)
        end    
    error_i=abs.(tv - v)
    crit = maximum(error_i)                
    copyto!(v,tv)
   global v 
   global crit 
   global j += 1    
end 
j

# Solución Final Numérica

copt  = k.^alpha + (1 - delta)*k - k[indice_i]   #calculo numerico
kk_0  = k[indice_i] 
b1    =   (log.(k[:,:]).- log(k0))\(log.(kk_0[:,:]).- log(k0)) 
b11   =   b1[1,1] 
typeof(b11)


kkk_0 = exp.(b11*(log.(k[:,:]).- log(k0)).+log(k0)) 

copte = k.^alpha + (1 - delta)*k - kkk_0

using Plots
plot(k,k, label = "Línea de 45 grados") #mismo consumo y capital en todos los periodos
plot!(k, kk_0, label = "Capital óptimo")
plot!(k, copt, label = "Consumo óptimo")
xlabel!("Capital")
ylabel!("Consumo")
title!("Consumo óptimo")
nombre_vars = ["Linea de 45 grados" "Capital óptimo" "Consumo óptimo"]
plot1 = plot(k, [k kk_0 copt], xlabel="Capital", ylabel="Consumo", label=nombre_vars, title="Consumo óptimo", linewidth = 2, 
line=[:solid :dot :solid], color=[:red :blue :green])
savefig(plot1, "cons_opt.png") 


## Simulación Shock MIT
N        = 20;
kk       = zeros(N,1)
AA       = zeros(N,1)
kk[1]  = 1.5*k0
log(k0)
kkk      = zeros(N,1)
#kkk[1,:] .= exp(b11*(log(kk[1,;]). - log(k0)). + log(k0))
kkk[1]   = exp.(b11*(log(kk[1]) - log(k0)) + log(k0))
cosim    = zeros(N,1)
#cosim(1) = (kk(1)^alpha) + (1 - delta)*kk(1) - kkk(1);
cosim[1] = kk[1].^alpha + (1 - delta)*kk[1] - kkk[1]
AA[1]    = (cosim[1] + kkk[1] - (1 - delta)*kk[1])/(kk[1]^alpha)
   
for j=2:N
    
        kkk[j]   = exp.(b11*(log(kkk[j-1]) - log(k0)) + log(k0));
        cosim[j] = kkk[j-1]^alpha + (1 - delta)*kkk[j-1] - kkk[j];
        
 end
#gráfico
 t=1:1:N
nomb_vars1 = ["Consumo" "Capital"]
plt1= plot(t, [cosim kkk], ylabel ="Unidades de consumo y capital", xlabel="Tiempo", label=nomb_vars1, title="IRF", linewidth = 2)
savefig(plt1, "IRF_MIT.png") 

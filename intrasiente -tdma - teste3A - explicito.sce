clear
funcprot(0)


function malha_1d = malha_simetrica(L,div_malha)
    malha_1d = zeros(1, div_malha);
    malha_1d (1) = 0
    for i=2:div_malha
        malha_1d(i)= L*(i-1)/(div_malha-1)
            end
endfunction

function criar_malha = malha_1d(a, l1)       //função que cria a malha 
criar_grid = zeros(1, l1+1);
    alp = 0.2d0;                             // quanto maior o alp maior refino da malha nas extremidades e menor no centro
    for i = 1:l1
      criar_grid (i) = a * (((-0.5d0 * tanh(alp * ((2.d0 * (i - 1.d0) / (l1 - 1.d0)) - 1.d0))) / (tanh(-alp))) + 0.5d0);
//  criar_grid (i) = a * (((-0.5d0 * tanh(alp * ((2.d0 * (i - 0) / (l1 - 0)) - 1.d0))) / (tanh(-alp))) + 0.5d0);
    end
    for i=1:l1
    criar_malha(i) = criar_grid(i)
    end
endfunction

//----DADOS DE MONTAGEM DA MALHA-------------------//
L = 0.02                                           // comprimento [m]
div_malha = 10+1                                   // numero de elementos [adm]
malha = malha_simetrica (L,div_malha)              // chamada da função que monta a malha
//malha = malha_1d (L,div_malha)                   // chamada da função malha assimetrica
//----FIM MONTAGEM DA MALHA------------------------//

k = 10              // [W/m.K]
p = 1               //massa especifica [Kg/m³]
Cp = 1e07          // condutibilidade térmica [J/Kg.K]
tempo_final = 2
t_div = 0.2           //segundos [s]


x(1) = malha(1) 
for i=2:div_malha
    x(i) = malha(i)-malha(i-1)  //matriz da distancia entre pontos (distancia L)
end

//--------MATRIZ PERFIL INICIAL DE TEMPERATURA-----//
for i=1:div_malha-1                                          //
matriz_Temp_Tp0(i)= 200                            // utilizar graus Kelvin
    end                                            //
                                                   //
for i=div_malha:div_malha                                  //
matriz_Temp_Tp0(i)= 0                              // utilizar graus Kelvin
end                                                //
//--------FIM  PERFIL INICIAL DE TEMPERATURA-------//

disp(matriz_Temp_Tp0)
array_TEMP = zeros(1,div_malha)

TEMPO = t_div


while TEMPO <= tempo_final

erro=10000
while erro > 1e-10

erro_a = array_TEMP(1)

//--------------------Array_Temp
for i=1:1
    ae = k / x(i+1)
ap0 = p*Cp* x(i+1) /t_div
array_TEMP(i) =( ae* matriz_Temp_Tp0(i+1)+ (ap0-ae)* matriz_Temp_Tp0(i))/ap0
end

for i=2:div_malha-1
    ae = k / x(i+1)
    aw = k / x(i)
ap0 = p*Cp* x(i) /t_div
array_TEMP(i) = (ae* matriz_Temp_Tp0(i+1)+ aw * matriz_Temp_Tp0(i-1)+(ap0-ae-aw)* matriz_Temp_Tp0(i))/ap0
end

for i=div_malha:div_malha
    aw = k / x(i-1)
ap0 = p*Cp* x(i) /t_div
array_TEMP(i) = 0
end
disp(array_TEMP)
//--------------------Array_Temp end

erro_b = array_TEMP(1)
erro = abs(erro_a - erro_b)/erro_b

end  // while (y) L64

disp(array_TEMP)
matriz_Temp_Tp0 = array_TEMP
TEMPO = TEMPO+t_div
printf ("tempo = %.2f \n", TEMPO)
end     // end do loop de tempo L61




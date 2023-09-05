clear
funcprot(0)

//--------------------------------------------------------
function malha_1d = malha_simetrica(L,div_malha)       //função que cria a malha simétrica (habilitar L48)
    malha_1d = zeros(1, div_malha);
    malha_1d (1) = 0
    for i=2:div_malha
        malha_1d(i)= L*(i-1)/(div_malha-1)
            end
endfunction
//--------------------------------------------------------
function criar_malha = malha_1d(a, l1)       //função que cria a malha assimetrica. ()habilitar L49) 
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
//-------------------------------------------------------
function condutividade_termica= calc_K(temperatura)
    condutividade_termica = 0.010*temperatura           // [W/mK] equação da variação da condutividade termica
    //condutividade_termica = 10
endfunction
//---------------------------------------------------------
function rho_cp = calc_p_Cp(Pressao, Temperatura)       // capacidade calorifica volumetrica  (ρ*Cp) [J/m³.K]
    W = 28.966                                          // [g/mol]
    W = W/1000                                          // conversão para Kg/mol
    R = 8.3144598                                       // [J/K.mol]   constante universal dos gases perfeitos
    massa_gas = (W/R)*(Pressao/Temperatura);            // [Kg/m³]   RESULTADO 
       
    E = 0.95                                            // proporção gás/sólido  (valor entre 0 e 1))
    Cp_gas =  1298                                      // unidade [J/Kg .K]   calor especifico do gas
    massa_solid = 8900                                  // unidade [Kg/m³]     massa especifica do sólido
    Cp_solido = 661                                     // unidade [J/Kg .K]   calor especifico do sólido
     rho_cp = E*massa_gas*Cp_gas+(1-E)*massa_solid*Cp_solido;
    // rho_cp = 1e7
    endfunction


//----DADOS DE MONTAGEM DA MALHA-------------------//
L = 0.02                                           // comprimento [m]
div_malha = 10+1                                   // numero de elementos [adm]
malha = malha_simetrica (L,div_malha)              // chamada da função malha simétrica
//malha = malha_1d (L,div_malha)                   // chamada da função malha assimetrica
//----FIM MONTAGEM DA MALHA------------------------//


tempo_final = 1    //segundos [s] -> tempo total da simulação
t_div = 0.001           //segundos [s] -> passo de tempo da simulação


//--------MATRIZ DE PERDA DE PRESSÃO ----------------------//
Pressao_inicial = 101325+8963.18                           // Pascal [N/m²]
Pressao_final = 101325                                     // Pascal [N/m²]
//perda_pressao = malha/L * (Pressao_inicial - Pressao_final)//
//----FIM MATRIZ DE PERDA DE PRESSÃO ----------------------//


x(1) = malha(1) 
for i=2:div_malha
    x(i) = malha(i)-malha(i-1)  //matriz da distancia entre pontos (distancia dx)
end

//--------MATRIZ PERFIL INICIAL DE TEMPERATURA-----//
for i=1:1                                           //
matriz_Temp_Tp0(i)= 2564                            // utilizar graus Kelvin
    end                                            //
                                                   //
for i=2:div_malha                                  //
matriz_Temp_Tp0(i)= 293                           // utilizar graus Kelvin 
end                                                //
//--------FIM  PERFIL INICIAL DE TEMPERATURA-------//



for i=1:div_malha
    array_TEMP(i) = 273
end



TEMPO = t_div


while TEMPO <= tempo_final

erro=10000
while erro > 1e-10

erro_a = array_TEMP(1)

for i=1:1
    k = calc_K(array_TEMP(i))
    ae = k / x(i+1)
    pressao_elemento = Pressao_inicial - malha(i)/L * (Pressao_inicial - Pressao_final)
    rho_cp = calc_p_Cp(pressao_elemento, array_TEMP(i))
    ap0 = rho_cp* x(i+1) /t_div
array_TEMP(i) = matriz_Temp_Tp0(1) //( ae* array_TEMP(i+1)+ (ap0)* matriz_Temp_Tp0(i))/(ap0+ae)
end

for i=2:div_malha-1
    k = calc_K(array_TEMP(i+1))
    ae = k / x(i+1)
    k = calc_K(array_TEMP(i))
    aw = k / x(i)
    pressao_elemento =Pressao_inicial - malha(i)/L * (Pressao_inicial - Pressao_final)
    rho_cp = calc_p_Cp(pressao_elemento, array_TEMP(i))
ap0 = rho_cp* x(i) /t_div
array_TEMP(i) = (ae* array_TEMP(i+1)+ aw * array_TEMP(i-1)+(ap0)* matriz_Temp_Tp0(i))/(ap0+ae+aw)
end

for i=div_malha:div_malha
    k = calc_K(array_TEMP(i-1))
    aw = k / x(i-1)
    pressao_elemento = Pressao_inicial - malha(i)/L * (Pressao_inicial - Pressao_final)
    rho_cp = calc_p_Cp(pressao_elemento, array_TEMP(i))
ap0 = rho_cp* x(i) /t_div
array_TEMP(i) = ( aw* array_TEMP(i-1)+ (ap0)* matriz_Temp_Tp0(i))/(ap0+aw) // 2570
end
//disp(array_TEMP)  // habilite a linha para exibir as iteração passo a passo;

erro_b = array_TEMP(1)
erro = abs(erro_a - erro_b)/erro_b

end  // while (y) L64

disp(array_TEMP)
matriz_Temp_Tp0 = array_TEMP
TEMPO = TEMPO+t_div
printf ("tempo = %.2f \n", TEMPO)
end     // end do loop de tempo L61

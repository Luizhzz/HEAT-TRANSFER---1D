clear
funcprot(0)

function x = TDMA(a, b, c, d)
    n = size(a, 1);
    cl = zeros(n, 1);
    dl = zeros(n, 1);
    x = zeros(n, 1);

    cl(1) = c(1) / b(1);
    for i = 2:n-1
        cl(i) = c(i) / (b(i) - a(i) * cl(i-1));
    end

    dl(1) = d(1) / b(1);
    for i = 2:n
        dl(i) = (d(i) - a(i) * dl(i-1)) / (b(i) - a(i) * cl(i-1));
    end

    x(n) = dl(n);
    for i = n-1:-1:1
        x(i) = dl(i) - cl(i) * x(i+1);
    end
endfunction


div_malha = 5 
k = 1000
x = 0.1

TEMP_entrada = 500
TEMP_saida = 100

fonte = 0
ae = k/x
aw = k/x
ap = ae + aw + fonte
ap0 = 0


//--------------------- Ap1
array_ap (1) = ae + aw*2
for i=2:div_malha-1
       array_ap (i) = ap
end
array_ap (div_malha)= ae*2 + aw
//--------------------- Ap1 end

//----------------------Aw
for i=1:div_malha-1
       array_aw (i) = -aw
end
array_aw(div_malha)= 0
//---------------------Aw end

//---------------------Ae
array_ae (1) = 0
for i=2:div_malha
       array_ae (i) = -ae
end
//---------------------Ae end

//--------------------Ap0
array_ap0 (1) = aw*2*TEMP_entrada
for i=2:div_malha-1
       array_ap0 (i) = ap0
end
array_ap0 (div_malha) = ae*2*TEMP_saida
//--------------------Ap0 end

array_TEMP = TDMA (array_ae, array_ap, array_aw, array_ap0)

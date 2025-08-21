function dw=lorenz_HGO(w1,w2,w3)

%%Parámetros


a=7.5; %%default a=10; %tambien funciona con 5
b=8/3;
c=28;

 w=[w1;w2;w3];

S=[-a a 0
    c -1 -w1
    w2 0 -b];

dw=S*w;
end

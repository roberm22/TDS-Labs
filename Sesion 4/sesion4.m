
% Sesión 4: Filtros digitales, compensación y fase lineal

% Soluciones

% Autor
% Roberto Martin

% Filtros digitales
%% Programación de filtros digitales

% Apartado 1.2
% Comprobación de la función sso.m

clear

% Datos
B=[1 0 -0.25]; % esto son los coeficientes del numerador
A=[1 -0.9*sqrt(2) 0.81]; % esto son los coeficientes del denominador
N=50;

% Impulso unidad
n=0:N-1;        % Eje temporal n
d=double(n==0); % Impulso unidad

% Respuesta impulsiva con sso(. . .)
ci=zeros(2,1);  % Reposo inicial
[hsso cf]=sso(B,A,d,ci); % No funcionará si no completa sso

% Respuesta impulsiva con filter para comparar con sso
hfilter=filter(B,A,d);

syms z;
H=(1-0.25*(z^-2))/(1-0.9*sqrt(2)*(z^-1)+0.81*(z^-2));

% Gráficas
stem(n,hsso,'g')
  grid on;
  hold on;
  xlabel('n');
stem(hfilter,'r') % añada código
grid, xlabel('n')
legend('h_{sso}[n]', 'h_{filter}[n]','Location','northeast')
hold off

% Comparación numérica
disp(' ')
disp(['max(abs(hsso[n]-hfilter[n])) = ' num2str(max(abs(hsso-hfilter)))])


%% Apartado 1.3
% Filtrado de una señal

clear

% Datos
load tds
Tv= 12.5e-3; % Añada código
Fs= 8e3; % Añada código
Ts=1/Fs;
Fc= 1e3; % Añada código
[B A]=ellip(2,0.5,20,Fc/(Fs/2)); % Añada código

% Filtrado localizado
Lv=fix(Tv/Ts);
Ld=Lv;
ci=zeros(2,1);
tdsf=[];
xl=slocal(tds,Lv,Ld,0);
while 1
    xl=slocal(tds,Lv,Ld,1);
    if isempty(xl)
        break
    end
    [y,cf]=sso(B,A,xl,ci);
    tdsf = [tdsf; y]; % Añada código
    ci=cf;
end

% Gráficas
n=7250:7330;
plot(n,tds(n),'g',n,tdsf(n),'r')
grid, xlabel('n')
legend('Entrada del filtro', 'Salida del filtro')

% Audición
%wavplay([int16(tds') zeros(1,Fs) int16(tdsf')],Fs)
p=audioplayer([int16(tds') zeros(1,Fs) int16(tdsf')],Fs);
play(p);

% COMENTARIOS
% Se aprecia la señal filtrada retrasada 1-2 muestras y ligeramente
% atenuada (como era de esperar). La señal filtrada se aprecia en la
% audición más de "baja frecuencia", pero perfectamente inteligible.

%% Ejercicio 3
% Problema 5.53 de:
% A. V. Oppenheim & R. W. Schafer: “Discrete-Time Signal Processing”, (3rd
% Edition), Ed. Prentice Hall 
%(se suministran los enunciados de ambos problemas en español)
clear
%(se suministran los enunciados de ambos problemas en español)
clear
% Datos
z1=0.9*exp(1i*0.6*pi); 
z2=1.25*exp(1i*0.8*pi);

% Diagramas polos-ceros
Z1=[z1; z1'; z2; z2'];
Z2=[1/z1'; 1/z1; z2; z2']; % reflejar sólo z1
Z3=[z1; z1'; 1/z2'; 1/z2]; % reflejar sólo z2
Z4=[1/z1'; 1/z1; 1/z2'; 1/z2]; % reflejar ambos z1 y z2

P=[0;0;0;0];  %vector de polo 4-druple en el origen
subplot(2,2,1), zplane(Z1,P) % Alternativamente zplane(Z1)
subplot(2,2,2), zplane(Z2,P) % Alternativamente zplane(Z2)
subplot(2,2,3), zplane(Z3,P) % Alternativamente zplane(Z3)
text(-2,-1,'Sistema de fase mínima')
subplot(2,2,4), zplane(Z4,P) % Alternativamente zplane(Z4)
subplot
pause

% Resp. impulsivas normalizadas a la ganancia en continua del 1er sistema.
h1=zp2tf(Z1,P,1); G1=sum(h1);
h2=zp2tf(Z2,P,1); G2=sum(h2);h2=h2*G1/G2;
h3=zp2tf(Z3,P,1); G3=sum(h3);h3=h3*G1/G3;
h4=zp2tf(Z4,P,1); G4=sum(h4);h4=h4*G1/G4;

n=0:length(h1)-1;
plot(n,h1,'-or',n,h2,'-og',n,h3,'-ob',n,h4,'-ok')
title('Respuestas impulsivas de los cuatro sistemas')
xlabel('n'), grid
text(0.5,0.75,'El sistema de fase mínima en azul')
pause

% Energías parciales
E1=cumsum(h1.^2);
E2=cumsum(h2.^2);
E3=cumsum(h3.^2);
E4=cumsum(h4.^2);
plot(n,E1,'-or',n,E2','-og',n,E3,'-ob',n,E4,'-ok') %añadir código
title('Energías parciales de los cuatro sistemas')
xlabel('n'), grid
text(0.5,2,'El sistema de fase mínima en azul')
pause

% Respuestas en frecuencia en magnitud y fase
[H1 F]=freqz(h1,1,1000); MH1=20*log10(abs(H1)); FH1=unwrap(angle(H1)); Gd1=grpdelay(h1,1,1000);
[H2 F]=freqz(h2,1,1000); MH2=20*log10(abs(H2)); FH2=unwrap(angle(H2)); Gd2=grpdelay(h2,1,1000);
[H3 F]=freqz(h3,1,1000); MH3=20*log10(abs(H3)); FH3=unwrap(angle(H3)); Gd3=grpdelay(h3,1,1000);
[H4 F]=freqz(h4,1,1000); MH4=20*log10(abs(H4)); FH4=unwrap(angle(H4)); Gd4=grpdelay(h4,1,1000);
plot(F/pi,MH1,'-r',F/pi,MH2,'-g',F/pi,MH3,'-b',F/pi,MH4,'-k')% Añadir código
title('Respuestas en magnitud de los cuatro sistemas')
xlabel('\omega (x \pi rad/muestra)'), ylabel('dB'), grid
text(0.5,-12,'El sistema de fase mínima en azul')
pause


plot(F/pi,FH1,'-r',F/pi,FH2,'-g',F/pi,FH3,'-b',F/pi,FH4,'-k')
title('Respuestas en fase de los cuatro sistemas')
xlabel('\omega (x \pi rad/muestra)'), 
ylabel('\angleH('') (x \pi rad)'), grid %añadir código
text(0.5,-13,'El sistema de fase mínima en azul')
pause

plot(F/pi,Gd1,'-r',F/pi,Gd2,'-g',F/pi,Gd3,'-b',F/pi,Gd4,'-k')
title('Retardo de grupo de los cuatro sistemas')
text(0.5,-7,'El sistema de fase míonima en azul')
xlabel('\omega (x \pi rad/muestra)'), ylabel('\tau_g(n) (muestras)'), grid

% EL SISTEMA DE FASE MÍNIMA EN AZUL


%% Ejercicio 4
% Problema 5.38 de:
% A. V. Oppenheim & R. W. Schafer: “Discrete-Time Signal Processing”, (3rd
% Edition), Ed. Prentice Hall 
%(se suministran los enunciados de ambos problemas en español)
clear

% Sistema H(z)
z1= -1/3 ; z2= -3; z3=3; % Añadir código
p1=1/3;
Z = [z1; z2; z3];
P = [p1; 0; 0];
[B,A]= zp2tf(Z,P,1); % Añadir código

pause(1)

zap1=3;
pap1=1/3;
pap2=-1/3;
Zap=[zap1;0];
Pap=[pap1;pap2];

zmin1=1/3;
zmin2=1/3;
zmin3=-1/3;
pmin1=1/3;
Zmin=[zmin1;zmin2;zmin3];
Pmin=[pmin1;0;0];

% Sistemas Hmin y Hap
[Bap, Aap]=zp2tf(Zap,Pap,1);
[Bmin, Amin]=zp2tf(Zmin,Pmin,1);
% Ahora se cumple que H(z)= Hmin(z)*Hap(z), pero |Hap(w)|= cte != 1. Para
% normalizar la constante a 1, escribimos H(z)= ((1/k)*Hmin(z))*(k*Hap(z)) 
% y forzamos a que k*Hap(z=1)=1, con lo cual k= sum(Aap)/sum(Bap)
%k=Bap(end);
k= sum(Aap)/sum(Bap);
% a) Corrijo Hap(z)
Bap=Bap*k; 
% b) corrijo la ganancia del sistema de fase mínima con el inverso
Bmin=Bmin/k;
% Diagramas polos ceros de la factorización H(z)=Hmin(z)*Hap(z)
subplot(2,1,1), zplane(Bmin,Amin)
text(0.5,0.75,'Sistema de fase mínima')

subplot(2,1,2), zplane(Bap,Aap)% Añadir código
text(0.5,0.75,'Sistema paso todo')
subplot

%Respuesta en frecuencia
[H1 F]=freqz(Bap,Aap,1000); MH1=20*log10(abs(H1)); 
FH1=unwrap(angle(H1)); Gd1=grpdelay(Bap,Aap,1000);

subplot (3,2,1)
plot(F/pi,MH1)
title('Paso-todo')
ylabel('|H(\omega)| dB'), grid
axis([0 1 -10 10])
subplot(3,2,3)
plot(F/pi,FH1/pi)
ylabel('\angle H(\omega)(x\pi)'), grid
subplot(3,2,5)
plot(F/pi,Gd1)
xlabel('\omega (x \pi rad/muestra)'), ylabel('\tau_g(\omega)'), grid
pause (1)

%[H2 F]=freqz(h2, 1 ,1000); MH2=20*log10(abs(H2));% Añadir código
%FH2=unwrap(angle(H2)); Gd2=grpdelay(h2, 1 ,1000);% Añadir código

[H2 F]=freqz(Bmin, Amin ,1000); MH2=20*log10(abs(H2));% Añadir código
FH2=unwrap(angle(H2)); Gd2=grpdelay(Bmin, Amin ,1000);% Añadir código

subplot (3,2,2)
plot(F/pi,MH2)
ylabel('|H(\omega)| dB'), grid
title('Fase mínima')
axis([0 1 5 25])
subplot(3,2,4)
plot(F/pi,FH2/pi)
ylabel('\angle H(\omega)(x\pi)'), grid
subplot(3,2,6)
plot(F/pi,Gd2)
xlabel('\omega (x \pi rad/muestra)'), ylabel('\tau_g(\omega)'), grid
pause

% Apartado d) EL SISTEMA Hmin(z) ES FIR PERO NO DE FASE LINEAL GENERAL
% Si se puede hacer la factorización H(z)=Hlin(z)*Hap2(z)

z=tf(z);
Hlin=(1+3*z^-1)*(1+(1/3)*z^-1);
Hap2=(1-3*z^-1)/(1-(1/3)*z^-1);
[plin,clin] = pzmap(Hlin);
%Blin=poly([z2 z3]); Alin=1;
[Blin,Alin]=zp2tf(clin,plin,1);% Añadir código
%Bap2=poly(z1); Aap2=poly(p1);
[Bap2,Aap2]=zp2tf([z1],[p1],1);
%Ajusto las ganancias para que Hap tenga ganancia constante unidad
k=Bap2(end);
Bap2=Bap2/k;Blin=Blin*k;
subplot(2,1,1), zplane(Blin,Alin)
text(0.5,0.75,'Sistema de fase lineal')

subplot(2,1,2), zplane(Bap2,Aap2)
text(0.5,0.75,'Sistema paso todo')

subplot
%
% EJERCICIOS ADICIONALES
%% Programación de un flujograma
% Ejercicio 1

clear

% Datos
B=[0 1/4 0 1];
A=[1 0 -1/4];
N=20;

% Impulso unidad
n=0:N-1;        % Eje temporal n
d=double(n==0); % Impulso unidad

% Respuesta impulsiva con flujog(. . .) Respuesta impulsiva con filter
% Gráficas

h_filter=filter(B,A,n,d);
[y,mf]=flujog(n,d);
h=y;
 figure('Name','Flujograma','NumberTitle','off');
    stem(h,'Color','g');
        grid on;
        hold on;
        xlabel('n');
    stem(h_filter,'Color','r');
        legend('h_f_l_u_j_o_g (n)', 'h_f_i_l_t_e_r (n)','Location','northeast');
        hold off;

% Comparación numérica
disp(['max(abs(hflujog(n)-hfilter(n))) = ' num2str(max(abs(hflujog-hfilter)))])


%% Ejercicio 2
% Filtrado de una señal cardiográfica

clear

% Datos
load ml2r
Ts=3e-3; Fs=1/Ts;
F0=50;
a=0.96;

% Diseño del filtro anti interferencias
w0=2*pi*F0*Ts;

% Diagrama polos-ceros

% Respuesta en amplitud y fase

% Filtrado de las interferencias

% Gráfica



%% Ejercicio 3
% Diseño de filtros con MATLAB

clear

% Datos
Fs=8000; 
Fc=1000;
wc=Fc/(Fs/2);

% FIR

s1=fir1(32, wc);
[H1,w]=freqz(s1,1,512);

% Butterworth

[s2, m2] = butter(6, wc);
[H2,w]=freqz(s2,m2,512);

% Chebyshef tipo 1

[s3,m3]=cheby1(6, 0.5, wc);
[H3,w]=freqz(s3,m3,512);

% Chebyshef tipo 2

[s4,m4]=cheby2(6, 60, wc);
[H4,w]=freqz(s4,m4,512);

% Elíptico

[s5,m5]=ellip(6, 0.5, 60, wc);
[H5,w]=freqz(s5,m5,512);

% Gráficas
plot(F,MH1,'r',F,MH2,'g',F,MH3,'k',F,MH4,'b',F,MH5,'m')
title('Respuesta en amplitud de los diferentes filtros')
xlabel('Hz'), ylabel('dB')
grid, axis([F(1) F(end) -125 5])
legend('FIR', 'Butterworth', 'Chebyshef 1', 'Chebyshef 2', 'Elíptico', 'Location', 'SouthWest')
pause

plot(F,FH1,'r',F,FH2,'g',F,FH3,'k',F,FH4,'b',F,FH5,'m')
grid
title('Respuesta en fase de los diferentes filtros')
xlabel('Hz'), ylabel('rad')
legend('FIR', 'Butterworth', 'Chebyshef 1', 'Chebyshef 2', 'Elíptico', 'Location', 'SouthWest')


%% Comprobación de la relación de filtrado
% Ejercicio 4
%
 
load tds
Fs=8000;
% Extracción de un segmento de tds(n)
n=6750:7915;
x=tds(n);
plot(n,x), title('Segmento de la señal tds(n)')
grid, xlabel('n')
axis tight
pause
 
% Diseño de un filtro de tipo Chebychef tipo 1
N=8; R=0.5;
Fc=500;
wp=Fc/(Fs/2);
[B A]=cheby1(N,R,wp);
 
% Respuesta en frecuencia
M=1e4;
[H F]=freqz(B,A,M,Fs);
MH=20*log10(abs(H));
plot(F,MH), title('20log_{10}(|H(F)|)')
grid, xlabel('Hz'), ylabel('dB')
pause
 
% Filtrado de x(n)
y=filter(B,A,x);
 
% Espectros de x(n) e y(n)
[X F]=freqz(x,1,M,Fs); MX=20*log10(abs(X));
[Y F]=freqz(y,1,M,Fs); MY=20*log10(abs(Y));
 
% Gráfica conjunta de la respuesta en amplitud del filtro y del módulo del
% espectro de las señales de entrada y salida del filtro (todas las
% cantidades en dB)
plot(F,MH,'k',F,MX,'g',F,MY,'r')
grid, xlabel('Hz'), ylabel('dB')
axis([0 Fs/2 -200 150])
title('Comprobación gráfica de la relación 20log_{10}(|Y(F)|)=20log_{10}(|X(F)|)+20log_{10}(|H(F)|)')
legend('|H|','|X|','|Y|', 'Location', 'SouthWest')
pause
 
% Enfocamos la gráfica a la zona de estudio
axis([0 2000 -150 125])
 
% Tomamos medidas con ginput
disp('Pique en las verticales correspondientes a los picos de abs(X) en torno a las frecuencias 700 y 810 Hz')
[Frec Amp]=ginput(6);
 
% Veamos si 20*log10(abs(X))+20*log10(abs(H))=20*log10(abs(Y))
disp('Para el pico en la frecuencia próxima a 700 Hz') 
disp(['20*log10(abs(X))+20*log10(abs(H) = ' num2str(Amp(1)+Amp(3)) ' ;  20*log10(abs(Y)) = ' num2str(Amp(2))])
disp('Para el pico en la frecuencia próxima a 810 Hz') 
disp(['20*log10(abs(X))+20*log10(abs(H) = ' num2str(Amp(4)+Amp(6)) ' ;  20*log10(abs(Y)) = ' num2str(Amp(5))])
 
% ¡NO SON APROXIMADAMENTE IGUALES!


%% Ejercicio 5
% Efecto de la aritmética finita en los filtros digitales
% Datos

clear

% Datos
Fs=8000; Fc=1000;
wc=Fc/(Fs/2);

% Normalización y cuantificación de los coeficientes B y A con 32, 24 y 16 bits
MaxB=max(abs(B)); MaxA=max(abs(A));
B=B/MaxB; A=A/MaxA;
G=MaxB/MaxA;
B32=fix(B*2^32)/2^32; A32=fix(A*2^32)/2^32;
B24=fix(B*2^24)/2^24; A24=fix(A*2^24)/2^24;
B16=fix(B*2^16)/2^16; A16=fix(A*2^16)/2^16;

[sos G]=tf2sos(B,A);

[hZ hP hl]=zplane(B,A);
set(hZ,'Color','r')
set(hP,'Color','r')
hold on
axis([-1.1 1.1 -1.1 1.1])
text(-0.9,0.05,'Secciones de segundo orden')
text(-0.9,-0.05,'cuantificadas con 16 bits')
text(-0.9,-0.15,'en azul, Con 64 bits en rojo')
for n=1:6
    pause(2)
    zplane(Bs(n,:),As(n,:))
    axis([-1.1 1.1 -1.1 1.1])
end
hold off


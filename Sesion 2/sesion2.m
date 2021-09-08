% Sesión 2 del laboratorio de TDSÑ (Solución Parcial)
%
% Sistemas LTI racionales y procesado local de señales

%% Apartado 1.1
%Diseño de un filtro elíptico de orden 6, Rp=0.5 dB, Rs=50dB y w_corte= pi/4.
[B,A]=ellip(6,0.5,50,1/4);

% Presentación por pantalla de los coeficientes
disp(['B =    ' num2str(B)])
disp(['A =    ' num2str(A)])

%% Apartado 1.2
%Diagrama polos-ceros
zplane(B,A)

%% Apartado 1.3
%Respuesta impulsiva
[h,n]=impz(B,A,30); % aquí hay que calcular n y h
stem(n,h,'k','MarkerfaceColor','k'), title('Respuesta impulsiva'),
xlabel('n');ylabel('h[n]')
axis([-0.5 29.5 min(h)-.05 max(h)+.05]), grid

% Medida de la frecuencia de oscilación de h(n)
[x,~]=ginput(5); % aquí hay que calcular x e y
disp(['Frecuencia de oscilación de h(n) = ' num2str(2*pi/(x(2)-x(1))) ' rad'])
disp(['pi/4 = ' num2str(pi/4)])
abs(roots(A))
angle(roots(A))

%% Apartado 1.4    
%Respuesta en frecuencia H(w)
[H, w]=freqz(B,A,n); % aquí hay que calcular H y w
% Dibujo
subplot(3,1,1)
plot(w/pi,20*log10(abs(H)))
grid, axis tight
title('Respuesta en frecuencia'), ylabel('|H(\omega)|_d_B')
subplot(3,1,2)
plot(w,angle(H))% falta código
grid, axis tight
ylabel('ARG[H(e^{j(\omega)})]'), xlabel('w'),% falta código
title('Valor principal de la fase')
subplot(3,1,3)
plot(w,unwrap(angle(H)))
ylabel('arg[H(e^{j(\omega)})]'), xlabel('w'),% falta código
title('Fase continua')
subplot


%% Apartado 1.5
%Filtrado de la vocal /a/
load a  % falta código
n=0:length(a)-1;
y=filter(B,A,a); % falta código
plot(n,a,'g',n,y,'r')
title('Filtrado de la señal correspondiente a la vocal /a/')
xlabel('n')
legend('Entrada','Salida','Location','NorthEast')% falta código
axis([0 n(end) min(a)-100 max(a)+100]), grid

%% Apartado 1.6
%Cálculo de la respuesta impulsiva filtrando el impulso unidad

N =30;% falta código  
n=0:N-1;
d=double(n==0);
h=impz(B,A,n); % falta código 
f=filter(d,1,h);

stem(n,f,'k','MarkerfaceColor','k'), 
title('Respuesta impulsiva con filter'),
xlabel('n');ylabel('h[n]')
axis([-0.5 n(end)+0.5 min(h)-.05 max(h)+.05]), grid




%% Apartado 2.1
% Datos del ejercicio
load tds

player = audioplayer(int16(tds), 8000);
play(player)%para escucharla, buscar ayuda sobre "audioplayer"

Tv=60e-3; Ts=125e-6;
Lv=Tv/Ts;   
Ld=Lv/2;
Ploc=[];

% Inicialización de la función slocal
xl=slocal(tds,Lv,Ld,0);

% Bucle de cálculo de la potencia localizada
while 1
    xl=slocal(tds,Lv,Ld,1);
    if isempty(xl)
        break,   % Se alcanzó el final de la señal tds
    end
    % Procesado local
     %Potencia local calculada como producto escalar
    Ploc(end+1)=1/Lv*(xl'*xl) ; %ponga su código aquí (no confundir 1(uno) con l(ele))
    
end

% Dibujo
n=0:length(tds)-1; % Eje temporal para tds
m=0:length(Ploc)-1; % Eje temporal para Ploc
subplot(2,1,1), plot(n,tds)
title('Señal tds'), xlabel('n'), ylabel('tds[n]'), grid
axis tight
subplot(2,1,2), plot(m,Ploc)
title('Señal potencia localizada'), xlabel('m'), ylabel('Ploc[m]'), grid
axis tight
subplot

%% Apartado 2.2
% Datos

SNR=[-10 0 10 20];
A=5;
L=100;
n=0:L-1;
x=A*cos((2*pi/20)*n);
sumx = 0;
sumr = 0;

% Bucle de cálculo
i=1;  % Índice de la gráfica
for c = 1:1:4 % falta código
    r=randn(1,L);
    for p = 0:L-1        
       sumx= sumx + (A*cos((2*pi/20)*p))^2;
       sumr = sumr + r(p+1).^2;
    end
    g= sqrt(sumx/(10^(SNR(c)/10)*sumr)); %falta código
    y=x+g*r;
    subplot(2,2,i)
    plot(n,y), 
    title(['Señal y[n] con SNR = ' num2str(SNR(c)) ' dB']), 
    xlabel('n')
    i=i+1;
end
subplot

%% Apartado 2.3
% Datos

load x; load r;
Lv=100; Ld=50;

% Cálculo de la energía localizada de la señal x
Plocx=[];
xl=slocal(x,Lv,Ld,0);
while 1
    xl=slocal(x,Lv,Ld,1);
    if isempty(xl)
        break,   % Se alcanzó el final de la señal tds
    end
    Plocx(end+1)=sum(xl'*xl); %Alternativamente
end

% Cálculo de la energía localizada del ruido r(n)
Plocr=[];
rl=slocal(r,Lv,Ld,0);
while 1
    rl=slocal(r,Lv,Ld,1); % segmento local de ruido
    if isempty(rl)
        break,   % Se alcanzó el final de la señal tds
    end
    Plocr(end+1)= sum(rl'*rl); %falta código
end

SNRseg=10*log10(Plocx./Plocr);

SNRglobal=10*log10(sum(x.^2)/sum(r.^2)); %falta código

% Dibujo
n=0:length(x)-1; % Eje temporal para x
m=0:length(SNRseg)-1; % Eje temporal para SNRseg
subplot(3,1,1), plot(n,x)
title('Señal x[n]'), xlabel('n'), grid
axis tight
subplot(3,1,2), plot(n,r)
title('Señal r(n)'), xlabel('n'), grid
axis tight
subplot(3,1,3), plot(m,SNRseg,'g',[0 m(end)],[SNRglobal SNRglobal],'r')
title('Señal SNR segmental'), xlabel('m'), grid
legend('SNR segmental', 'SNR global')
axis tight
subplot

%% Apartado 3.1
%Diseñe un filtro FIR de tipo 1, orden del filtro es 28 y 
%la frecuencia de corte pi/4

h=fir1(28,1/4);
A=1;
B=h;

%% Apartado 3.1.2

zplane(B,A)

%% Apartado 3.1.3

%Respuesta impulsiva
[h,n]=impz(B,A,30); % aquí hay que calcular n y h
stem(n,h,'k','MarkerfaceColor','k'), title('Respuesta impulsiva'),
xlabel('n');ylabel('h[n]')
axis([-0.5 29.5 min(h)-.05 max(h)+.05]), grid

% Medida de la frecuencia de oscilación de h(n)
[x,~]=ginput(5); % aquí hay que calcular x e y
disp(['Frecuencia de oscilación de h(n) = ' num2str(2*pi/(x(2)-x(1))) ' rad'])
disp(['pi/4 = ' num2str(pi/4)])

%% Apartado 3.1.4

%Respuesta en frecuencia H(w)
[H, w]=freqz(B,A,n); % aquí hay que calcular H y w
% Dibujo
subplot(3,1,1)
plot(w/pi,20*log10(abs(H)))
grid, axis tight
title('Respuesta en frecuencia'), ylabel('|H(\omega)|_d_B')
subplot(3,1,2)
plot(w,angle(H))% falta código
grid, axis tight
ylabel('ARG[H(e^{j(\omega)})]'), xlabel('w'),% falta código
title('Valor principal de la fase')
subplot(3,1,3)
plot(w,unwrap(angle(H)))
ylabel('arg[H(e^{j(\omega)})]'), xlabel('w'),% falta código
title('Fase continua')
subplot


%% Apartado 3.1.5

%Filtrado de la vocal /a/
load a  % falta código
n=0:length(a)-1;
y=filter(B,A,a); % falta código
plot(n,a,'g',n,y,'r')
title('Filtrado de la señal correspondiente a la vocal /a/')
xlabel('n')
legend('Entrada','Salida','Location','NorthEast')% falta código
axis([0 n(end) min(a)-100 max(a)+100]), grid

%% Apartado 3.1.6

N =30;% falta código  
n=0:N-1;
d=double(n==0);
h=impz(B,A,n); % falta código 
f=filter(d,1,h);

stem(n,f,'k','MarkerfaceColor','k'), 

title('Respuesta impulsiva con filter'),
xlabel('n');ylabel('h[n]')
axis([-0.5 n(end)+0.5 min(h)-.05 max(h)+.05]), grid

%% Apartado 3.2

% Datos del ejercicio
load tds

Tv=20e-3; Ts=125e-6;
Lv=Tv/Ts;   
Ld=Lv/2;
Ploc=[];

% Inicialización de la función slocal
xl=slocal(tds,Lv,Ld,0);

sum(abs(diff(tds>0)))




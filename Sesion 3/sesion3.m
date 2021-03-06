% Sesi?n 3: Muestreo

% Autor: Roberto Mart?n Luengo

%% Muestreo de sinusoides

% Apartado 1.1

% Datos
Frecuencias=[200 850 1850 3800]; A=1;
Tv=10e-3;
Fs=8000; Ts=1/Fs;

% Rejilla temporal para el muestreo
t=0:Ts:Tv-Ts;

% C?lculo y dibujo de las sinusoides
n=1;
for F0=Frecuencias
    subplot(2,2,n)
    plot(t,A*cos(2*pi*F0*t))
    title(['cos(2\pi' num2str(F0,4) 't)']), xlabel('t'), grid
    n=n+1;
end
pause

% Apartado 1.2
% Datos
Frecuencias=[200 850 1850 5800];
Fs=6000; Ts=1/Fs;

% Rejilla temporal para el muestreo
t=0:Ts:Tv-Ts;

% C?lculo y dibujo de las sinusoides
n=1;
for F0=Frecuencias
    subplot(2,2,n)
    plot(t,A*cos(2*pi*F0*t))
    title(['cos(2\pi' num2str(F0,4) 't)']), xlabel('t'), grid
    n=n+1;
end
pause

% Apartado 1.3

% C?lculo y dibujo de las sinusoides
n=1;
for F0=Frecuencias
    subplot(2,2,n)
    plot(t,A*sin(2*pi*F0*t)),
    title(['sin(2\pi' num2str(F0,4) 't)']), xlabel('t'), grid
    n=n+1;
end

subplot

% a.- ?Por qu? aparece un efecto parecido a una modulaci?n de amplitud?
%      Se modula entre 0 y 1 porque a A le damos el valor 1, asi es
%      parecido a una modulaci?n de amplitud
%
% b.- ?Por qu? hay un cambio de signo de una de las se?ales en los apartados 1.2 y 1.3?
%      El cambio de signo de una de las se?ales en los apartados 1.2 y 1.3
%      se debe a que se esta calculando y dibujando una se?al distinta,
%      concretamente 



%% Interpolaci?n

% Apartado 2
load tds
%wavplay(int16(tds), 24000)
p=audioplayer(int16(tds),24000);
play(p);
pause

% Apartado 2.1
L=3;
tdsL = zeros(L*length(tds), 1);
tdsL(1:L:length(tdsL)) = tds;
%wavplay(int16(tdsL), 24000)
p=audioplayer(int16(tdsL),24000);
play(p);

N=30;
Ncomienzo=1250; 
%n=L*Ncomienzo:L*Ncomienzo+N-1; %si el indexado de vectores comenzase en n=0
n=L*Ncomienzo+1:L*Ncomienzo+N; %Pero en MATLAB los ?ndices comienzan en n=1
stem(n,tdsL(n)), title('Se?al expandida por un factor L=3'), xlabel('n'), grid
pause

% Apartado 2.2
nh = -100:100;
h = L*sin(pi/L*nh)./(pi*nh);
h(isnan(h)) = 1;
plot(nh,h), title('Respuesta impulsiva del filtro interpolador'), xlabel('n'), grid
pause
stem(-20:20,h( (nh>=-20) & (nh<=20) ))
title('Respuesta impulsiva del filtro interpolador en torno a n=0'), xlabel('n'), grid
axis([-20.5 20.5 -0.3 1.1])
pause

%   Explique la operaci?n del tercer comando MATLAB.
%   isnan(h) devuelve una matriz del mismo tama?o que A que contiene 1 l?gicos (verdadero) 
%   donde los elementos de h son NaN simb?licos y 0 l?gicos (falso) donde
%   no lo son. 
%   Es equivalente a un filtro paso bajo

% Apartado 2.3
tdsLi=filter(h,1,tdsL);
%wavplay(int16(tdsLi),24000)
p=audioplayer(int16(tdsLi),24000);
play(p);


% Visualizaci?n de un segmento de 30 puntos de todas las se?ales
Npuntos=30;
Ncomienzo=1250;
%n=L*Ncomienzo:L*Ncomienzo+Npuntos-1; %si el indexado de vectores comenzase en n=0
n=L*Ncomienzo+1:L*Ncomienzo+Npuntos; %Pero en MATLAB los ?ndices comienzan en n=1
nnn=L*Ncomienzo+1:L:L*Ncomienzo+Npuntos;
nn=Ncomienzo+1:Ncomienzo+length(nnn);
[hmax,nmax]=max(h);
N=abs(nh(nmax)-nh(1));
plot(nnn,tds(nn),'b')
hold
stem(n,tdsLi(n+N),'g')
stem(n,tdsL(n),'r')
hold
grid, xlabel('n')
legend('tds(t)','tdsLi[n]','tdsL[n]','Location','NorthWest')

%   Se introduce +N en la octava y d?cima l?neas del c?digo para compensar el retardo
%   de N muestras que introduce el filtro interpolador, ?cu?nto vale ese retardo?.
%   El valor de N es 100. N es la semilongitud del filtro, la mitad de la
%   longitud del filtro -1. El retardo de N es una fraccion exacta de la
%   anchura del filtro

%% Diezmado

% Apartado 3.1

% Fs = 1/(2*Fn) = 1/(2*(2F0)) = 1/4F0 ---> fs1 = 1/(3*4F0) ---> M=3
% El valor de M es 3 porque la Ts que hay es 1/12Fo y la frecuencia m?xima es
% Fs= 1/4Fo; por lo tanto, el factor M que est? multiplicando a 12F0 tiene
% que ser 3.

% Datos
F0=1e3;
T1=-1; T2=1;

% Datos derivados
Ts=1/(12*F0);
t=T1:Ts:T2;

% C?lculo y dibujo de la se?al
x=(sin(2*pi*F0*t)./(pi*t)).^2;
x(isnan(x))=4*F0^2;
n=round(t/Ts);
stem(-30:30,x( (n>=-30) & (n<=30))), title('(sin(2*\piF_0nT_s)/\pinT_s)^2'), xlabel('n'), grid
axis([-30.5 30.5 -0.05*4*F0^2 1.05*4*F0^2])
pause

% Apartado 3.2
[X w]=freqz(x,1,'whole');
plot(w/pi,abs(X)), title('|X(\omega)|'), xlabel('\omega (x \pi rad/muestra)'), grid
axis tight
pause

% Apartado 3.3
L=3;
xe=zeros(L*length(x),1);
xe(1:L:length(xe))=x;
[Xe w]=freqz(xe,1,'whole');
plot(w/pi,abs(Xe)), title('|Xe(\omega)|'), xlabel('\omega (x \pi rad/muestra)'), grid
axis tight
pause

% Aliasing, porque no se cumple el criterio de Nyquist.

% Apartado 3.4
M=2;
xd=x(1:M:length(x));
[Xd w]=freqz(xd,1,'whole');
plot(w/pi,abs(Xd)), title('|Xd(\omega)|'), xlabel('\omega (x \pi rad/muestra'), grid
axis tight


%% Cuantificaci?n

% Apartado 4.1

% Datos
Fs=8000; Ts=1/Fs;
N=8;
Xm= 2^(16-1);% El fondo de escala es el mayor nivel que se puede reproducir
             % de la se?al original que est? cuantificada con 16 bits.
load tds

d = 2*Xm/2^N; % ?ste es delta, el tama?o del escal?n para N bits y Xm
tdsQ = d*round(tds/d); %este es el cuantificador por redondeo al n. m. c.
e= tdsQ-tds;  % la se?al de error
n=0:length(tds)-1;
subplot(3,1,1), plot(n,tds), title('Se?al tds[n]'), grid, xlabel('n')
axis tight
subplot(3,1,2), plot(n,tdsQ), title('Se?al cuantificada tdsQ[n]'), grid, xlabel('n')
axis tight
subplot(3,1,3), plot(n,e), title('Se?al e[n]'), grid, xlabel('n')
axis([0 n(end) -d/2*1.1 d/2*1.1])
subplot

%wavplay([int16(tds)' zeros(1,4000) int16(tdsQ)' zeros(1,4000) int16(e)'],Fs)
p=audioplayer([int16(tds)' zeros(1,4000) int16(tdsQ)' zeros(1,4000) int16(e)'],Fs);
play(p)
pause

% Apartado 4.2
%SNRglobal=10*log10(sum(tds.^2)/sum(e.^2));
SNRglobal=10*log10((tds'*tds)/(e'*e));
fprintf('SNRglobal = %5.2f dB\n',SNRglobal)
pause

% Falta el termino de adaptaci?n

% Apartado 4.3
% Datos
Tv=20e-3;

Lv=Tv/Ts; Lsol= 0.4*Lv; Ld=Lv-Lsol;


% C?lculo de la energ?a localizada de la se?al x
Ploctds=[];
tdsl=slocal(tds,Lv,Ld,0);
while 1
    tdsl=slocal(tds,Lv,Ld,1);
    if isempty(tdsl)
        break,   % Se alcanz? el final de la se?al tds
    end
    Ploctds(end+1)=sum(tdsl.^2)/Lv;
end

% C?lculo de la energ?a localizada del ruido e(n)
Ploce=[];
el=slocal(e,Lv,Ld,0);
while 1
    el=slocal(e,Lv,Ld,1); % segmento local de ruido
    if isempty(el)
        break,   % Se alcanz? el final de la se?al tds
    end
    Ploce(end+1)=sum(el.^2)/Lv;
end

SNRseg=10*log10(Ploctds./Ploce);

t=(0:length(tds)-1)*Ts;
subplot(211), plot(t,tds), title('Se?al original tds'), xlabel('t'), grid
axis tight
tSNRseg=(0:length(SNRseg)-1)*Ld*Ts;
subplot(212), plot(tSNRseg,SNRseg,'g',[tSNRseg(1) tSNRseg(end)],[SNRglobal SNRglobal],'r')
title('SNR segmental en verde y SNRglobal en rojo'), xlabel('t'), grid
axis tight
subplot
pause


%   Deberia de ser una distribuci?n uniforme. Es porque hay un solapamiento
%   del 40%

% Apartado 4.4
M = 50;
[h, b] = hist(e, M);
bar(b, h/(length(e)*d/M));
hold on
plot(b, (1/d)*ones(size(b)),'g')
title('Diagrama de barras del histograma de amplitudes y densidad de distribuci?n de amplitudes te?rica')
xlabel('e')
hold off
pause

%


% Apartado 4.5
ree=xcorr(e)/length(e);
[Mx Ind]=max(ree);
m=-30:30;
subplot(2,1,1)
stem(m,ree(m+Ind))
title('\phi_{ee}[m]'), xlabel('m'), grid

%Estimaci?n del espectro de potencia:
%[Pee w]=periodogram(e);
%Pee=pi*Pee; % Suprimimos las unidades que introduce MATLAB (units of power per radians per sample)

%Estimaci?n alternativa del espectro de potencia:
%(considerando que la autocorrelaci?n es una secuencia determinista)
[Pee w]=freqz(ree,1,'whole');%Transformada de Fourier de la autocorrelaci?n
Pee=abs(Pee); 

Pote=10*log10(d^2/12);
subplot(2,1,2)
plot(w/pi,10*log10(Pee),'g',[w(1) w(end)]/pi,[Pote Pote],'r')
title('Estimaci?n de la densidad espectral de potencia de e[n] en verde, valor te?rico en rojo')
xlabel('\omega (x \pi rad/muestra)'), ylabel('dB'), grid
subplot
axis tight

%   h/(longitud de e/M)*1/delta. Esta dividido por el numero de intervalos.
%   Quieres que valga 1/delta, lo dejas en torno a 1/delta
%   M no es la longitud del vector. h es el numero de muestras


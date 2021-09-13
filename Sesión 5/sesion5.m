
% Sesión 5: La DFT


%% Muestreo de la transformada de Fourier
% Apartado 1.1 y 1.2
% Muestreo con N=L y N=5L

clear

% Datos
L=49;
Nvec=[L, 5*L];
%Ventana rectangular de L muestras
for N=Nvec,
    % Eje de frecuencias
    k=(0:N-1); wk=2*pi*k/N;

    % Muestreo de W(w)
    W=exp(-1i*wk*(L-1)/2).*sin(wk*L/2)./sin(wk/2); W(1)=L;
    
    % Gráficas
    subplot(2,1,1), plot(wk,20*log10(abs(W))), grid
    title(['Muestreo de ventana rectangular (N=',num2str(N),')']);
    ylabel('|W(2\pik/L)|')
    axis ([min(wk),max(wk),-80, max(20*log10(abs(W)))])
    subplot(2,1,2), plot(wk,angle(W)), grid
    xlabel('\omega'), ylabel('<W(2\pik/L)(rad)')
    axis tight
    subplot
    pause(2)
end

%% Apartado 1.3 y 1.4
% Se repita el apartado 1.1 y 1.2 para la ventana de Hamming
% Muestreo con N=L y N=5*L

clear

% Datos
L=49;
Nvec=[L,5*L];
% Ventana de Hamming
n=0:L-1;
h=0.54-0.46*cos(2*pi*n/(L-1));

for N=Nvec
    % Muestreo de W(w)
    H=zeros(1,N);
    for k=0:N-1
        H(k+1)=sum(h.*exp(-1i*2*pi*k*n/N)); %DFT de la ventana de Hamming
                 %código aquí
    end
    % Eje de frecuencias discretas
    k=(0:N-1); wk=2*pi*k/N;
    % Gráficas
    subplot(2,1,1), plot(wk,20*log10(abs(H))), grid
    title(['Muestreo de ventana de Hamming (N=',num2str(N),')']);
    ylabel('|H(2\pik/L)|_d_B') %código aquí
    axis ([min(wk),max(wk),-80, max(20*log10(abs(H)))])
    subplot(2,1,2), plot(wk,angle(H)), grid
    xlabel('\omega'), ylabel('<H(2\pik/L)(rad)')
    axis tight
    subplot
    pause(2)
end

%% Muestreo de la transformada de Fourier con la DFT
% Apartado 2.1 y 2.2
% Muestreo con N=L y N=5L

clear

% Datos
L=49;
Nvec=[L,5*L];
% Ventana Rectangular
w=ones(1,L);
for N=Nvec
    % Muestreo de W(w)
    W=fft(w,N);%%%código aquí
    
    % Eje de frecuencias discretas
    k=(0:N-1); wk=2*pi*k/N;
    
    % Gráficas
    subplot(2,1,1), plot(wk/pi,20*log10(abs(W))), grid 
    title(['Muestreo de ventana rectangular (N=',num2str(N),')']);
    ylabel('|H(2\pik/L)|_d_B')
    axis ([min(wk/pi),max(wk/pi),-80, max(20*log10(abs(W)))])
    subplot(2,1,2), plot(wk/pi,angle(W)), grid
    xlabel('\omega'), ylabel('<H(2\pik/L)(rad)')
    axis tight
    subplot
    pause(2)
end

%% Apartado 2.3 y 2.4
% Se repita el apartado 2.1 y 2.2 para la ventana de Hamming
% Muestreo con N=L y N=5L

clear

% Datos
L=49;
N=L;
Nvec=[L,5*L];
% Ventana de Hamming
n=0:L-1;
h=0.54-0.46*cos(2*pi*n/(L-1)); %%%código aquí

for N=Nvec
    % Cálculo de las muestras de H(w) con la fft
    H=fft(h,N);
    % Eje de frecuencias discretas
    k=(0:N-1); wk=2*pi/N*k;%%%código aquí
    % Gráficas
    subplot(2,1,1), plot(wk,20*log10(abs(H))), grid
    title(['Muestreo de ventana de Hamming (N=',num2str(N),')']);
    ylabel('|H(2\pik/L)|_d_B')
    axis ([min(wk),max(wk),-80, max(20*log10(abs(H)))])
    subplot(2,1,2), plot(wk,angle(H)), grid
    xlabel('\omega'), ylabel('<H(2\pik/L)(rad)')
    axis tight
    subplot
    pause(2)
end

%% Comando fftshift
% Apartado 3.1 y 3.2

clear

% Datos
L=49;
N=5*L;

n=0:L-1;
% Ventana de Hamming
h=0.54-0.46*cos(2*pi*n/(L-1));
% Ventana de rectangular
wr=ones(1,L);
% Eje de frecuencias
k=(0:N-1); wk=2*pi*k/N; wk=wk-pi;
%Tipos de ventana
ventanas={'Rectangular','Hamming'};
w=wr; %se comienza con la ventana rectangular
for ventana=ventanas
    if strcmp('Hamming',ventana{1})
        w=h;
    end
    % Muestreo de W(w)
    W=fft(w,N); W=fftshift(W);
    Menos3dB=max(20*log10(abs(W)))-3;
    
    % Gráficas
    subplot(2,1,1)
    plot(wk,20*log10(abs(W)),'r',[-pi pi],[Menos3dB Menos3dB],'g'), grid
    title(['Ventana ',ventana,' entre [-\pi,\pi]']),
    xlabel('\omega'), ylabel('|W(2\pik/L)|_{dB}')
    axis ([min(wk),max(wk),-30, max(20*log10(abs(W)))])
    
    subplot(2,1,2), plot(wk,angle(W)), grid
    xlabel('\omega'), ylabel('<W(2\pik/L) (rad)')
    axis tight
    subplot
    
    % Medida de la atenuación del lóbulo secundario
    disp('Toque sucesivamente en intersección 3dB con lóbulo, la otra intersección y el máximo del lóbulo secundario')

    [x y]=ginput(3);
    fprintf(['\nVentana ',ventana{1},', atenuación del lóbulo secundario = %6.4f dB'], y(1)-y(3)+3)
    fprintf(['\nVentana ',ventana{1},', ancho del lóbulo principal = %6.4f rad\n\n'], abs(x(2)-x(1)))
end
fprintf('Lea resultados en la consola')

%% Transformada de Fourier de una señal senoidal
% Apartado 4.1, ventana rectangular

clear

% Datos
L=87;
N=9*L;

% Datos derivados
n=0:L-1;
x=10*cos((pi/4)*n);
%Eje de frecuencias
k=(0:N-1); wk=2*pi*k/N; wk=wk-pi;

% Muestras de X(w)
X=fft(x,N); X=fftshift(X); %%%código aquí

% Gráfica
plot(wk/pi,20*log10(abs(X))), grid
title('T. de Fourier de un coseno (enventanado rectangular)')
xlabel('\omega (x \pi) rad/muestra'), ylabel('|X(\omega)|_{dB}')
axis([-1 1 (.95)*min(20*log10(abs(X))) (1.05)*max(20*log10(abs(X)))])

% Medidas de la amplitud y frecuencia
[frec amp]=ginput(1);

% Presentación por pantalla
fprintf('\nFrecuencia del coseno = %6.4f rad\n', frec*pi)
fprintf('Amplitud del coseno = %6.4f\n',2/L*10^(amp/20)) %%%código aquí

%% Apartado 4.2, ventana de Hamming

clear

% Datos
L=87;
N=9*L;

% Datos derivados
n=0:L-1;
x=10*cos(pi/4*n);
k=(0:N-1); wk=2*pi*k/N; wk=wk-pi;

% Muestras de X(w)
X=fft(x.*hamming(L).',N); X=fftshift(X);

% Gráfica
plot(wk/pi,20*log10(abs(X))), grid
title('T. de Fourier de un coseno (enventanado Hamming)')
xlabel('\omega (x \pi) rad/muestra'), ylabel('|X(\omega)|_{dB}')
axis([-1 1 (.95)*min(20*log10(abs(X))) (1.05)*max(20*log10(abs(X)))])

% Medidas de la amplitud y frecuencia
[frec amp]=ginput(1);

% Presentación por pantalla
fprintf('\nFrecuencia del coseno = %6.4f rad\n', frec*pi)
fprintf('Amplitud del coseno = %6.4f\n',2/sum(hamming(L))*10^(amp/20))


%% Transformada de Fourier de una señal exponencial
% Apartado 5.1, ventana rectangular

clear

% Datos
L=87;
N=9*L;

% Datos derivados
n=0:L-1;
x=10*exp(-1i*pi/4*n);
k=(0:N-1); wk=2*pi*k/N; wk=wk-pi;

% Muestras de X(w)
X=fft(x,N); X=fftshift(X);

% Gráfica
plot(wk/pi,20*log10(abs(X))), grid
title('T. de Fourier de una señal exponencial (enventanado rectangular)')
xlabel('\omega (x \pi) rad/muestra'), ylabel('|X(\omega)|_{dB}')
axis([-1 1 10 60])

% Medidas de la amplitud y frecuencia
[frec amp]=ginput(1);

% Presentación por pantalla
fprintf('\nFrecuencia del coseno = %6.4f rad\n', frec*pi)
fprintf('Amplitud del coseno = %6.4f\n',1/L*10^(amp/20))


%% Apartado 5.2

clear

% Datos
L=87;
N=9*L;

% Datos derivados
n=0:L-1;
x=10*exp(-1i*pi/4*n);
k=(0:N-1); wk=2*pi*k/N; wk=wk-pi;

% Muestras de X(w) enventanadas con hamming(L)
X=fft(x.*hamming(L).',N); X=fftshift(X);%%%código aquí

% Gráfica
plot(wk/pi,20*log10(abs(X))), grid
title('T. de Fourier de una señal exponencial (enventanado rectangular)')
xlabel('\omega (x \pi) rad/muestra'), ylabel('|X(\omega)|_{dB}')
axis([-1 1 -20 55])

% Medidas de la amplitud y frecuencia
[frec amp]=ginput(1);

% Presentación por pantalla
fprintf('\nFrecuencia del coseno = %6.4f rad\n', frec*pi)
fprintf('Amplitud del coseno = %6.4f\n',1/sum(hamming(L))*10^(amp/20)) %%%código aquí

%% Análisis espectral con la DFT

% Apartado 6

clear

% Datos
load x1
L=length(x);
Fs=2e4;
N=1e4;

% Datos derivados
k=(0:N/2);
ejeF=Fs*k/N;
% Espectros de x(n)
Xr=fft(x,N); Xr=Xr(1:N/2+1);
Xh=fft(x.*hamming(L),N); Xh=Xh(1:N/2+1);

% Gráficas
subplot(2,1,1), plot(ejeF,20*log10(abs(Xr))), grid
title('Espectro de x(n) con ventana rectangular')
xlabel('F Hz'), ylabel('|X(F)|_{dB}')
subplot(2,1,2), plot(ejeF,20*log10(abs(Xh))), grid
title('Espectro de x(n) con ventana de Hamming')
xlabel('Hz'),ylabel('|X(F)|_{dB}')
% Medidas de las amplitudes y frecuencias

% Con la ventana rectangular se aprecian 3 cosenos (se resuelven dos muy
% próximos en frecuencia) y con la ventana de hamming se aprecia un coseno
% de baja amplitud. Se "pica de izquierda a derecha y de arriba abajo en
% los máximos correspondientes a cosenos. Por tanto, P=4

P=4;
[Frec Amp]=ginput(P);

% Cambio de amplitudes en dB a lineales y compensación del efecto de la
% ventana
Amp(1:3)=2/L*10.^(Amp(1:3)/20);
Amp(4)=2/sum(hamming(L))*10^(Amp(4)/20);

% Presentación de los resultados
fprintf('\n\nNúmero de componentes cosenoidales = %2.0f\n',P)
fprintf('Amplitudes\t\t\t%4.1f\t%4.1f\t%4.1f\t%4.1f\n',Amp)
fprintf('Frecuencias (Hz)\t%4.0f\t%4.0f\t%4.0f\t%4.0f\n',Frec)

%% Análisis espectral con la DFT
% Apartado 7

clear

% Datos
load x2
L=length(x);
Fs=1e3;
N=1e4;

% Datos derivados
k=(0:N-1);
ejeF=(Fs/N)*k;ejeF=ejeF-Fs/2;

% Espectros de x(n)
Xr=fft(x,N); Xr=fftshift(Xr);
Xh=fft(x.*hamming(L),N); Xh=fftshift(Xh);

% Gráficas
subplot(2,1,1), plot(ejeF,20*log10(abs(Xr))), grid
title('Espectro de x(n) con ventana rectangular')
xlabel('F Hz'), ylabel('|X(F)|_{dB}')
subplot(2,1,2), plot(ejeF,20*log10(abs(Xh))), grid
title('Espectro de x(n) con ventana de hamming')
xlabel('F Hz'), ylabel('|X(F)|_{dB}')

% Medidas de las amplitudes y frecuencias

% Con la ventana rectangular se aprecian x cosenos en frecuencias próximas
% a 200 y 350 Hz y con la ventana de hamming se aprecia una exponencial
% de baja amplitud y frecuencia en torno a -250 Hz. Se hace click de
% izquierda a derecha y de arriba abajo en los máximos correspondientes 
% a las componentes. Por tanto, P=2 y Q=1

P=2;
Q=1;
[Frec Amp]=ginput(P+Q);

% Cambio de amplitudes en dB a lineales y compensación del efecto de la
% ventana

Amp(1:2)=2/L*10.^(Amp(1:2)/20);
Amp(3)=2/sum(hamming(L))*10^(Amp(3)/20);

% Presentación de los resultados
fprintf('\n\nNúmero de componentes cosenoidales = %2.0f\n',P)
fprintf('Amplitudes\t\t\t%4.1f\t%4.1f\n',Amp(1:2))
fprintf('Frecuencias (Hz)\t%4.0f\t%4.0f\n',Frec(1:2))

fprintf('\nNúmero de componentes exponenciales = %2.0f\n',Q)
fprintf('Amplitudes\t\t\t%4.1f\n',Amp(3))
fprintf('Frecuencias (Hz)\t%4.0f\n',Frec(3))

%% Análisis espectral con la DFT

clear

% Apartado 8

% Datos
load x3
L=length(x);
Fs=2e4;
N=Fs/2;   % Para que el error de medida sea inferior o igual a 1Hz

% Datos derivados
k=(0:N-1);
ejeF=(Fs/N)*k;ejeF=ejeF-Fs/2;

% Espectros de x(n)
Xr=fft(x,N); Xr=fftshift(Xr);
Xh=fft(x.*hamming(L),N); Xh=fftshift(Xh);
Xk=fft(x.*kaiser(L,12),N); Xk=fftshift(Xk);


% Gráficas
subplot(3,1,1), plot(ejeF,20*log10(abs(Xr))), grid
title('Espectro de x(n) con ventana rectangular'), xlabel('Hz'), ylabel('dB')
subplot(3,1,2), plot(ejeF,20*log10(abs(Xh))), grid
title('Espectro de x(n) con ventana de hamming'), xlabel('Hz'), ylabel('dB')
subplot(3,1,3), plot(ejeF,20*log10(abs(Xk))), grid
title('Espectro de x(n) con ventana de kaiser'), xlabel('Hz'), ylabel('dB')

% Medidas de las amplitudes y frecuencias

% Con la ventana rectangular se resuelven x exponenciales en torno a -500 Hz
% Con la ventana de hamming no se aprecian componentes adicionales.
% Con la ventana de Kaiser (beta=12) se aprecian x cosenos en torno a las 
% frecuencias de 1000 y 9000 Hz. Por tanto, P= y Q=

P=2;
Q=2;
[Frec Amp]=ginput(P+Q);

% Cambio de amplitudes en dB a lineales y compensación del efecto de la
% ventana

Amp(1:2)=2/L*10.^(Amp(1:2)/20);
Amp(3:4)=2/sum(kaiser(L,12))*10.^(Amp(3:4)/20);

% Presentación de los resultados
fprintf('\n\nNúmero de componentes cosenoidales = %2.0f\n',P)
fprintf('Amplitudes\t\t\t%6.4f\t%6.4f\n',Amp(3:4))
fprintf('Frecuencias (Hz)\t%4.0f\t%4.0f\n',Frec(3:4))

fprintf('\nNúmero de componentes exponenciales = %2.0f\n',Q)
fprintf('Amplitudes\t\t\t%4.1f\t%4.1f\n',Amp(1:2))
fprintf('Frecuencias (Hz)\t%4.0f\t%4.0f\n',Frec(1:2))

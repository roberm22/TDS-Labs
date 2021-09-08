% Tratamiento digital señales
% 
% Exploración local de una señal
%
% xl=slocal(x,Lv,Ld,start)
% Entradas:
%    x,  señal que se explora
%    Lv, duración de la señal local
%    Ld, desplazamiento de la señal local
%    start=0 en una llamada inicial antes de comenzar la exploración
% Salidas
%    xl, señal local
%        si xl=[], se ha alcanzado el final de la señal
%         
% Ejemplo de uso
%    Lv=100;
%    Ld=50;
%    xl=slocal(x,Lv,Ld,0);
%    while 1
%        xl=slocal(x,Lv,Ld,1);
%        if isempty(xl)
%            break
%        end
%         . . . EL PROCESADO
%         . . . LOCAL AQUÍ
%    end
%
function xl=slocal(x,Lv,Ld,start)

persistent puntero    % Puntero que marca la muestra de exploración actual

if start==0           % Llamada inicial a la función para la puesta a cero del puntero
    puntero=1;
    xl=0;
    return
end

if (puntero+Lv-1) > length(x)   % Comprobación de que se ha alcanzado el final de la señal
    xl=[];
    return
end

xl=x(puntero:puntero+Lv-1);     % Señal local actual
puntero=puntero+Ld;             % Actualización del puntero
    
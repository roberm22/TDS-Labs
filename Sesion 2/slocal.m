% Tratamiento digital se�ales
% 
% Exploraci�n local de una se�al
%
% xl=slocal(x,Lv,Ld,start)
% Entradas:
%    x,  se�al que se explora
%    Lv, duraci�n de la se�al local
%    Ld, desplazamiento de la se�al local
%    start=0 en una llamada inicial antes de comenzar la exploraci�n
% Salidas
%    xl, se�al local
%        si xl=[], se ha alcanzado el final de la se�al
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
%         . . . LOCAL AQU�
%    end
%
function xl=slocal(x,Lv,Ld,start)

persistent puntero    % Puntero que marca la muestra de exploraci�n actual

if start==0           % Llamada inicial a la funci�n para la puesta a cero del puntero
    puntero=1;
    xl=0;
    return
end

if (puntero+Lv-1) > length(x)   % Comprobaci�n de que se ha alcanzado el final de la se�al
    xl=[];
    return
end

xl=x(puntero:puntero+Lv-1);     % Se�al local actual
puntero=puntero+Ld;             % Actualizaci�n del puntero
    
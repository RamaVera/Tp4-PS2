function plotloc3(RAD, pest, p, map)
%  
% Función: 
% plotloc3(RAD, pest, p, map)
% 
% Descripción: 
% Función para gaficar en 3D la trayectoria estimada.
%
% Parámetros:
% RAD   Matriz de Lx3 cuyas columnas son las coordenadas (px, py, pz) de
%       cada radar y sus filas corresponden a los radares R1, R2, ..., RL.
% pest  Matriz de Nx3, cuyas columnas son las coordenadas (px, py, pz) para
%       cada posición estimada y sus filas son dichas posiciones estimadas  
%       en cada tiempo de muestreo siendo N muestras totales.
% p     Matriz de Nx3, cuyas columnas son las coordenadas (px, py, pz) para
%       cada posición real y sus filas son dichas posiciones para cada 
%       tiempo de muestreo siendo N muestras totales.
% map   Matriz que representa un mapa correspondiente a la superficie sobre
%       la que está sobrevolando la aeronave. Su uso es opcional y para no
%       incluirla debe ingresarse la opción 'none' en lugar de la imagen.

%figure('PaperSize',[24 11])

imsize = size(map);
xmin = -8000;
xmax = 5000;
ymin = -3000;
ymax = 2000;

% Trayectoria 3D

hold on

if strcmp(map,'none')
else
X = linspace(xmin, xmax, imsize(1));
Y = linspace(ymin , ymax, imsize(2));
[X,Y] = meshgrid(X,Y);
s = pcolor(X',Y', map);
set(s, 'EdgeColor', 'none');
colormap('gray')
end

plot3(RAD(:,1), RAD(:,2), RAD(:,3),'yo','LineWidth',2,'MarkerFaceColor','k');
plot3(p(:,1), p(:,2), p(:,3), 'Color',[0,0.9,0], 'LineWidth', 1.2);
plot3(pest(:,1), pest(:,2), pest(:,3), 'r-', 'LineWidth', 1.2);

grid on
box on
ylabel('Y [m]')
xlabel('X [m]')
zlabel('Z [m]')
title('Trayectoria estimada')

if strcmp(map,'none')
legend('Radares','Real', 'FKE')
else
legend('Mapa','Radares','Real', 'FKE')
end
hold off

xlim([xmin xmax])
ylim([ymin ymax])
view([150 30])

end
clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Definiciones Previas Obligatorias
plotInitData = 'false';
plotObservability = 'false';
use_extended_system = false;
use_square_root_algorithm = true;
params.cantRadares = 4;

if params.cantRadares > 8
    Error: 'Cantidad Maxima de Radares'
end
%/-------------------------------------/
%Auxiliares
data = load("tp4_fke");
X = [];
E = [];
I = eye(3);
O = zeros(3);
params.posicionesRadares = data.RP(1:params.cantRadares,:);
params.tiempoFinal = length(data.T);
params.T = data.T(:,1:params.cantRadares)' ;
params.MuMed = 0;
params.VarMed= 3.6*10^-13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parametros del modelo continuo
syms px py pz vx vy vz ax ay az tb
q = data.q;
%u  = 0;
X_sym = [px;py;pz;vx;vy;vz;ax;ay;az];
X_e_sym = [X_sym ; tb];
A = [O,I,O;O,O,I;O,O,O]; %Matriz de Estados
A_e = [A zeros(size(A, 1), 1) ; zeros(1, size(A, 2) + 1)];
Q = [O,O,O;O,O,O;O,O,q]; %Cov ruido de proces
Q_e = [Q zeros(size(Q, 1), 1) ; zeros(1, size(Q, 2) + 1)];
R = eye(params.cantRadares)*params.VarMed; %Cov ruido de medicion

%Linealizacion
y_sym = alinealFunc(params);
y_e_sym = [y_sym(1) + tb; y_sym(2:end)];  %si tb suma el modelo te queda t-tb ojo
C_sym = jacobian(y_sym ,X_sym); 
C_e_sym = jacobian(y_e_sym ,X_e_sym); 

%Discretizacion
h=1;
Ad = expm(A*h);
Ad_e = expm(A_e*h);
Qd = [q*h^5/20,q*h^4/8,q*h^3/6;q*h^4/8,q*h^3/3,q*h^2/2;q*h^3/6,q*h^2/2,q*h];
Qd_e = [Qd zeros(size(Qd, 1), 1) ; zeros(1, size(Qd, 2) + 1)];

%Condiciones Iniciales
x0_0 = zeros(size(X_sym));
x0_0(1:3)=data.p(1,:)';
x0_0(4:6)=data.v(1,:)';
x0_0(7:9)=data.a(1,:)';
x0_0_e = [x0_0 ; 0];

P0_0 = diag([1 1 1 10^3 10^3 10^3 10 10 10]);
P0_0_e = diag([1 1 1 10^3 10^3 10^3 10 10 10 10^-3]);

%System selection
if use_extended_system == true
    A = A_e;
    X_sym = X_e_sym;
    y_sym = y_e_sym;
    C_sym = C_e_sym;
    Q = Q_e;
    P0_0 = P0_0_e;
    x0_0 = x0_0_e;
    Qd = Qd_e;
    Ad = Ad_e;
end

% Genero todas las muestras con ruido
Ym = getMeasurement(params);

%Algoritmo de Kalman
for k = 1:params.tiempoFinal
    
    %Inicializacion
    if k == 1
        X_kminus_kminus = x0_0;
        P_kminus_kminus = P0_0;
    else
        X_kminus_kminus = X_k_k;
        P_kminus_kminus = P_k_k;
    end
    
    %Obtengo muestra del tiempo k
    Yk=Ym(:,k);
     
    %Prediccion
    X_k_kminus = Ad * X_kminus_kminus ;
    P_k_kminus = Ad * P_kminus_kminus * Ad' + Qd;
    
    %Obtengo C
    values.px = X_k_kminus(1);
    values.py = X_k_kminus(2);
    values.pz = X_k_kminus(3);
    if use_extended_system == true
        values.tb = X_k_kminus(10);
    end
    C = double(subs(C_sym,values));
    
    %Test de Observabilidad
    if strcmp(plotObservability,'true')
        L = isObsv(Ad,C);
    end

    %Actualizacion
    if use_square_root_algorithm == true
        R_chol = chol(R);
        P_k_kminus_chol = chol(P_k_kminus);
        Qd_chol = chol(Qd);
        n = length(X_k_kminus);
        p = length(Yk);
        q = size(Qd, 1);
        M = [R_chol' C*P_k_kminus_chol' zeros(p, q) ; ...
            zeros(n, p) Ad*P_k_kminus_chol' Qd_chol' ; ...
            -Yk'*inv(R_chol) X_k_kminus'*inv(P_k_kminus_chol) zeros(1, q)];
        [Q_QRM, R_QRM] = qr(M');
        R_QRM_tranposed = R_QRM';
        Z = R_QRM_tranposed(p + 1 : p + n, p + 1 : p + n);
        W2 = R_QRM_tranposed(end, p + 1 : p + n);
        X_k_k = Z*W2';
        P_k_k = Z*Z';
    else
        K_k =  P_k_kminus * C' * inv( C * P_k_kminus * C' + R);
        X_k_k =  X_k_kminus + K_k * (Yk - double(subs(y_sym,values)) );
        P_k_k = (eye(size(K_k*C)) - K_k*C) * P_k_kminus ;
    end
    
    X = [X (X_kminus_kminus) ];
    E = [E (Yk - double(subs(y_sym,values)) )];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparacion de Estados y Estimaciones
%------------------------------------%
% Graficos de Datos Iniciales
if strcmp(plotInitData,'true')
    newfig
        subplot(1,2,1)
            hold on
            grid on
            plot(p(1,1),p(1,2),'-o','Color','g','MarkerSize',7)
            plot(p(:,1),p(:,2),'-o','Color','r','MarkerSize',5)
            plot(p(end,1),p(end,2),'-o','Color','b','MarkerSize',7)
            plot(0,0,'o','MarkerSize',10);
            plot(RP(:,1),RP(:,2),'rx','MarkerSize',10); 
            legend({'Inicio','Trayectoria Real','Final','Origen','Radares'})
            title('Posiciones de los Radares en Z= 0')
            xlabel('x(t)')
            ylabel('y(t)')
            xlim([-7000 2200])
            ylim([-3000 1000])
        subplot(1,2,2)
            hold on
            grid on
            plot(0,p(1,3),'-o','Color','g','MarkerSize',7)
            plot(p(:,3),'-o','Color','r','MarkerSize',5)
            plot(length(p(:,1)),p(end,3),'-o','Color','b','MarkerSize',7)
            plot(0,0,'o','MarkerSize',10);
            legend({'Inicio','Trayectoria Real','Final','Origen'})
            xlabel('t')
            ylabel('z(t)')
    set(gcf, 'Position', get(0, 'Screensize'));
end
%------------------------------------%  
%Resultados de Correr el Algoritmo
newfig
%Posicion en X
    subplot(3,1,1)
    hold on
    plot(1:params.tiempoFinal,data.p(:,1)','-o','Color','r')
    plot(1:params.tiempoFinal,X(1,:),'-x','Color','k')
    legend({'Real','Estimacion'})
    title('Posicion en x')
%Posicion en Y
    subplot(3,1,2)
    hold on
    plot(1:params.tiempoFinal,data.p(:,2)','-o','Color','r')
    plot(1:params.tiempoFinal,X(2,:),'-x','Color','k')
    legend({'Real','Estimacion'})
    title('Posicion en y')
%Posicion en Z
    subplot(3,1,3)
    hold on
    plot(1:params.tiempoFinal,data.p(:,3)','-o','Color','r')
    plot(1:params.tiempoFinal,X(3,:),'-x','Color','k')
    legend({'Real','Estimacion'})
    title('Posicion en z')
    set(gcf, 'Position', get(0, 'Screensize'));
 saveas(gcf, 'posicion.png')
%------------------------------------%
newfig
%Velocidad en X
    subplot(3,1,1)
    hold on
    plot(1:params.tiempoFinal,data.v(:,1)','-o','Color','r')
    plot(1:params.tiempoFinal,X(4,:),'-x','Color','k')
    legend({'Real','Estimacion'})
    title('Velocidad en x')
%Velocidad en Y
    subplot(3,1,2)
    hold on
    plot(1:params.tiempoFinal,data.v(:,2)','-o','Color','r')
    plot(1:params.tiempoFinal,X(5,:),'-x','Color','k')
    legend({'Real','Estimacion'})
    title('Velocidad en y')
%Velocidad en Z
    subplot(3,1,3)
    hold on
    plot(1:params.tiempoFinal,data.v(:,3)','-o','Color','r')
    plot(1:params.tiempoFinal,X(6,:),'-x','Color','k')
    legend({'Real','Estimacion'})
    title('Velocidad en Z')
    set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'velocidad.png')
%------------------------------------%
newfig
%Aceleracion en X
    subplot(3,1,1)
    hold on
    plot(1:params.tiempoFinal,data.a(:,1)','-o','Color','r')
    plot(1:params.tiempoFinal,X(7,:),'-x','Color','k')
    legend({'Real','Estimacion'})
    title('Aceleracion en x')
%Aceleracion en Y
    subplot(3,1,2)
    hold on
    plot(1:params.tiempoFinal,data.a(:,2)','-o','Color','r')
    plot(1:params.tiempoFinal,X(8,:),'-x','Color','k')
    legend({'Real','Estimacion'})
    title('Aceleracion en y')
%Aceleracion en Z
    subplot(3,1,3)
    hold on
    plot(1:params.tiempoFinal,data.a(:,3)','-o','Color','r')
    plot(1:params.tiempoFinal,X(9,:),'-x','Color','k')
    legend({'Real ','Estimacion'})
    title('Aceleracion en Z')
    set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'aceleracion.png')
%------------------------------------%
%Trayectoria
newfig
    hold on
    plot(data.p(:,1),data.p(:,2),'-o','Color','r')
    plot(X(1,:),X(2,:),'x','Color','k')
  %  plot(p(:,1)' + etha(1,:),p(:,2)' + etha(2,:), 'linestyle', 'none', ...
  %     'marker', '.')
    legend({'Real','Estimacion'})
    title('Trayectoria X e Y')
    set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'trayectoria.png')
%------------------------------------%\
% Trayectoria en tres dimensiones
newfig
    hold on
    grid on
    plot3(data.p(1,1),data.p(1,2),data.p(1,3),'-o','Color','g','MarkerSize',7)
    plot3(data.p(:,1),data.p(:,2),data.p(:,3),'-o','Color','r','MarkerSize',5)
    plot3(X(1,:),X(2,:),X(3,:),'x','Color','k','MarkerSize',5)
    plot3(data.p(end,1),data.p(end,2),data.p(end,3),'-o','Color','b','MarkerSize',7)
    plot3(0,0,0,'ko','MarkerSize',10);
    plot3(params.posicionesRadares(:,1),params.posicionesRadares(:,2),...
        params.posicionesRadares(:,3),'rx','MarkerSize',10); 
    xlabel('x(t)')
    ylabel('y(t)')
    zlabel('z(t)')
    legend({'Inicio','Trayectoria Real','Trayectoria Estimada','Final','Origen','Radares'})
    view(-115,42)
    set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'trayectoriaRealCompleta.png')

%------------------------------------%
% Trayectoria en tres dimensiones ECM
delta = (data.p - X(1:3, :)').^2;
ecm = mean(sqrt(sum(delta, 2)));
fprintf('Error cuadratico medio de la trayectoria %.2f\n', ecm)
%------------------------------------%
% Trayectoria sobre el mapa
%newfig
%    plotloc3(params.posicionesRadares,X(1:3,:),data.p(:,1:3)',data.map)
%    set(gcf, 'Position', get(0, 'Screensize'));
    
%Autocorrelacion de la innovacion
newfig
    if rem(size(E,1), 2) == 0
        for i=1:size(E,1)
            subplot(size(E,1) / 2, 2 ,i)
            [c,lags] = xcov(E(i,:));
            stem(lags,c)
            title_str = sprintf('Autocorrelacion de la innovacion del radar %d', i);
            title(title_str);
        end
    else
        for i=1:size(E,1)
            subplot(size(E,1),1,i)
            [c,lags] = xcov(E(i,:));
            stem(lags,c)
            title_str = sprintf('Autocorrelacion de la innovacion del radar %d', i);
            title(title_str);
        end
    end
    set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'innovacion.png')
%------------------------------------%
% Bias in radar one measurements
if use_extended_system == true
    newfig
    
    yyaxis left
    y_axis = X(10,:)*1e9;
    x_axis = 1:params.tiempoFinal;
    hPlot = plot(x_axis,y_axis,'-x','Color','k');
    title('Bias en la medición de tiempos del radar 1')
    xlabel('k')
    ylabel('tiempo [ns]')
    
    [max_val,index] = max(x_axis);
    h = gcf;
    cursorMode = datacursormode(h);
    hDatatip = cursorMode.createDatatip(hPlot);
    pos = [x_axis(index) y_axis(index) 0];
    set(hDatatip, 'Position', pos)         
    updateDataCursors(cursorMode)
    
    yyaxis right
    plot(1:params.tiempoFinal,data.T(:,1)*1e6','-o','Color','r')
    xlabel('k')
    ylabel('tiempo [us]')
    
    legend('Bias', 'Mediciones del radar 1')
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, 'bias.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Funciones Invocadas en el programa
%------------------------------------% 
function newfig()
    persistent i;
    if isempty(i)
        i = 0;
    end
    i=i+1;
    figure(i);
end
%------------------------------------% 
function L = isObsv(Ad,C)
    L = rank(obsv(Ad,C));
    if (L ~= length(Ad))
        disp('El sistema NO es completamente observable')
        disp(['La cantidad de estados observables es: ',num2str(L)])
    %disp('El sistema es completamente observable')
    %else
    end 
end
%------------------------------------% 
function y = alinealFunc(params)
    cantRadares = params.cantRadares;
    posicionesRadares = params.posicionesRadares;
    global px py pz
    syms px py pz
        c = 3*10^8;
        for i=1:cantRadares
            px_r = posicionesRadares(i,1);
            py_r = posicionesRadares(i,2);
            pz_r = posicionesRadares(i,3);
            %Y siempre es positiva no hace falta modulo
            y(i,:) = (2/c)*sqrt((px-px_r)^2 + (py-py_r)^2 +(pz-pz_r)^2);
        end
end 
%------------------------------------% 
function Ym = getMeasurement(params)
    cantRadares = params.cantRadares;
    tiempoFinal = params.tiempoFinal;
    mu=params.MuMed;
    sigma=sqrt(params.VarMed);
    etha=normrnd(mu,sigma,cantRadares,tiempoFinal);
    Ym = params.T + etha;
    %Ym(1,:) = Ym(1,:) + 1e-6; %Manual bias for testing purpose
end
% Author: Christopher Reyes
% October 2024
close all;
clc;

% Definición de la función Lagrangiana, su gradiente (Jacobiano) y Hessiano
L = @(u1,u2) 10*u1.^2 + u2.^2; % Función Lagrangiana
g = @(u) [20*u(1); 2*u(2)]; % Gradiente de la función (Jacobiano)
G = @(u) [[20, 0]; [0, 2]]; % Matriz Hessiana de la función

% Inicialización de variables
H = eye(2,2); % Matriz Hessiana inicial (identidad 2x2)

num_iter = 10; % Número de iteraciones del algoritmo
tol = 1e-6; % Tolerancia para el criterio de convergencia
uk = zeros(2, num_iter); % Matriz que almacena los valores de u_k en cada iteración
gk = zeros(2, num_iter); % Matriz que almacena los gradientes en cada iteración
sk = zeros(2, num_iter); % Matriz que almacena las direcciones de búsqueda en cada iteración

Hk = H; % Inicialización de la matriz Hessiana 

% Punto inicial u_0
uk(:,1) = [2,2];%[0.1; 1.0]; 

fprintf('u(1) = \n');
disp(uk(:,1)); % Valor inicial de u

% Configuración y creación de la gráfica de la función Lagrangiana
figure
fsurf(L, [-2 2], 'FaceAlpha', 0.5, 'ShowContours', 'on') % Superficie de la función
hold on 

% Parámetros para la búsqueda de línea
rho = 0.25;
sigma = 0.5;
tau1 = 9;
tau2 = 0.1;
tau3 = 0.5;

% Bucle principal de optimización
for k = 1:num_iter - 1
    gk(:,k) = g(uk(:,k)); % Evaluación del gradiente en el punto u_k

    % Verifica si el criterio de convergencia se cumple (norma del gradiente)
    if norm(gk(:,k)) < tol
        fprintf('Convergencia alcanzada en la iteración %d\n', k);
        break; % Detiene el bucle si el gradiente es suficientemente pequeño
    end

    sk(:,k) = -Hk(:,:)*gk(:,k); % Cálculo de la dirección de búsqueda con la aproximación de la Hessiana

    % Definición de f(alpha) y df(alpha) para la búsqueda de línea
    f_alpha = @(alpha) L(uk(1,k) + alpha*sk(1,k), uk(2,k) + alpha*sk(2,k)); % Función unidimensional de alpha
    df_alpha = @(alpha) transpose(g(uk(:,k) + alpha*sk(:,k)))*sk(:,k); % Derivada de f respecto a alpha

    % Cálculo de alpha usando la búsqueda de línea
    num_iter_line_search = 10;
    alpha0 = 0;
    alpha1 = 1;
    alpha = line_search(num_iter_line_search, alpha0, alpha1, rho, sigma, tau1, tau2, tau3, f_alpha, df_alpha); % Obtiene el alpha óptimo 

    disp(alpha) % Alpha obtenido

    % Actualización del nuevo punto u_{k+1}
    uk(:,k+1) = uk(:,k) + alpha*sk(:,k); 
    
    % Cálculo del nuevo gradiente en el punto u_{k+1}
    gk(:,k+1) = g(uk(:,k+1)); 

    % Variables intermedias para la actualización de la Hessiana
    deltak = uk(:,k+1) - uk(:,k); % Cambio en u
    gammak = gk(:,k+1) - gk(:,k); % Cambio en el gradiente

    % Actualización de la matriz Hessiana
    vk = deltak - Hk*gammak;
    Hk = Hk + vk*transpose(vk)/(transpose(vk)*gammak);
   
    % Punto nuevo calculado
    fprintf('u(%d) = \n', k+1);
    disp(uk(:,k+1));

    % Se grafica Dibuja el nuevo punto
    plot3([uk(1,k) uk(1,k+1)], [uk(2,k) uk(2,k+1)],...
        [L(uk(1,k),uk(2,k)) L(uk(1,k+1),uk(2,k+1))], 'k-o', 'LineWidth', 2)
end

% Función de búsqueda de línea
function alpha = line_search(num_iter, alpha0, alpha1, rho, sigma, tau1, tau2, tau3, f, df)
    alpha = zeros(1, num_iter); % Inicialización del vector de alphas
    a = zeros(1, num_iter); % Inicialización de a
    b = zeros(1, num_iter); % Inicialización de b
    alpha(1) = alpha0; % Primer valor de alpha
    alpha(2) = alpha1; % Segundo valor de alpha
    fmin = 0; % Valor mínimo de f
    mu = (fmin - f(0)) / (rho * df(0)); % Cálculo inicial de mu

    % Bracketing
    for iter = 2:num_iter
        val = f(alpha(iter)); % Evaluación de f(alpha_i)
        if val <= fmin
            alpha = alpha(iter); % Si encuentra un valor que satisface la condición, devuelve alpha
            return;
        end
        if (val > f(0) + alpha(iter) * rho * df(0)) || (f(alpha(iter)) >= f(alpha(iter - 1)))
            a(iter) = alpha(iter - 1); % Actualización de los límites de alpha
            b(iter) = alpha(iter);
            break;
        end
        val2 = df(alpha(iter)); % Evaluación de la derivada df(alpha_i)
        if abs(val2) <= -sigma * df(0)
            alpha = alpha(iter); % Si la derivada cumple con la condición, devuelve alpha
            return;
        end
        if val2 >= 0
            a(iter) = alpha(iter - 1); % Actualización de los límites con base en la derivada
            b(iter) = alpha(iter);
            break;
        end
        if mu <= (2 * alpha(iter) - alpha(iter - 1))
            alpha(iter + 1) = mu; % Actualización del valor de alpha
        else
            a_inter = 2 * alpha(iter) - alpha(iter - 1);
            b_inter = min(mu, alpha(iter) + tau1 * (alpha(iter) - alpha(iter - 1)));
            dfz = (b_inter - a_inter) * df(alpha(iter));
            zmin = -dfz / (2 * (f(b_inter) - f(a_inter) - dfz));
            alpha(iter + 1) = a_inter + zmin * (b_inter - a_inter); % Cálculo del nuevo alpha
        end
    end

    % Sectioning
    for iter2 = iter:num_iter
        a_inter = a(iter2) + tau2 * (b(iter2) - a(iter2));
        b_inter = b(iter2) - tau3 * (b(iter2) - a(iter2));
        dfz = (b_inter - a_inter) * df(alpha(iter2));
        zmin = -dfz / (2 * (f(b_inter) - f(a_inter) - dfz));
        alpha(iter2) = a_inter + zmin * (b_inter - a_inter);
        
        % Evaluación de f(alpha_j)
        val1 = f(alpha(iter2));
        if (val1 > f(0) + rho * alpha(iter2) * df(0)) || (f(alpha(iter2)) >= f(a(iter2)))
            a(iter2 + 1) = a(iter2);
            b(iter2 + 1) = alpha(iter2);
        else
            val2 = df(alpha(iter2));
            if abs(val2) <= -sigma * df(0)
                alpha = alpha(iter2);
                return;
            end
            a(iter2 + 1) = alpha(iter2);
            if (b(iter2) - a(iter2)) * df(alpha(iter2)) >= 0
                b(iter2 + 1) = a(iter2);
            else
                b(iter2 + 1) = b(iter2);
            end
        end
    end
    alpha = alpha(iter2); % Se devuelve el mejor alpha encontrado
end

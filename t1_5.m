function t1_5()

  addpath(['uteis']);

  available_graphics_toolkits();
  graphics_toolkit gnuplot;
  pkg load statistics;


  %Item 1
  FDP_Pareto = @(x, sigma, mu, mi) (1./sigma).*(1 + (mu.*(x - mi))./sigma).^ (-1./mu - 1);

  %Item 2
  x = linspace(0, 5, 100);

  sigma = 1;     % Parâmetro de escala
  mu_values = [1, 5, 20];  % Valores diferentes para o parâmetro de forma
  mi = 0;        % Parâmetro de localização
  cores = ['r', 'g', 'b'];

  figure(8);
  hold on;

  for i = 1:length(mu_values)
      mu = mu_values(i);
      y = FDP_Pareto(x, sigma, mu, mi);
      plot(x, y, 'LineWidth', 2, 'Color', cores(i));
  end

  sigma = 2;

  for i = 1:length(mu_values)
      mu = mu_values(i);
      y = FDP_Pareto(x, sigma, mu, mi);
      plot(x, y, 'LineWidth', 2, 'Color', cores(i), 'LineStyle', '--');
  end



  leg = {'sigma=1, eps=1','sigma=1, eps=5','sigma=1, eps=20','sigma=2, eps=1','sigma=2, eps=5','sigma=2, eps=20'};
  legend(leg, 'Location', 'northeast');
  xlabel('x');
  ylabel('FDP Pareto');
  title('Distribuição Generalizada de Pareto - Função de Densidade de Probabilidade');
  hold off;

  %item 3
  FDA_Pareto = @(x, sigma, mu, mi) 1 - (1 + (mu.*(x - mi))./sigma).^(-1./mu);

  figure(9);
  hold on;

  sigma = 1;

  for i = 1:length(mu_values)
      mu = mu_values(i);
      y = FDA_Pareto(x, sigma, mu, mi);
      plot(x, y, 'LineWidth', 2, 'Color', cores(i));
  end

  sigma = 2;

  for i = 1:length(mu_values)
      mu = mu_values(i);
      y = FDA_Pareto(x, sigma, mu, mi);
      plot(x, y, 'LineWidth', 2, 'Color', cores(i), 'LineStyle', '--');
  end

  leg = {'sigma=1, eps=1','sigma=1, eps=5','sigma=1, eps=20','sigma=2, eps=1','sigma=2, eps=5','sigma=2, eps=20'};
  legend(leg, 'Location', 'northwest');
  xlabel('x');
  ylabel('FDA Pareto');
  title('Distribuição Generalizada de Pareto - Função de Densidade Acumulada');

  hold off;


  %Item 4

n = 1000;  % Tamanho da amostra
mu = 2;  % Parâmetro de forma da distribuição
sigma = 1;  % Parâmetro de escala da distribuição
mi = 0;  % Parâmetro de localização da distribuição

% Gerar a amostra aleatória
x = gprnd(mu, sigma, mi, [n 1]);

% Estimar os parâmetros da distribuição
params = gpfit(x);

% Exibir os parâmetros estimados
disp("Parâmetros estimados:");
disp(["mu = " num2str(params(1))]);
disp(["sigma = " num2str(params(2))]);


endfunction

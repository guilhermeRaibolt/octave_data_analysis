function t1_4()
  addpath(['uteis']);

  available_graphics_toolkits();
  graphics_toolkit gnuplot;
  pkg load statistics;

  disp('                                                                                        ==== Item 1.4.1 ====');
  fdp = @(x,lamb,k) k./lamb.*(x./lamb).^(k-1).*exp(-(x./lamb).^k)

  disp('                                                                                        ==== Item 1.4.2 ====');
  lamb = 2;  % Parâmetro de forma
  k = 1;  % Parâmetro de escala

  x = linspace(0, 5, 10);

  pdf_anonima = fdp(x, lamb, k);
  pdf_wblpdf = wblpdf(x, lamb, k);

  disp('Diferença entre as FDPs:');
  disp(pdf_anonima-pdf_wblpdf);

  x0 = 0;
  xmax = 2.5;
  h = 0.01;
  x = x0:h:xmax;

  disp('                                                                                        ==== Item 1.4.3 ====');
  param = {
    {1.0, 0.5},
    {1.0, 1.0},
    {1.0, 1.5},
    {1.0, 5.0},
    {2.0, 1.0},
    {2.0, 5.0}
  };

  figure();
  titulo = 'Distribuicao de Weibull: FDP';
  ymax = 2.5;
  label_y = 'p(x)';

  leg = plot_param(fdp, x, param, titulo, ymax,label_y);
  legend(leg, 'Location', 'northeast');

  disp('                                                                                        ==== Item 1.4.4 ====');

  fda = @(x,lamb,k) 1 - exp(-(x./lamb).^k);

  figure();
  titulo = 'Distribuicao de Weibull: FDA';
  ymax = 1;
  label_y = 'F(x)';

  leg = plot_param(fda,x,param,titulo,ymax,label_y);
  legend(leg, 'Location', 'southeast');

  disp('                                                                                        ==== Item 1.4.5 ====');
  phi = @(k, x) (mean(x.^k.*log(x))/mean(x.^k)-mean(log(x)))^-1;

  lamb_real = 1.0;
  k_real = 1.0;
  n = 1000;

  x = wblrnd(lamb_real, k_real, [n , 1]);

  disp('                                                                                        ==== Item 1.4.6 ====');
  k = 1;
  for i = 1:100
    k = phi(k, x);
  endfor
  lamb = mean(x.^k).^(1./k);

  fprintf('--Estimativa dos Parametros da Distribuicao de Weibull--\n');
  fprintf('Valor de k real: %.1f\n',k_real)
  fprintf('Valor de lambda real: %.1f\n',lamb_real)
  fprintf('Valor de k estimado: %.4f\n',k)
  fprintf('Valor de lambda estimado: %.4f\n',lamb)

  params = wblfit(x);
  lamb_wblfit = params(1);
  k_wblfit = params(2);

  fprintf('--Comparação dos Valores Estimados--\n');
  fprintf('Valor de lamb estimado (phi): %.4f\n', lamb);
  fprintf('Valor de k estimado (phi): %.4f\n', k);
  fprintf('Valor de lamb estimado (wblfit): %.4f\n', lamb_wblfit);
  fprintf('Valor de k estimado (wblfit): %.4f\n', k_wblfit);

  disp('                                                                                        ==== Item 1.4.7 ====');
  figure();
  hold on;

  k = ceil(1 + 3.322 * log(n));
  [nn_hist xx_hist] = hist(x, k);
  delta_x = xx_hist(2) - xx_hist(1);
  [yy_hist xx_hist] = hist(x, k, 1/delta_x);
  hist(x, k, 1/delta_x)

  xx_fit = 0:0.01:8;
  plot(xx_fit, wblpdf(xx_fit, lamb_wblfit, k_wblfit), 'Color','r', 'LineWidth', 2);

  y_estimado = wblpdf(xx_hist, lamb_wblfit, k_wblfit);
  rrs = sumsq(y_estimado - yy_hist)

  fprintf('--Ajuste de Bondade (Soma Residual de Quadrados)--\n');
  fprintf('Soma Residual de Quadrados: %.4f\n', rrs);
  fprintf('Valor de lamb estimado (wblfit): %.4f\n', lamb_wblfit);
  fprintf('Valor de k estimado (wblfit): %.4f\n', k_wblfit);

  xlabel('Valores');
  ylabel('Frequência');
  title('Histograma com 24 barras');

  hold off;
  disp('                                                                                        ==== Item 1.4.8 ====');
  mu = -2;
  sigma = 0.7;

  figure();
  hold on;

  hist(x, k, 1/delta_x)
  mu_estimate = 1/n * sum(x)
  sigma_estimate = sqrt(1/(n-1)*sumsq(x - mu_estimate))
  printf("|mu_estimate - mu| = %f\n",  abs(mu_estimate - mu));
  printf("|sigma_estimate - sigma| = %f\n",  abs(sigma_estimate - sigma));
  plot(xx_fit, normpdf(xx_fit, mu_estimate, sigma_estimate), 'Color','r', 'LineWidth', 2);
  yy_hat = normpdf(xx_hist, mu_estimate, sigma_estimate);
  rss = sumsq(yy_hat - yy_hist)
  hold off;


endfunction

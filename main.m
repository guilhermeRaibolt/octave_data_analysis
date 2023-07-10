addpath(['./' 'Support_functions']);

% ---- Item 1.3.2.1 ----
pkg load statistics;
warning('off', 'Octave:shadowed-function');



x = -3:0.1:3;
mu = 0;
sigma = 1;

% ---- Item 1.3.2.2 ----
FDP_Normal = @(x, mu, sigma) (1 / sqrt(2 * pi * sigma^2)) * exp(-((x - mu).^2) / (2 * sigma^2));


fdp_caseira = FDP_Normal(x, mu, sigma);
fdp_sistema = normpdf(x, mu, sigma);

diferenca = fdp_caseira - fdp_sistema;

disp('                                                                                        ==== Item 1.3.2.3 ====')
disp(diferenca);




u_values_fdp = [0, 0, 0, -2];
sigma2_values_fdp = [0.2, 1, 5, 0.5];


x = -10:0.1:10;

disp('                                                                                        ==== Item 1.3.2.4 ====')

figure;

% Loop para plotar os gráficos
hold on;
for i = 1:4
    fdp = FDP_Normal(x, u_values_fdp(i), sqrt(sigma2_values_fdp(i)));

    plot(x, fdp, 'LineWidth', 2);
end
hold off;


title('Distribuição Normal (Gaussiana): FDP ');
xlabel('x');
ylabel('FDP');
legend('u=0, sigma^2=0.2', 'u=0, sigma^2=1', 'u=0, sigma^2=5', 'u=-2, sigma^2=0.5');
grid on;


disp('                                                                                        ==== Item 1.3.2.5 ====')



FDA_Normal = @(x, mu, sigma) 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2))));

% Valores de exemplo para x, mu e sigma
x_values = [-1, 0, 1];
mu = 0;
sigma = 1;

fda_caseira = FDA_Normal(x_values, mu, sigma);


fda_sistema = normcdf(x_values, mu, sigma);


diferenca = fda_caseira - fda_sistema;


disp(diferenca);


u_values_fda = [0, 0, 0, -2];
sigma2_values_fda = [0.4, 1, 2.2, 0.7];




x = -10:0.1:10;


figure;

% Loop para plotar os gráficos
hold on;
for i = 1:4
    fda = FDA_Normal(x, u_values_fda(i), sqrt(sigma2_values_fda(i)));
    plot(x, fda, 'LineWidth', 2);
end
hold off;

title('Distribuição Normal (Gaussiana): FDA ');
xlabel('x');
ylabel('FDA');
legend('u=0, sigma^2=0.4', 'u=0, sigma^2=1', 'u=0, sigma^2=2.2', 'u=-2, sigma^2=0.7');
grid on;


disp('                                                                                        ==== Item 1.3.2.6 ====')


mu = -2;
sigma = sqrt(0.5);

percentil = 0.99;

# Manipulando a equação obtemos:
x = mu + sigma * sqrt(2) * erfinv(2 * percentil - 1);

disp(x);




disp('                                                                                        ==== Item 1.3.2.7 ====')

FDP_Normal = @(x, mu, sigma2) (1 / sqrt(2 * pi * sigma2)) * exp(-((x - mu).^2) / (2 * sigma.^2));
% Parâmetros da distribuição normal
mu = -2;
sigma2 = 0.5;
sigma = 0.5;
% Limites de integração
a = mu - 10 * sqrt(sigma2);
b = -2 + erfinv(49/50);

% Valor real da integral
valor_real = normcdf(b, mu, sqrt(sigma2)) - normcdf(a, mu, sqrt(sigma2));

% Matriz de erros
erros = zeros(12, 4);

func = @(x) FDP_Normal(x, mu, sigma)

for n = 1:12

    resultado_trap = integralTrapeziosRepetidaFunc(func, a, b, n, 0);
    erro_trap = abs(valor_real - resultado_trap);
    erros(n, 1) = erro_trap;

    printf("Integral Simpson Repetida:\n");
    resultado_simp = integralSimpsonRepetidaFunc(func, a, b, n, 0);
    erro_simp13 = abs(valor_real - resultado_simp);
    erros(n, 2) = erro_simp13;

    printf("Integral Simpson 3/8 Repetida:\n");
    resultado_simp38 = integralSimpson38RepetidaFunc(func, a, b, n, 0);
    erro_simp38 = abs(valor_real - resultado_simp38);
    erros(n, 3) = erro_simp38;

    printf("Integral Quadratura Gaussiana:\n");
    C = coefGaussLegendre(n + 1);
    [T, A] = tabelaAbcissasPesosGaussLegendre(C);
    resultado_gauss = integralGaussLegendreFunc(func, a, b, n, T, A, 0);
    erro_gauss = abs(valor_real - resultado_gauss);
    erros(n, 4) = erro_gauss;
end

% Impressão da tabela de erros
printf("\nErros dos métodos diferentes de integração\n");
printf("Subdivisões  Trapézios       Simpson 1/3     Simpson 3/8     Quad. Gauss.\n");
printf("---------------------------------------------------------------------\n");
for n = 1:12
    printf("%-13d  %-13.6e  %-13.6e  %-13.6e  %-13.6e\n", n, erros(n, 1), erros(n, 2), erros(n, 3), erros(n, 4));
end

% Plotagem da evolução dos erros em um gráfico de barras
figure;
bar(1:12, erros);
legend('Trapezios', 'Simpson 1/3', 'Simpson 3/8', 'Quadratura Gaussiana');
xlabel('Quantidade de Subdivisões');
ylabel('Erro');
title('Evolução do Erro em Relação à Quantidade de Subdivisões');





disp('                                                                                        ==== Item 1.3.2.8 ====')

% Definir a semente para reproduzibilidade
semente = 2023;
randn("seed", semente);

% Gerar as amostras
n = 1000;
amostras = randn(1, n) * 0.7 - 2;

fprintf('x1: %.4f\n', amostras(1));
fprintf('x2: %.4f\n', amostras(2));
fprintf('x3: %.4f\n', amostras(3));



disp('                                                                                        ==== Item 1.3.2.9 ====')
mu_estimado = sum(amostras) / n;

var_estimada = sum((amostras - mu_estimado).^2) / n;

% Calcular o desvio padrão σ estimado
sigma_estimado = sqrt(var_estimada);

% Calcular a diferença em relação aos valores verdadeiros
mu_diferenca = mu_estimado - (-2);
sigma_diferenca = sigma_estimado - 0.7;

% Imprimir os resultados
fprintf('Estimativa de mu: %.4f\n', mu_estimado);
fprintf('Estimativa de sigma: %.4f\n', sigma_estimado);
fprintf('Diferença em relação a mu verdadeiro: %.4f\n', mu_diferenca);
fprintf('Diferença em relação a sigma verdadeiro: %.4f\n', sigma_diferenca);




disp('                                                                                        ==== Item 1.3.2.10 ====')

mu_verdadeiro = -2
sigma_verdadeiro = 0.7;

FDA_Normal = @(x, mu, sigma) 0.5 * (1 + erf((x - mu) / (sqrt(2) * sigma)));

% Valor de probabilidade desejado (99%)


funcao = @(x) FDA_Normal(x, mu, sigma) - 0.99;
bisec = raizBisecPosFalsa( 0, funcao, -1, 0, 0, 10000);
solucao_analitica = mu + sigma * sqrt(2) * erfinv(2 * 0.99 - 1);

%fprintf('\nSolução analítica: %.4f',solucao_analitica);

fprintf('\nSolução com bisec: %.4f\n',bisec);
fprintf('Solução analítica: %.4f\n',solucao_analitica);


disp('                                                                                        ==== Item 1.3.3.1 ====')


% Usando regra de Sturges
num_barras = ceil(1 + log2(n));
fprintf('Quantidade de barras pela regra de Sturges: %.0f\n',num_barras);

figure;
hold on;

[frequencias, bordas] = hist(amostras, num_barras);  % Calcular o histograma

dist_barras = bordas(2) - bordas(1);  % Largura da barra
[y_normalizado x_normalizado] = hist(amostras, num_barras, 1/dist_barras);
hist(amostras, num_barras, 1/dist_barras);
xlabel('Valores');
ylabel('Frequência Normalizada');
title('Histograma Normalizado');

disp('                                                                                        ==== Item 1.3.3.2 ====')

x = -7:0.01:3;
plot(x, normpdf(x, mu_estimado, sigma_estimado));
y_previsto = normpdf(-7:1:3, mu_estimado, sigma_estimado);
RSS = sumsq(y_previsto - y_normalizado)
fprintf('RSS: %.4f\n',RSS);



t1_4();
t1_5();





disp('                                                                                        ==== Item 2.1 ====')
nome_arquivo = '41002h2022.txt';


disp('                                                                                        ==== Item 2.3 ====')
fid = fopen(nome_arquivo, 'r');


cabecalho = fgetl(fid);


velocidade_vento = [];
altura_ondas = [];


while ~feof(fid)
    linha = fgetl(fid);


    if isempty(linha)
        break;
    end

    valores = strsplit(linha, ' ');

    velocidade_vento = [velocidade_vento; str2double(valores{7})];
    altura_ondas = [altura_ondas; str2double(valores{9})];
end


fclose(fid);

save('velocidade_vento.bin', 'velocidade_vento', '-mat');
save('altura_ondas.bin', 'altura_ondas', '-mat');









disp('                                                                                        ==== Item 2.4 ====')

tamanho_vento_original = length(velocidade_vento);
tamanho_ondas_original = length(altura_ondas);


altura_ondas(altura_ondas==0)=0.1;
velocidade_vento(velocidade_vento==0)=0.1;






%velocidade_vento(velocidade_vento == 99) = NaN;
velocidade_vento = velocidade_vento(velocidade_vento ~= 99);

tempo_ondas = 1:length(altura_ondas);
tempo_vento = 1:length(velocidade_vento);



disp('                                                                                        ==== Item 2.6 ====')
plot(tempo_vento, velocidade_vento);

titulo_vento = sprintf('Velocidade vento - Valores válidos: %d do total (%d)', sum(~isnan(velocidade_vento)), tamanho_vento_original);

title(titulo_vento);
xlabel('t [10 min]');
ylabel('Velocidade [m/s]');
legend('Velocidade do vento');



figure;
% Filtrar os valores não-99
indices_validos = altura_ondas ~= 99;
tempo_valido = tempo_ondas(indices_validos);
altura_ondas_valida = altura_ondas(indices_validos);

plot(tempo_valido, altura_ondas_valida, '.','MarkerSize', 0.2);

titulo_ondas = sprintf('Altura - Valores válidos: %d do total (%d)', length(altura_ondas_valida), tamanho_ondas_original);

title(titulo_ondas);
xlabel('t [10 min]');
ylabel('Altura [m]');
legend('Altura das ondas');


disp('                                                                                        ==== Item 2.7 ====')

figure()
hold on
num_barras_velocidade = ceil(1 + 3.322 * log(length(velocidade_vento)));

[y_hist x_hist] = hist(velocidade_vento,num_barras_velocidade);
hist(velocidade_vento,num_barras_velocidade);
mu_vento = sum(altura_ondas_valida) / length(altura_ondas_valida);
std_vento = std(altura_ondas_valida);


x = altura_ondas_valida(1):0.01:altura_ondas_valida(end);
fprintf('%.2f\n',x);
figure()
pdf = normpdf(x, mu_vento, std_vento);
plot(x, pdf);

xlabel('Velocidade do Vento');
ylabel('PDF');
title('Distribuição Normal da Velocidade do Vento');
hold off;




[frequencias, bordas] = hist(amostras, num_barras);  % Calcular o histograma

dist_barras = bordas(2) - bordas(1);  % Largura da barra
[y_normalizado x_normalizado] = hist(amostras, num_barras, 1/dist_barras);
hist(amostras, num_barras, 1/dist_barras);
xlabel('Valores');
ylabel('Frequência Normalizada');
title('Histograma Normalizado');

disp('                                                                                        ==== Item 2.8 ====')

x = -7:0.01:3;
plot(x, normpdf(x, mu_estimado, sigma_estimado));
y_previsto = normpdf(-7:1:3, mu_estimado, sigma_estimado);
RSS = sumsq(y_previsto - y_normalizado)
fprintf('RSS: %.4f\n',RSS);



k_vento = ceil(1 + 3.322 * log(length(velocidade_vento)));

figure()
[n x] = hist(velocidade_vento, k_vento);
hist(velocidade_vento, k_vento);


k_onda = ceil(1 + 3.322 * log(length(altura_ondas_valida)));
figure()
[n x] = hist(altura_ondas_valida, k_onda);
hist(altura_ondas_valida, k_onda);

























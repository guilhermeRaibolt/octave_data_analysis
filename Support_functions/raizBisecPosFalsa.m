%%=========================================================%%
%% Determinação de raiz de função: Bissecção e Posição Falsa
%% R.-L. p. 41 e p. 47
%% Autor: Thomas W. Rauber trauber@gmail.com , 2008 a 2021
%%
%% Input:
%%		1)	Variável lógica 'posfalsa'. Se for verdadeira,
%%			use Posição Falsa, caso contrário, use Bissecção
%%		2)	Função de qual se quer obter (um)a raiz, um argumento f(x)
%%		3)	Intervalo inicial [a,b]
%%		4)	Precisão de aproximação. Se for atingida, pare a execução
%%		5)	Número máximo de iterações
%%
%% Output:
%%		1)	raiz
%%
%%=========================================================

function r = raizBisecPosFalsa( posfalsa, func, a, b, eps, maxiter )
	if posfalsa
	else
	end
	fa = func(a);
	fb = func(b);
	trocaSinal = fa * fb < 0.0;
	if trocaSinal
	else
  r = NaN;
		return;
	end

	k = 0; d = eps + 1; fx = fa; x = a;
	if ~posfalsa
		if eps <= 0.0
		else
			kmin = (log2(b-a) - log2(eps)) / log2(2);
		end
	end

	while k < maxiter && d > eps && abs(fx) > eps
		if ~posfalsa
			x = (a+b) / 2.0;			%%% BISSECÇÂO
		else
			x = (a*fb - b*fa) / (fb - fa);		%%% POSIÇÂO FALSA
		end
		fx = func(x);
		d = b-a;
		if fa * fx < 0.0  % Troca de sinal acontece no subintervalo esquerdo
			b = x;
			fb = fx;
		else				% fx * fb < 0.0
			a = x;
			fa = fx;
		end
		k = k + 1;
	end
	r = x;
end



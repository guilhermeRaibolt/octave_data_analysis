%%
%% Integração numérica usando a regra 3/8 de Simpson repetida
%%
%% Input: a função a ser integrada, os limites inferior 'a' e superior 'b' de integração
%%        o número 'n' de subdivisões (=número de nós menos 1).
%%
%% Output: Integral numérica
%% Restrições: número 'n' de subdivisões tem que ser divisível por três
%%
function ISR38 = integralSimpson38RepetidaFunc( func, a, b, n, verbose  )

	%% Decide se número é múltiplo de três: Uso de uma função anônima
	% http://www.mathworks.com/help/matlab/matlab_prog/anonymous-functions.html
	ismult3 = @(x) (x - 3*floor(x/3) == 0);

	if ~ismult3( n )
		fprintf('Numero de subdivisoes n=%d. Tem que ser numero multiplo de três !\n', n); % wait();
		ISR38 = NaN;
		return;
	end
	n3 = n/3;
	h = (b-a)/n;
	nos = a:h:b;
	if verbose
		fprintf('Integracao pela Regra 3/8 de Simpson ');
		if n > 4 fprintf('repetida '), end;
		fprintf('no intervalo [%.2f,%.2f] com %d nos\n', a, b, n+1 );
		fprintf('a=%.2f, b=%.2f, n=%d, ===> h=%.2f\n', a, b, n, h);
		printTabXY( nos, 'Nos de interpolacao', func(nos), 'f(nos)', '%9.4f', 'x, f(x)' );
		xidx = 0:1:n;
		midx = ones(n+1, 1);
	end

	limits = func(a) + func(b);
	x = a + h;
	sum2 = 0.0;
	sum3 = 0.0;
	for i=1:n-1
		if ~ismult3(i)
			sum3 = sum3 + func(x); midx(i+1) = 3;
		else
			sum2 = sum2 + func(x); midx(i+1) = 2;
		end
		x = x + h;
	end
	ISR38 = 3*h/8 * (limits + 2 * sum2 + 3 * sum3);

    if ~verbose
      return;
    end
	printTabXY( xidx, 'Indice do no', midx, 'Multiplicador', '%3d', 'Multiplicadores dos nos' );
end


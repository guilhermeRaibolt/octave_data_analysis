function leg = plot_param(f,x,param,titulo,ymax,label_y)
  clf;
  hold on;
  leg = {};
  xlabel( 'x' );
  ylabel(label_y);
  title(titulo);
  axis([0 2.5 0 ymax]);

  for i=1:length(param)
    lamb = param{i}{1};
    k = param{i}{2};
    y = f(x, lamb, k);
    plot(x, y);
    leg{end+1} = sprintf('%s=%3.1f, k=%3.1f', '\lambda' , lamb, k);
  end
  hold off;
  shg; 
endfunction

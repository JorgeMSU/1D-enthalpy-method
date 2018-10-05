% This code gives the analytical solution for the ABT and Shoreline
% trajectories by solving the two nonlinear algebraic equations described
% in Lorenzo-Trueba and Voller 2013.
function F = NonlineSys(y,Rab,a)

F = [(Rab*exp(-(y(2).^2))-(a/2)*(pi.^(1/2))*(erf(y(1))+erf(y(2))))/(exp(-(y(2).^2))+(pi.^(1/2))*y(2)*(erf(y(1))+erf(y(2)))) - 2*(y(2)+a/2)*y(2);
    ((y(2)*Rab+a/2)*exp(-(y(1).^2)))/(exp(-(y(2).^2))+(pi.^(1/2))*y(2)*(erf(y(1))+erf(y(2)))) - y(1)*(1-Rab);];

end



%y(1) = lambda_ab 
%y(2) = lambda_sh
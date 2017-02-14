for i = 2:5
    
    
    % Cell size
    h = (10)^(-i);
      
    Xc = [(h/2):h:(1-(h/2))]';
    
    Xn = [h:h:1-h]';
    
    sigma = 1+(Xc).^2;
    
    pexact = sin(2*pi*Xn);
    
    %mu = 0.1;
    mu = 10;
    
    q = (2*pi*(cos(2*pi*Xn).*(2*Xn+mu)-2*pi*sin(2*pi*Xn).*(1+Xn.^2)));
    
    pa = convAdvec(h,sigma,mu,q);
    
    %error
    e = norm(pexact-pa,inf);
    
    fprintf('%3.2e\n',e);
    
    
    plot(Xn,pexact,Xn,pa);
    
    
    
end

% Discussion
% The plot that is produced shows close to identical graphs for both values of mu. 
%
%My error for mu = 0.1 is 
%3.60e-04
%3.60e-06
%3.61e-08
%3.97e-08
%
% Error for mu = 10 is
%7.97e-04
%7.97e-06
%7.97e-08
%2.17e-08
% These error results make sense, because they decrease for every value of
% i. i runs from 2:5 and effects the h = 10^(-i) so as i gets larger the
% error in h should decrease as it does. 
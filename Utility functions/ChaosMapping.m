function ChaosValue = ChaosMapping(value,method)
%

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
if nargin < 2
    method = 'Gaussian';
end
ChaosValue = feval(method,value);
end

function ChaosValue=Chebyshev(value)
    chebyshev=4;%混沌系数
    ChaosValue=cos(chebyshev.*acos(value));
end

function ChaosValue=Logistic(value)
%0<=mu<=4
    mu=4;
    ChaosValue=mu.*value.*(1-value);
end

function ChaosValue=Tent(value)
%0<value<1
beta=0.6;
ChaosValue=zeros(size(value,1),size(value,2));
ChaosValue(value<beta)=value(value<beta)./beta;
ChaosValue(value>=beta)=(1-value(value>=beta))./(1-beta);

end

function ChaosValue=Cubic(value)
    rho=2.59;
    ChaosValue=rho*value*(1-value.^2);
end

function ChaosValue=Bernoulli(value)
%0<value<1
    lambda=0.4;
    if value<=(1-lambda) && value>0
        ChaosValue=value/(1-lambda);
    else
        ChaosValue=(value-1+lambda)/lambda;
    end
end

function ChaosValue=PWLCM(value)
%0<value<1
    p=0.7;
    if value<p && value>0
        ChaosValue=value/p;
    else
        ChaosValue=(1-value)/(1-p);
    end
end

function ChaosValue=Singer(value)
%0.9<=mu<1.08
    mu=1.073;
    ChaosValue=mu*(7.86*value-23.31*value.^2+28.75*value.^3-13.302875*value.^4);
end

function ChaosValue=Sine(value)
    a=4;
    ChaosValue=a/4.*sin(pi.*value);
end

function ChaosValue=Gaussian(value)
    mu=1;
    if value==0
        ChaosValue=0;
    else
        ChaosValue=rem(mu./value,1);
    end
end

function ChaosValue=Circle(value)
    a=0.5;
    b=2.2;
    ChaosValue=value+a-mod(b/(2*pi).*(sin(2*pi.*value)),1);
end

function ChaosValue=Sinusoidal(value)
    a=2.3;
    ChaosValue=a.*(value.^2).*sin(pi.*value);
end

function ChaosValue=ICMIC(value)
    a=70;
    ChaosValue=sin(a./value);
end


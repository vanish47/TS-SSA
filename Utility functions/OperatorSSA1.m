function Offspring = OperatorSSA1(Problem,Parent,iterator,r2,Parameter)
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

if nargin > 4
    [ST,SD] = deal(Parameter{:});%ST（0.5，1）:警戒值；SD:警戒者比例，一般不设置太高
else
    [ST,SD] = deal(0.6,0.1);
end
%通过非支配排序，对种群进行分层，处于F0的为生产者集群，即最优集群，处于Fmax的为最差集群。
[FrontNo,MaxFNo] = NDSort(Parent.objs,Parent.cons,Problem.N);
N=Problem.N;
D=Problem.D;
M=Problem.M;
lower=Problem.lower;
upper=Problem.upper;
if isa(Parent(1),'SOLUTION')
    evaluated  = true;
    ParentDecs = Parent.decs;
    ParentDecs(:,D+1)=FrontNo;
    ParentDecs=sortrows(ParentDecs,D+1);
    ParentObjs=Parent.objs;
    ParentObjs(:,M+1)=FrontNo;
    ParentObjs=sortrows(ParentObjs,M+1);
    best=sum(ParentDecs(:,D+1)==1);
    worst=sum(ParentDecs(:,D+1)==MaxFNo);
else
    evaluated = false;
end

Offspring = zeros(size(Parent,1),size(Parent,2));
%对种群进行生产者和跟随者的位置更新。
A=randi(2,1,D);
A(A==2)=-1;
A=A'*inv(A*A');
A=A';

iterMax=Problem.maxFE/N;

% a=rand;
% l=exp(3*cos(pi*(iterator-iterMax)^2/(iterMax^2)));
% beta=exp(a*l)*cos(2*pi*a);

% c=0.9;
% omega=c^(iterator^2);

weight=1-((exp(iterator/iterMax)-1)/exp(1)-1)^2;%惯性权重
ST=(1+2*0.5)/3+(1-0.5)/3*tanh(-5+10*iterator/iterMax);%前期比较危险，安全阈值低；后期安全，阈值提高
if best==N
    best=floor(0.4*N);
    worst=N-best;
end
if r2<ST

    % % 二进制交叉
    % [proC,disC] = deal(1,20);
    % Parent1   = ParentDecs(1:floor(best/2),1:D);
    % Parent2   = ParentDecs(floor(best/2)+1:floor(best/2)*2,1:D);
    % OffspringBest = zeros(2*size(Parent1,1),size(Parent1,2));
    % [B,D] = size(Parent1);
    % beta  = zeros(B,D);
    % mu    = rand(B,D);
    % beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    % beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    % beta = beta.*(-1).^randi([0,1],B,D);
    % beta(rand(B,D)<0.5) = 1;
    % beta(repmat(rand(B,1)>proC,1,D)) = 1;
    % OffspringBest = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
    %     (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
    % ParentDecs(1:2*B,1:D)=OffspringBest;
    %%常规麻雀算法搜索
    % for i=1:best
    %     ParentDecs(i,1:D)=ParentDecs(i,1:D)*exp(-i/(iterMax.*rand))*beta;
    % end

    % %差分变异
    % F=0.5;
    % %rand
    % for i=1:best
    %     ParentDecs(i,1:D)=ParentDecs(i,1:D)+randn*(ParentDecs(randi(best),1:D)-ParentDecs(randi(best),1:D));
    % end
    %柯西变异
    for i=1:best
        ParentDecs(i,1:D)=ParentDecs(i,1:D)*exp(-i/(iterMax.*rand))*(1+tan((rand-1/2)*pi));
    end
else
    for i=1:best
        ParentDecs(i,1:D)=ParentDecs(i,1:D)+(1+tan((rand-1/2)*pi))*(upper-lower);
    end


end
for i=best+1:N
    if i>N/2
        % %余弦相似度
        % Cosine   = 1 - pdist2(ParentObjs(i,1:D),ParentObjs(N-worst+1:N,1:D),'cosine');%obj为距离
        % Distance = repmat(sqrt(sum(ParentObjs(i,1:D).^2,2)),1,best).*sqrt(1-Cosine.^2);
        % [~,w] = min(Distance',[],1);
        % ParentDecs(i,1:D)=weight*randn*exp((ParentDecs(w,1:D)-ParentDecs(i,1:D))/(i*i));

        %常规方法
        ParentDecs(i,1:D)=weight*randn*exp((ParentDecs(randi([N-worst+1,N]),1:D)-ParentDecs(i,1:D))/(i*i));%有惯性权重
        % ParentDecs(i,1:D)=randn*exp((ParentDecs(randi([N-worst+1,N]),1:D)-ParentDecs(i,1:D))/(i*i));%无惯性权重
    else
        % %余弦相似度
        % Cosine   = 1 - pdist2(ParentObjs(i,1:M),ParentObjs(1:best,1:M),'cosine');%obj为距离
        % Distance = repmat(sqrt(sum(ParentObjs(i,1:M).^2,2)),1,best).*sqrt(1-Cosine.^2);
        % [~,b] = min(Distance',[],1);
        % ParentDecs(i,1:D)=ParentDecs(b,1:D)+abs(ParentDecs(i,1:D)-ParentDecs(b,1:D)).*A*weight;%有惯性权重

        %常规方法
        ParentDecs(i,1:D)=ParentDecs(randi([1,best]),1:D)+abs(ParentDecs(i,1:D)-ParentDecs(randi([1,best]),1:D)).*A*weight;%有惯性权重
        % ParentDecs(i,1:D)=ParentDecs(randi([1,best]),1:D)+abs(ParentDecs(i,1:D)-ParentDecs(randi([1,best]),1:D)).*A;%无惯性权重
    end
end

SD=fix(N*SD);

for i=1:SD
    index=randi(N);
    if index<best
        % %余弦相似度
        % Cosine   = 1 - pdist2(ParentObjs(i,1:M),ParentObjs(N-worst+1:N,1:M),'cosine');%obj为距离
        % Distance = repmat(sqrt(sum(ParentObjs(i,1:M).^2,2)),1,worst).*sqrt(1-Cosine.^2);
        % [~,w] = min(Distance',[],1);
        % f=ParentObjs(index,1:M)-ParentObjs(w,1:M);
        % ParentDecs(index,1:D)=ParentDecs(index,1:D)+(rand*2-1)*abs(ParentDecs(index,1:D)-ParentDecs(w,1:D))*sum(1./f)/M;

        %常规方法
        w=randi([N-worst+1,N]);
        f=ParentObjs(index,1:M)-ParentObjs(w,1:M);
        ParentDecs(index,1:D)=ParentDecs(index,1:D)+(rand*2-1)*abs(ParentDecs(index,1:D)-ParentDecs(w,1:D))*sum(1./f)/M;
        %ParentDecs(index,1:D)=upper-lower-ParentDecs(index,1:D);%反向学习
    else
        % %余弦相似度
        % Cosine   = 1 - pdist2(ParentObjs(i,1:M),ParentObjs(1:best,1:M),'cosine');%obj为距离
        % Distance = repmat(sqrt(sum(ParentObjs(i,1:M).^2,2)),1,best).*sqrt(1-Cosine.^2);
        % [~,b] = min(Distance',[],1);
        % ParentDecs(index,1:D)=ParentDecs(b,1:D)+rand*abs(ParentDecs(index,1:D)-ParentDecs(b,1:D));
        %常规方法
        b=randi([1,best]);
        ParentDecs(index,1:D)=ParentDecs(b,1:D)+rand*abs(ParentDecs(index,1:D)-ParentDecs(b,1:D));
    end
end


Offspring=ParentDecs(:,1:D);
%反向种群
%Offspring=[Offspring;upper+lower-Offspring];
if evaluated
    Offspring = Problem.Evaluation(Offspring);
end
end





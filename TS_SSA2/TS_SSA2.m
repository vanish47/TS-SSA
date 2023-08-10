classdef TS_SSA2 < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained/none>
% Many objects sparrow algorithm

%------------------------------- Reference --------------------------------
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            [Z,Problem.N]        = UniformPoint(Problem.N,Problem.M);
            Population           = Problem.chaosInitialization();
            Zmin                 = min(Population(all(Population.cons<=0,2)).objs,[],1);
            iterator=0;
            r2=0.6;
            %% Optimization
            while Algorithm.NotTerminated(Population)
                iterator=iterator+1;                
                r2=ChaosMapping(r2,'Sine');
                MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));
                Offspring  = OperatorSSA2(Problem,Population(MatingPool),iterator,r2);
                Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin);
            end

            % Archive    = UpdateArchive(Population(NDSort(Population.objs,1)==1),[],Problem.N);
            % Archive = Population;
            % while Algorithm.NotTerminated(Archive)
            %     iterator=iterator+1;
            %     r2=ChaosMapping(r2,'Logistic');
            %     MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));
            %     Offspring  = OperatorSSA(Problem,Population(MatingPool),iterator,r2);
            %     Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
            %     Population = EnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin);
            %     Archive    = UpdateArchive(Archive,Population,Problem.N);
            % end
        end
    end
end
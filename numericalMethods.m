classdef numericalMethods
    
    properties (Access = public)
       solExact
       y1Exact
       y2Exact
       
       solExplicit
       y1Explicit
       y2Explicit
       
       solImplicit
       y1Implicit
       y2Implicit
       
       solCrankNic
       y1CrankNic
       y2CrankNic
    end % end of properties
    
    methods
        function soln = numericalMethods()
            
        end
        function soln = computeExact(soln, prob)
           soln.solExact = zeros(length(prob.t),2);
           soln.solExact(1,:) = prob.yexact(prob.t(1));
           for i = 1 : prob.nt+1
               soln.solExact(i,:) = (prob.yexact(prob.t(i)))';
           end
           soln.y1Exact = soln.solExact(:,1);
           soln.y2Exact = soln.solExact(:,2);
        end % end of computeExact
        
        function soln = computeExplicit(soln, prob)
            soln.solExplicit = zeros(length(prob.t),2);
            soln.solExplicit(1,:) = prob.y0';
            for i = 2 : prob.nt +1
               soln.solExplicit(i,:) = ((eye(2) + prob.h * prob.matrix)*...
                   soln.solExplicit(i-1,:)'); 
            end
            soln.y1Explicit = soln.solExplicit(:,1);
            soln.y2Explicit = soln.solExplicit(:,2);

        end % end  of computeExplicit
        
        function soln = computeImplicit(soln, prob)
            soln.solImplicit = zeros(length(prob.t),2);
            soln.solImplicit(1,:) = prob.y0';
            for i = 2 : prob.nt +1
               soln.solImplicit(i,:) = inv(eye(2)- prob.h * prob.matrix)...
                           *((soln.solImplicit(i-1,:))'); 
            end
            soln.y1Implicit = soln.solImplicit(:,1);
            soln.y2Implicit = soln.solImplicit(:,2);
        end
        
        function soln = computeCrankNicolson(soln ,prob)
            soln.solCrankNic = zeros(length(prob.t),2);
            soln.solCrankNic(1,:) = prob.y0';
            for i = 2: prob.nt + 1
                soln.solCrankNic(i,:) = (inv(eye(2) - 0.5 * prob.h * prob.matrix)...
                           *(eye(2) + 0.5 * prob.h * prob.matrix) * (soln.solCrankNic(i-1,:))'); 
            end
            soln.y1CrankNic = soln.solCrankNic(:,1);
            soln.y2CrankNic = soln.solCrankNic(:,2);
        end
        function soln = solve(soln,prob)
           % Exact
           soln = soln.computeExact(prob); 
           % Explicit
           soln = soln.computeExplicit(prob);
           % Implicit
           soln = soln.computeImplicit(prob);
           % Crank Nicolson
           soln = soln.computeCrankNicolson(prob);
           
           % plotting y1 vs t
           figure(1)
           plot(prob.t ,soln.y1Exact,...
               prob.t ,soln.y1Explicit,...
               prob.t ,soln.y1Implicit,...
               prob.t ,soln.y1CrankNic)
           xlabel('t'); ylabel('y1');
           legend('Exact','Explicit Euler',...
               'Implicit','Crank-Nicolson');
           % plotting y2 vs t
           figure(2)
           plot(prob.t ,soln.y2Exact,...
               prob.t ,soln.y2Explicit,...
               prob.t ,soln.y2Implicit,...
               prob.t ,soln.y2CrankNic)
           xlabel('t'); ylabel('y2');
           legend('Exact','Explicit Euler',...
               'Implicit','Crank-Nicolson');
           % plotting y1 vs y2
           figure(3)
           plot(soln.y1Exact ,soln.y2Exact,...
               soln.y1Explicit ,soln.y2Explicit,...
               soln.y1Implicit,soln.y2Implicit,...
               soln.y1CrankNic ,soln.y2CrankNic)
           xlabel('y1'); ylabel('y2');
           legend('Exact','Explicit Euler',...
               'Implicit','Crank-Nicolson');
           % plotting 3d y1 vs y2 vs t
           figure(4)
           plot3(soln.y1Exact ,soln.y2Exact,prob.t,...
               soln.y1Explicit ,soln.y2Explicit,prob.t,...
               soln.y1Implicit ,soln.y2Implicit,prob.t,...
               soln.y1CrankNic ,soln.y2CrankNic,prob.t)
           xlabel('y1'); ylabel('y2'); zlabel('t');
           legend('Exact','Explicit Euler',...
               'Implicit','Crank-Nicolson');
        
        end
    end % end of methods
    
end % end of classdef
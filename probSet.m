classdef probSet
   properties ( Access = public )
       f;
       yexact;
       y0;
       matrix;
       
       t0;
       tf;
       nt;
       h;
       t;
       
   end
   
   methods
       function prob = probSet(type)
           switch type
               case 1
                    prob.f = @(t)  [-0.9 -6.3; 6.3 -0.9 ] .* y(t);
                    prob.y0 = [ 1; 0];
                    prob.yexact = @(t) exp(-0.9 * t) * [cos(6.3 * t);sin(6.3 * t)] ; 
                    prob.matrix = [-0.9 -6.3; 6.3 -0.9 ];
                    
                    % time
                    prob.t0 = 0;
                    prob.tf = 10;
                    prob.nt = 1000;
                    prob.h = (prob.tf - prob.t0) /prob.nt;
                    prob.t = [prob.t0 : prob.h : prob.tf];       

                  
           end
       end
   end
end
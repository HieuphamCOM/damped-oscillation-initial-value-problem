%% Solve the damped oscillation problem
% Written by : Pham Nguyen Minh Hieu
% Write Matlab to solve the problem by using:
% 1/ Explicit Euler Method
% 2/ Implicit Euler Method
% 3/ Cank-Nicolson

clear all
close all
clc 

type = 1;
prob = probSet(type);

soln = numericalMethods();

soln = soln.solve(prob);
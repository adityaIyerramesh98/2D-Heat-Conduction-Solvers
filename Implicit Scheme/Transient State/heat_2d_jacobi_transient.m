%Solving 2D Heat Conduction Equation by Gauss Jacobi Method
%(Trasient State Implicit Scheme)
%Author: Aditya Iyer Ramesh

close all
clc

%Initializing parameters
%grid points
nx = 10;
ny = nx;
nt = 1400;
Lx = 1; 
Ly = 1;

%For Mesh
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);

%Spacing of grids
dx = x(2) - x(1);
dy = dx;

%Error criteria & tolerance (EC > Tol)
error = 9e9;
tol = 1e-4;
dt = 1e-3;

%Thermal Diffusivity
alpha = 1.4;

%Temperature Matrix
Temp = 300*ones(nx,ny);

%Boundary Conditions
Temp(nx,2:ny-1) = 600; %Temp at Top
Temp(1,2:ny-1) = 900;   %Temp at Bottom
Temp(2:nx-1,1) = 400;   %Temp at Left
Temp(2:nx-1,ny) = 800;   %Temp at Right


%For Average Temperatures at Corners
Temp(nx,1) = (600 + 400)/2;
Temp(1,1) = (900 + 400)/2;
Temp(1,ny) = (900 + 800)/2;
Temp(nx,ny) = (600 + 800)/2;

%Assigning Meshes
[xx,yy] = meshgrid(x,y);

%Solver for Jacobi Iteration
iterative_solver = 1;

%Counter for iteration 
iteration_counter = 1;

%Declaring original values of Temp
Temp_old = Temp;
Temp_initial_dt = Temp;

%Computing of 2D Transient Heat Conduction Eqn by Jacobi Solver starts
tic;

if iterative_solver == 1
    dt = 1e-3;
    k1 = (alpha*dt)/(dy^2);
    
    for nt = 1:1400
        error = 9e9;
        
        while (error > tol)
            for i = 2:nx-1
                for j = 2:ny-1
                    Temp_orig = Temp_old(i-1,j) + Temp_old(i+1,j) + Temp_old(i,j-1) + Temp_old(i,j+1);
                    Temp(i,j) = Temp_initial_dt(i,j)/(1+(4*k1)) + ((k1)/(1+(4*k1)))*Temp_orig;
                end
            end
            error = max(max(abs(Temp_old - Temp)));
           
            %Updating value of temperatures after each iteration
            Temp_old = Temp;
            Temp_initial_dt = Temp;
            iteration_counter = iteration_counter + 1;
        end
    end
    
    toc;
    time = toc;

%Contour Plots
figure(1)
contourf(x,y,Temp)
clabel(contourf(x,y,Temp))
colorbar
colormap('jet')

%Labelling
xlabel("X-Axis")
ylabel("Y-Axis")
title(sprintf('No. of Trans.St Jacobi Iterations for IMPLICIT SCH = %d, Total Time = %f seconds',iteration_counter,time));
end    

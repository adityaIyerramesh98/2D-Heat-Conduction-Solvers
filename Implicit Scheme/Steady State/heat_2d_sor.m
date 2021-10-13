%Solving 2D Heat Conduction Equation by Sucessive Relaxation Method
%Implicit Scheme - (Steady State)
%Author: Aditya Iyer Ramesh

close all
clc

%Initializing parameters
%grid points
nx = 10;
ny = nx;

x = linspace(0,1,nx);
y = linspace(0,1,ny);

dx = x(2) - x(1);
dy = dx;

%Error criteria & tolerance (EC > Tol)
error = 9e9;
tol = 1e-4;
omega = 1.2;

%Boundary Conditions
Temp_Left = 400; 
Temp_Top = 600;
Temp_Right = 800;
Temp_Bottom = 900;

Temp = 300*ones(nx,ny);

Temp(2:ny - 1, 1) = Temp_Left;
Temp(2:ny - 1, nx) = Temp_Right;
Temp(1, 2:nx - 1) = Temp_Top;
Temp(ny, 2:nx - 1) = Temp_Bottom;

%For Average Temperatures at Corners
Temp(1,1) = (Temp_Top + Temp_Left)/2;
Temp(nx,ny) = (Temp_Right + Temp_Bottom)/2;
Temp(1,ny) = (Temp_Top + Temp_Right)/2;
Temp(nx,1) = (Temp_Left + Temp_Bottom)/2;

%Declaring original values of Temp
Temp_initial = Temp;
Temp_old = Temp;

%Computing of 2D SS Heat Conduction Eqn SOR Solver starts
iterative_solver = 3;

if iterative_solver == 3
    sor_iteration = 1;
    
    while (error > tol)
        
        for i = 2:nx - 1
            
            for j = 2:ny - 1
                Temp(i,j) = omega*(0.25*(Temp(i-1,j) + Temp_old(i+1,j) + Temp(i,j-1) + Temp_old(i,j+1))) + (1 - omega)*Temp_old(i,j);
            end
        end
        
        error = max(max(abs(Temp_old - Temp)));
        Temp_old = Temp;
        sor_iteration = sor_iteration + 1;
        
    end
end

%Contour Plots
figure(1)
contourf(x,y,Temp)
clabel(contourf(x,y,Temp))
colorbar
colormap(jet)
set(gca,'ydir','reverse')

%Labelling
xlabel("X-Axis")
ylabel("Y-Axis")
title(sprintf('No. of SS Sucessive Over Relaxation Iterations for IMPLICIT SCH = %d',sor_iteration));
pause(0.03);